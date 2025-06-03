package edu.ucsc.neurobiology.vision.anf;

import java.io.*;
import java.text.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.analysis.*;
import edu.ucsc.neurobiology.vision.calculations.*;
import edu.ucsc.neurobiology.vision.electrodemap.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.plot.*;
import static edu.ucsc.neurobiology.vision.util.VisionParams.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * Calculation that employs Expectation Maximization (EM) to cluster the spikes on every
 * electrode into clusters. This calculation saves the results in the .neurons-raw file.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 * @owner Matthew Grivich, University of California, Santa Cruz
 */
public class PCANFClustering extends AbstractCalculation {

    private String datasetFolder, reportDirectory;
    private ElectrodeMap electrodeMap;

    private int fromElectrode, toElectrode;
    private boolean generateReport;
    private int binsPerDimension;
    private double clusteringSignificance;
    private int minClusters;
    private int maxClusters;
    private int nEMSpikes;
    private int nClusteringThreads;
    private int densityClusteringMaxSpikeLoss;
    private int minEMIterations;
    private int maxEMIterations;
    private double emLikelihoodDelta;
    private Vision app;
    private DecimalFormat format = new DecimalFormat();
    private ArrayList<Integer> electrodePool = new ArrayList<Integer>();
    private ClusteringModelFile modelFile;
    private ProjectionsFile projectionsFile;
    private NeuronFile neuronFile;
    private CovarianceFile covarianceFile;
    private VisionHeader prjFileHeader;


    private synchronized int getNextElectrodeToBeAnalyzed() {
        if (electrodePool.isEmpty()) {
            return -1;
        } else {
            int N = toElectrode - fromElectrode + 1;
            int n = electrodePool.size();
            app.setProgress(100 * (N - n) / N);

            return electrodePool.remove(0);
        }
    }


    class ClusteringThread extends Thread {

        private PlotPanel pInitial, pFinal;
        private double[] contaminationIndex, nBadSpikes;
        private DoubleHistogram[] autocorrelations;
        private int[] electrodeSpikeTimes;
        private IntegerList[] neuronSpikeTimes;
        private int nSpikesOnElectrode;
        private float[][] projections;
        private ExpectationMaximization expectationMaximization;
        PrintWriter html = null;

        
        public ClusteringThread() {
            electrodeSpikeTimes = new int[projectionsFile.maxSpikesPerElectrode];
            projections = new float[prjFileHeader.nDimensions][projectionsFile.maxSpikesPerElectrode];
            
            expectationMaximization = new ExpectationMaximization1(
                    prjFileHeader.nDimensions, maxClusters,
                    projectionsFile.maxSpikesPerElectrode);

            neuronSpikeTimes = new IntegerList[maxClusters];
            for (int i = 0; i < neuronSpikeTimes.length; i++) {
                neuronSpikeTimes[i] = new IntegerList(100000, 2);
            }
        }
        
        
        public void run() {
            while (true) {
                int electrode = getNextElectrodeToBeAnalyzed();
                if (electrode == -1) {
                    return; // no more electrodes to cluster
                } else {
                    try {
                        if (clusterElectrode(electrode)) {
                            try {
                            saveNeurons(electrode);
                            } catch(NullPointerException ex) {
                                ex.printStackTrace();
                                System.out.println("electrode: " + electrode);
                            }
                        } else {
                            if (generateReport) {
                                if(!Vision.isGUIBased()) {
                                    throw new IllegalStateException("HTML report generation requires the Vision GUI to be active.");
                                }
                                html = initHTML(electrode);
                                closeHTML();
                            }

                        }
                    } catch (IOException ex) {
                        ex.printStackTrace();
                    }
                }
            }
        }


        private boolean clusterElectrode(int electrode) throws IOException {
            // load the data for this electode
            //			app.sendMessage("PCA on " + electrode + " : Loading projections...");

            nSpikesOnElectrode = projectionsFile.readProjections(electrode, projections,
                    electrodeSpikeTimes);

            if (nSpikesOnElectrode <= 0) {
                return false;
            }

            // density clustering
            ArrayList<double[]> means = new ArrayList<double[]>();
            ArrayList<double[]> sigmas = new ArrayList<double[]>();

            double[] min = new double[prjFileHeader.nDimensions];
            double[] max = new double[prjFileHeader.nDimensions];
            double[] binInterval = new double[prjFileHeader.nDimensions];
            densityClustering(
                    projections, nSpikesOnElectrode, binsPerDimension,
                    densityClusteringMaxSpikeLoss, clusteringSignificance, means, sigmas, min,
                    max, binInterval, true, minClusters, maxClusters);
            /*
                        // find the minimum and maximum limits on all the clustering dimensions
                        for (int d = 0; d < projectionsFile.nDimensions; d++) {
                            min[d] = Double.POSITIVE_INFINITY;
                            max[d] = Double.NEGATIVE_INFINITY;

                            // find the maximum and minimum
                            for (int i = 0; i < nSpikesOnElectrode; i++) {
                                if (projections[d][i] < min[d]) {
                                    min[d] = projections[d][i];
                                }
                                if (projections[d][i] > max[d]) {
                                    max[d] = projections[d][i];
                                }
                            }

                            binInterval[d] = (max[d] - min[d]) / binsPerDimension;
                        }

                        // fill the gaussians
                        while (means.size() < maxClusters) {
                            double[] mean = new double[projectionsFile.nDimensions];
                            double[] sigma = new double[projectionsFile.nDimensions];
                            for (int d = 0; d < projectionsFile.nDimensions; d++) {
//                    mean[d] = min[d] + (max[d] - min[d]) * random.nextDouble();
                                mean[d] = min[d] + (max[d] - min[d]) * 0.5;

                                sigma[d] = binInterval[d] * binInterval[d];
                            }
                            means.add(mean);
                            sigmas.add(sigma);
                        }
             */

            // do the EM fit
            expectationMaximization.reset(projections, nEMSpikes, nSpikesOnElectrode);
            for (int i = 0; i < means.size(); i++) {
                expectationMaximization.addGaussian(
                        (double[]) means.get(i), (double[]) sigmas.get(i));
            }

            if (generateReport) {
                pInitial = expectationMaximization.getInitialConditions();
            }

            try {
                int nEMIterations = expectationMaximization.fit(
                        emLikelihoodDelta, minEMIterations, maxEMIterations);
            } catch (FitFailedException e) {
                return false;
            }
            updateNeuronsInfo();

            if (generateReport) {
                pFinal = expectationMaximization.getResults(1);
            }

            return true;
        }


        private void updateNeuronsInfo() {
            final int nGaussians = expectationMaximization.getClustersCount();
            contaminationIndex = new double[nGaussians];
            nBadSpikes = new double[nGaussians];
            autocorrelations = new DoubleHistogram[nGaussians];
            for (int i = 0; i < neuronSpikeTimes.length; i++) {
                neuronSpikeTimes[i].clear();
            }

            // extract the neuron times and compute the autocorrelations
            for (int neuron = 0; neuron < nGaussians; neuron++) {
                for (int i = 0; i < nSpikesOnElectrode; i++) {
                    if (expectationMaximization.getCluster(i) == neuron) {
                        neuronSpikeTimes[neuron].add(electrodeSpikeTimes[i]);
                    }
                }

                autocorrelations[neuron] = AutocorrelationCalculator.calculateISI(neuronSpikeTimes[neuron], 50, 0.05);
                double N = neuronSpikeTimes[neuron].size();
                nBadSpikes[neuron] = autocorrelations[neuron].count(ACFT1, ACFT2);
                double T = neuronFile.getNumberOfSamples() / VisionParams.SAMPLES_PER_MILLISECOND;
                double dT = ACFT2 - ACFT1;
                contaminationIndex[neuron] = nBadSpikes[neuron] * T / (dT * N * N);
            }
        }


        public void saveNeurons(int electrode) throws IOException {
            int nNeurons = expectationMaximization.getClustersCount();

            ClusteringModelFile.EMModel m = new ClusteringModelFile.EMModel();
            m.extractionID = electrode;
            m.neuronIndex = new int[nNeurons];
            m.neuronID = new int[nNeurons];
            if (generateReport) {
                html = initHTML(electrode);
            }
            for (int neuron = 0; neuron < nNeurons; neuron++) {
                // write the neuron to the neuron file

                // Neuron ID is now assigned automatically according the cluster, electrode
                int neuronID = NeuronFile.getNeuronID(electrode, neuron);
                neuronFile.addNeuron(electrode, neuronID,
                        neuronSpikeTimes[neuron],
                        neuronSpikeTimes[neuron].size());

                m.neuronIndex[neuron] = neuron;
                m.neuronID[neuron] = neuronID;

                if (generateReport) {
                    writeNeuron(neuron, neuronID, electrode);
                }
            }
            if (generateReport) {
                closeHTML();
            }
            m.cleaningLevel = 1;
            m.electrodes = electrodeMap.getAdjacentsTo(
                    electrode, prjFileHeader.electrodeUsage.ordinal());
            m.threshold = 0; //FIXME
            m.nDimensions = expectationMaximization.getDimensionsCount();

            m.eigenvectors = new double[prjFileHeader.nDimensions][];
            float[] covarianceMatrix = covarianceFile.getCovarianceMatrix(electrode);
            // Calculate eigenvectors
            try {
                PCA pca = new PCA(covarianceMatrix);
                pca.doPCA();
                for (int d = 0; d < prjFileHeader.nDimensions; d++) {
                    m.eigenvectors[d] = pca.getEigenVector(d);
                }
            } catch (TooManyIterationsException e) {
                e.printStackTrace();
            }

            m.nClusters = expectationMaximization.getClustersCount();

            // specific to EM
            m.nGaussians = expectationMaximization.getGaussiansCount();
            m.probability = new double[m.nGaussians];
            m.means = new double[m.nGaussians][];
            m.covariances = new double[m.nGaussians][];
            for (int i = 0; i < m.nGaussians; i++) {
                m.probability[i] = expectationMaximization.getGaussianProbability(i);
                m.means[i] = expectationMaximization.getMeans(i);
                m.covariances[i] = expectationMaximization.getCovariances(i);
            }

            modelFile.addExtraction(m);
        }


        private PrintWriter initHTML(int electrode) throws IOException {
            // prepare the HTML document
            PrintWriter html = null;
            html = new PrintWriter(new FileWriter(
                    reportDirectory + format.format(electrode) + ".html"), true);
            html.write("<html><body>");
            html.write("EM Dimensions: " + prjFileHeader.nDimensions + ".        ");

            String next = format.format(electrode + 1) + ".html";
            String previous = format.format(electrode - 1) + ".html";

            html.write("<a href=\"" + previous + "\">Previous</a>");
            html.write(" : ");
            html.write("<a href=\"" + next + "\">Next</a>");
            if (nSpikesOnElectrode > -1) {
                html.println("<hr><table><tr>");

                pInitial.addToLegend("Electrode " + electrode);
                pInitial.addToLegend("Spikes: " + nSpikesOnElectrode);
                //				pInitial.addToLegend("Electrode " + electrode);
                pInitial.addToLegend("Components: " +
                        expectationMaximization.getClustersCount());
                //				html.write("<br>PCA: " + VisionUtilities.format(EMPercentage, 1) + "%");
                //				html.write("<br>Threshold: " + VisionUtilities.format(clusteringThreshold, 1));
                //				html.write("<br>MSPB: " + maxSpikesPerBin);
                //				html.write("<br>EM: " + nEMIterations + " iter.");

                String name = format.format(electrode) + "initial.png";
                pInitial.saveAsPNG(new File(reportDirectory + name));
                html.write("<td><img src=" + name + ">");

                name = format.format(electrode) + "result.png";
                pFinal.saveAsPNG(new File(reportDirectory + name));
                html.write("<td><img src=" + name + ">");
            } else {
                html.println("<br> No spikes found on electrode " + electrode + ".<br>");
            }
            return html;
        }


        private void writeNeuron(int neuronIndex, int id, int electrode) throws
        IOException {
            int nSpikes = neuronFile.getSpikeCount(id);

            if (neuronIndex % 2 == 0) {
                html.write("<tr>");
            }

            // save the Autocorrelation
            PlotPanel pAuto = new PlotPanel();
            pAuto.addToLegend("ID: " + id + " : " +
                    PlotUtil.getColorName(neuronIndex));
            pAuto.addToLegend("Spikes: " + nSpikes);
            pAuto.addToLegend("C: " +
                    StringUtil.format(contaminationIndex[neuronIndex], 4));
            pAuto.setAxisVisible(false);
            pAuto.addData(autocorrelations[neuronIndex], new HistogramStyle());
            pAuto.autoscale();
            pAuto.setSize(350, 150);

            String fileName = format.format(electrode) +
            PlotUtil.getColorName(neuronIndex) + "auto.png";
            pAuto.saveAsPNG(new File(reportDirectory + fileName));
            html.write("<td><img src=" + fileName + ">");
        }


        private void closeHTML() {
            html.write("</table></body></html>");
            html.close();
        }
    }


    public void startCalculation() throws Exception {
        format.setMaximumFractionDigits(0);
        format.setMinimumIntegerDigits(3);
        format.setMaximumIntegerDigits(3);
        app = Vision.getInstance();

        String datasetName = new File(datasetFolder).getName();

        reportDirectory = datasetFolder + File.separator + "Clustering Plots" +
        File.separator;
        if (generateReport) {
            new File(reportDirectory).mkdir();
        }

        // open the .prj file
        app.sendMessage("PCANF:Clustering, Open projections file...");
        String prjFileName = datasetFolder + File.separator + datasetName + ".prj";
        String covFileName = datasetFolder + File.separator + datasetName + ".wcov";
        File covF = new File(covFileName);
        //use standard covariance file if whitened covariance file was not calculated.
        if(!covF.exists()) {
            covFileName = datasetFolder + File.separator + datasetName + ".cov";
        }
        projectionsFile = new ProjectionsFile(prjFileName);
        prjFileHeader = projectionsFile.getHeader();
        covarianceFile = new CovarianceFile(covFileName);
        electrodeMap = ElectrodeMapFactory.getElectrodeMap(prjFileHeader.arrayID);

        if (toElectrode == -1) {
            toElectrode = electrodeMap.getNumberOfElectrodes() - 1;
        }

        //		SpikeFile spikes = new SpikeFile(
        //		datasetFolder + File.separator + datasetName + ".spikes");

        // create the neuron file
        app.sendMessage("PCANF:Clustering, Create neuron file...");
        String neuronsFileName = datasetFolder + File.separator + datasetName + ".neurons-raw";

        VisionHeader neuronFileHeader = null;
        File nF = new File(neuronsFileName);
        if (nF.exists() && nF.isFile() && nF.canRead()) {
            //  neuronFile = new NeuronFile(neuronsFileName);
            nF.delete();
        } 

        neuronFileHeader = projectionsFile.getHeader();

        neuronFileHeader.version = NeuronFile.INT_VERSION;
        neuronFileHeader.binsPerDimension = binsPerDimension;
        neuronFileHeader.clusteringSignificance = clusteringSignificance;
        neuronFileHeader.densityClusteringMaxSpikeLoss =
            densityClusteringMaxSpikeLoss;
        neuronFileHeader.minClusters = minClusters;
        neuronFileHeader.maxClusters = maxClusters;
        neuronFileHeader.nEMSpikes = nEMSpikes;
        neuronFileHeader.minEMIterations = minEMIterations;
        neuronFileHeader.maxEMIterations = maxEMIterations;
        neuronFileHeader.emLikelihoodDelta = emLikelihoodDelta;


        neuronFile = new NeuronFile(neuronsFileName, neuronFileHeader, VisionParams.NEURONS_HEADER_CAPACITY,
                projectionsFile.getTTLTimes());


        File modelF = new File(datasetFolder + File.separator + datasetName + ".model");
        if (modelF.exists() && modelF.isFile() && modelF.canRead()) {
            modelFile = new ClusteringModelFile(modelF.getAbsolutePath());
        } else {
            VisionHeader modelFileHeader =
                new VisionHeader(neuronFileHeader);

            modelFileHeader.nlPointsEI = 10;
            modelFileHeader.nrPointsEI = 20;
            //modelFileHeader.minimizationInterval = 1;
            modelFileHeader.minimizationError = 1e-6;
            modelFileHeader.maxElectrodePatternSize = -1;
            modelFileHeader.minThreshold = prjFileHeader.threshold;

            modelFile = new ClusteringModelFile(modelF.getAbsolutePath(), VisionParams.MAX_MODEL_SLOTS,
                    modelFileHeader);
        }

        app.sendMessage("PCANF:Clustering, Create various structures...");

        // find neurons for electrodes in the core region (not overlapping extras)
        for (int electrode = fromElectrode; electrode <= toElectrode; electrode++) {
//			if (electrodeMap.isCoreElectrode(electrode)) electrodePool.add(electrode);
            electrodePool.add(electrode);
        }

        ClusteringThread[] clusteringThreads = new ClusteringThread[nClusteringThreads];
        for (int i = 0; i < clusteringThreads.length; i++) {
            clusteringThreads[i] = new ClusteringThread();
            clusteringThreads[i].start();
        }

        app.sendMessage("Clustering...");
        app.startProgressBar();
        for (boolean threadsStillAlive = true; threadsStillAlive; ) {
            // sleep a while
            Thread.sleep(1000);

            // then check if the threads are still alive
            threadsStillAlive = false;
            for (int i = 0; i < clusteringThreads.length; i++) {
                threadsStillAlive |= clusteringThreads[i].isAlive();
            }
        }
        app.endProgressBar();

        // all threads died, close the files
        neuronFile.close();
        modelFile.close();

        Vision.getInstance().getCalculationManager().calculationDone();
    }


    public static void densityClustering(
            float[][] data, int nSpikesOnElectrode, int binsPerDimension, int maxSpikeLoss,
            double clusteringSignificance, ArrayList<double[]> means, ArrayList<double[]> sigmas, double[] min,
            double[] max, double[] binInterval, boolean minimizeSpace, int minClusters,
            int maxClusters) {

        means.clear();
        sigmas.clear();
        final int nDimensions = data.length;

        // find the minimum and maximum limits on all the clustering dimensions
        for (int d = 0; d < nDimensions; d++) {
            min[d] = Double.POSITIVE_INFINITY;
            max[d] = Double.NEGATIVE_INFINITY;

            // find the maximum and minimum
            for (int i = 0; i < nSpikesOnElectrode; i++) {
                if (data[d][i] < min[d]) {
                    min[d] = data[d][i];
                }
                if (data[d][i] > max[d]) {
                    max[d] = data[d][i];
                }
            }
            //System.out.println("D: " + d + " Max: " + max[d] + " Min: " + min[d]);
            
            binInterval[d] = (max[d] - min[d]) / binsPerDimension;
        }

        // histogram the data points
        if (minimizeSpace) {
            HistogramND.minimizeSpace(data, nSpikesOnElectrode, min, max, binInterval,
                    maxSpikeLoss);
        }
        HistogramND densityHist = new HistogramND(min, max, binInterval);

        // fill the N-D histogram
        double[] x = new double[nDimensions];
        for (int i = 0; i < nSpikesOnElectrode; i++) {
            for (int d = 0; d < nDimensions; d++) {
                x[d] = data[d][i];
            }
            densityHist.fill(x, 1);
        }
        final int nonZeroBins = densityHist.countNonZero();
        double clusteringThreshold = nSpikesOnElectrode / ( (double) nonZeroBins) *
        clusteringSignificance;

        if (densityHist.getBinSum() < 2) {
            return;
        }

        // do the actual density clustering
        int[] maxBin = new int[nDimensions];

        while (densityHist.getMaxValueBin(maxBin) > clusteringThreshold) {
            double[] mean = new double[nDimensions];
            double[] sigma = new double[nDimensions];
            int n = densityHist.formDensityCluster(
                    maxBin, mean, sigma, 0, clusteringThreshold);
            if (n > 0) {
                for (int d = 0; d < nDimensions; d++) {
                    mean[d] = min[d] + (maxBin[d] + 0.5) * densityHist.getBinInterval(d);
                    sigma[d] = binInterval[d] * binInterval[d];
                }
                means.add(mean);
                sigmas.add(sigma);
            }
        } // clustering loop

        // make sure the model size is at least one, check for insufficient model size
        while (means.size() < minClusters) {
            double[] mean = new double[nDimensions];
            double[] sigma = new double[nDimensions];
            for (int d = 0; d < nDimensions; d++) {
                mean[d] = (min[d] + max[d]) * 0.5;
                sigma[d] = binInterval[d] * binInterval[d];
            }
            means.add(mean);
            sigmas.add(sigma);
        }
    }


    public void setParameters(HashMap<String, String> _p) {
        HashMap<String, String> p = _p;

        datasetFolder = new File(p.get("Dataset_Folder")).getAbsolutePath();
        fromElectrode = Integer.parseInt(p.get("From"));
        toElectrode = Integer.parseInt(p.get("To"));
        generateReport = Boolean.valueOf(p.get("Generate Report")).booleanValue();

        binsPerDimension = Integer.parseInt(p.get("Bins Per Dimension"));
        clusteringSignificance = Double.parseDouble(p.get("Clustering Significance"));
        minClusters = Integer.parseInt(p.get("Minimum Clusters"));
        maxClusters = Integer.parseInt(p.get("Miximum Clusters"));
        nEMSpikes = Integer.parseInt(p.get("Spikes Used For EM"));

        densityClusteringMaxSpikeLoss = Integer.parseInt(p.get(
        "Density Clustering Spike Loss"));
        minEMIterations = Integer.parseInt(p.get("Min EM Iterations"));
        maxEMIterations = Integer.parseInt(p.get("Max EM Iterations"));
        emLikelihoodDelta = Double.parseDouble(p.get("EM Likelihood Delta"));

        nClusteringThreads = Integer.parseInt(p.get("Clustering Threads"));
    }

}
