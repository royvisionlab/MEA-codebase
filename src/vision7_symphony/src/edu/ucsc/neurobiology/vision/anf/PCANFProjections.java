package edu.ucsc.neurobiology.vision.anf;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.calculations.*;
import edu.ucsc.neurobiology.vision.electrodemap.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * Calculation that uses the .spikes file and the .cov file to project all spikes on the
 * first 5 PCA components obtained from the spike shapes on every electrode. Saves the
 * projections to the .prj file.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class PCANFProjections
    extends AbstractCalculation implements SampleListener {

    public int nDimensions;

    float[][][] spikeProjections;
    int[][] spikeTimes;
    private String spikeFileName, rawDataSource, covFileName, modelFileName,
        projectionsFileName;
    ElectrodeMap map;
    private int nlPoints, nrPoints, nPoints;
//    ElectrodeUsage electrodeUsage = null;
    RawDataHeader rawDataHeader;
    int nElectrodes;
    long nSpikes;
    long sampleIndex;
    LinkedList<short[]> sampleList;
    int[][] adjacent;
    int[] spikeIndex;
    ProjectionsFile projectionsFile;
    double[][][] eigenvectors;
    Vision app;
    private SpikeFile spikes;
    SpikeIterator spikeIterator;
    Spike currentSpike;
    float[] mean;

    private final double timeConstant = 0.010; // in seconds
    private final float alpha = (float) (1.0 / (timeConstant * 20000.0));
    private long n = 1; // DO NOT CHANGE, the algorithm needs a starting value of 1
    int percentage, oldPercentage = -1;
    int oldTime = -1;
    int[] data;
    ProjectingThread[] projThreads;
    ParallelUtil.StartEndIndices[] electrodeBins;
    int nThreads;
    int maxAdjacentElectrodes;

    VisionHeader header;


    public void startCalculation() throws Exception {
        app = Vision.getInstance();
        app.sendMessage("Preparing...");

        // get the raw data file
        MultipleCompressedSampleInputStream sampleInputStream;
        sampleInputStream = new MultipleCompressedSampleInputStream(rawDataSource);
        rawDataHeader = sampleInputStream.getHeader();
        this.nElectrodes = rawDataHeader.getNumberOfElectrodes();
        mean = new float[nElectrodes];

        // Get the spikes
        spikes = new SpikeFile(spikeFileName);
        map = ElectrodeMapFactory.getElectrodeMap(spikes.getArrayID());
        spikeIterator = spikes.iterator();
        nSpikes = spikes.getSpikesCount();
        nextSpike();

        // load eigenvectors (either from .cov or .model)
        eigenvectors = new double[nElectrodes][][];
        if (modelFileName != null && modelFileName.trim().length() != 0) {
            // a model file was provided, load covariances from it
            loadModel(modelFileName);

            // bug fix: shlens
            // correct header if mapping from different data set
            header.nSamples = spikes.getNumberOfSamples();
        } else {
            // load covariances from the .cov file
            loadCovariances(covFileName);
        }
        nPoints = nlPoints + nrPoints + 1;

        // get the maximum number of adjacents
        maxAdjacentElectrodes = 0;
        adjacent = new int[nElectrodes][];
        for (int electrode = 0; electrode < nElectrodes; electrode++) {
            adjacent[electrode] = map.getAdjacentsTo(electrode,
                header.electrodeUsage.ordinal());
            if (adjacent[electrode].length > maxAdjacentElectrodes) {
                maxAdjacentElectrodes = adjacent[electrode].length;
            }
        }

        sampleList = new LinkedList<short[]>();
        for (int i = 0; i < nPoints; i++) {
            sampleList.add(new short[nElectrodes]);
        }
        data = new int[nPoints * maxAdjacentElectrodes + 2];

        int[] spikeCount = new int[nElectrodes];
        for (int i = 0; i < nElectrodes; i++) {
            if (eigenvectors[i] != null) {
                spikeCount[i] = spikes.getSpikesCount(i);
            } else {
                spikeCount[i] = 0;
            }
        }

        // fill the header
        header.nDimensions = nDimensions;

        // create the .prj file
        projectionsFile = new ProjectionsFile(projectionsFileName, header,
                                              spikeCount, spikes.getTTLTimes());

        spikeProjections = new float[nElectrodes][ProjectionsFile.writeBufferSize][nDimensions];
        spikeTimes = new int[nElectrodes][ProjectionsFile.writeBufferSize];
        spikeIndex = new int[nElectrodes];

        projThreads = new ProjectingThread[nThreads];
        electrodeBins = new ParallelUtil.StartEndIndices[nThreads];
        for (int i = 0; i < nThreads; i++) {
            projThreads[i] = new ProjectingThread();
            electrodeBins[i] = ParallelUtil.subIndex(i, nThreads, nElectrodes);
        }

        sampleInputStream.addSampleListener(this);

        app.startProgressBar();
        app.sendMessage("Projecting...");
        
        sampleInputStream.start();
        for (ProjectingThread pt : projThreads) pt.start();
    }


    public final void processSample(short[] sample) {
        if (currentSpike == null) return;

        short[] s = (short[]) sampleList.removeFirst();
        System.arraycopy(sample, 0, s, 0, nElectrodes);
        for (int i = 0; i < nElectrodes; i++) {
            mean[i] += alpha * (s[i] - mean[i]);
            s[i] -= mean[i];
        }
        sampleList.addLast(s);

        // while because spikes on different electrodes may have the same time
        while (sampleIndex - currentSpike.time == nrPoints) {
            for (int i = 0; i < nPoints; i++) {
                s = (short[]) sampleList.get(i);
                for (int adj = 0; adj < adjacent[currentSpike.electrode].length; adj++) {
                    data[adj * nPoints + i] = s[adjacent[currentSpike.electrode][adj]];
                }
            }
            
            if (sampleIndex >= nPoints) {
                data[data.length - 1] = currentSpike.electrode;
                data[data.length - 2] = currentSpike.time;

                for (int j = 0; j <= electrodeBins.length; j++) {
                    if (currentSpike.electrode <= electrodeBins[j].end) {
                        projThreads[j].add(data);
                        break;
                    }
                    
                    // This should never happen
                    if (j == electrodeBins.length) throw new Error("Failed to assign electrode spike to ProjectingThread!");
                }
            }
            
            nextSpike();
            if (currentSpike == null) {
                return;
            }
        }

        sampleIndex++;
    }


    private void nextSpike() {
        currentSpike = null;

        // get the next spike, but ignore TTL signals
        do {
            currentSpike = spikeIterator.next();
            if (currentSpike == null) {
                return;
            }
            n++;

            if (oldTime != -1) {
                if (currentSpike.time < oldTime) {
                    System.out.println("Sorting problem:");
                    System.out.println("t: " + oldTime + ", " + currentSpike.time);
                }
            }
            oldTime = currentSpike.time;

            percentage = (int) (100L * n / nSpikes);
            if (percentage != oldPercentage) {
                app.setProgress(percentage);
                oldPercentage = percentage;
            }
        } while (currentSpike.electrode == 0);
    }


    public void finishSampleProcessing() throws IOException {
        app.endProgressBar();
        for (ProjectingThread pt : projThreads) pt.finish();
        while (anyAlive(projThreads)) {
            try {
                Thread.sleep(100);
            } catch (InterruptedException e) {}
        }

        // finish saving all data
        for (int i = 0; i < nElectrodes; i++) {
            if ( (eigenvectors[i] != null) && (spikeIndex[i] == 0)) {
                System.out.println("Warning: electrode " + i +
                                   " has zero spikes, check thresholds, they may be too high");
            }

            if ((eigenvectors[i] != null) && (spikeIndex[i] != 0)) {
                projectionsFile.saveData(i, spikeIndex[i], spikeTimes, spikeProjections);
            }
        }
        // close the file
        projectionsFile.close();

        app.getCalculationManager().calculationDone();
    }


    private void loadCovariances(String covariancesFileName) throws Exception {
        CovarianceFile covFile = new CovarianceFile(covariancesFileName);
        header = covFile.getHeader();
        nlPoints = header.nlPoints;
        nrPoints = header.nrPoints;

        ElectrodeMap m = ElectrodeMapFactory.getElectrodeMap(covFile.getArrayID());
        int nElectrodes = m.getNumberOfElectrodes();

        for (int electrode = 1; electrode < nElectrodes; electrode++) {
            float[] cov = covFile.getCovarianceMatrix(electrode);
            if (cov == null) continue; // We should have complained about this already in the NWCovariance step
            
            PCA pca = new PCA(cov);
            try {
                pca.doPCA();
            } catch (TooManyIterationsException e) {
                System.err.println("Skipping electrode " + electrode + "; too many iterations in PCA calculation");
                continue;
            }
            
            eigenvectors[electrode] = new double[nDimensions][];
            for (int d = 0; d < nDimensions; d++)
                eigenvectors[electrode][d] = pca.getEigenVector(d);
        }

        covFile.close();
    }


    private void loadModel(String fileName) throws Exception {
        ClusteringModelFile f = new ClusteringModelFile(fileName);

        header = f.getUserHeader();
        nlPoints = header.nlPoints;
        nrPoints = header.nrPoints;

        int[] elec = f.getExtractionIds();

        for (int i = 0; i < elec.length; i++) {
            int electrode = elec[i];
            ClusteringModelFile.Model m = f.getNeuronExtraction(electrode);

            eigenvectors[electrode] = new double[nDimensions][];
            for (int d = 0; d < nDimensions; d++) {
                eigenvectors[electrode][d] = m.eigenvectors[d];
            }
        }
    }


    static private boolean anyAlive(ProjectingThread[] projThreads) {
        for (ProjectingThread pt : projThreads)
            if (pt.isAlive()) return true;
        return false;
    }
    
    
    class ProjectingThread extends Thread {

        private SpikeAligner aligner;
        private LinkedList<int[]> fullList, emptyList;
        private boolean finished = false;

        public ProjectingThread() {
            aligner = new SpikeAligner(nPoints, nlPoints, adjacent,
                                         header.minimizationError);
            emptyList = new LinkedList<int[]>();
            fullList = new LinkedList<int[]>();
            for (int i = 0; i < 10000; i++) {
                emptyList.add(new int[nPoints * maxAdjacentElectrodes + 2]);
            }
        }


        public void add(int[] spike) {
            int[] d;

            synchronized (emptyList) {
                while (emptyList.isEmpty()) {
                    try {
                        emptyList.wait();
                    } catch (InterruptedException e) {}
                }
                d = emptyList.removeFirst();
            }

            System.arraycopy(spike, 0, d, 0, spike.length);

            synchronized (fullList) {
                fullList.addLast(d);
                fullList.notifyAll();
            }
        }


        public void finish() {
            synchronized (fullList) {
                finished = true;
                fullList.notifyAll();
            }
        }
        

        public void run() {
            int[] spike;
            int n = 0;

            while (true) {
                synchronized (fullList) {
                    while (fullList.isEmpty()) {
                        if (finished) {
                            return;
                        }
                        try {
                            fullList.wait();
                        } catch (InterruptedException e) {}
                    }
                    spike = fullList.removeFirst();
                }

                // the analysis of the spike
                final int electrode = spike[spike.length - 1];
                final int time = spike[spike.length - 2];

                try {
                    if (eigenvectors[electrode] != null) {
                        // align the spike
                        double[] alignedSpike = aligner.align(electrode, spike, true);
                        // project onto eigenvectors
                        float[] proj = spikeProjections[electrode][spikeIndex[electrode]];
                        for (int d = 0; d < nDimensions; d++) {
//                            proj[d] = (float) pca[electrode].project(allignedSpike, d);
                            proj[d] = (float) MathUtil.project(eigenvectors[electrode][d], alignedSpike);
                        }
                        spikeTimes[electrode][spikeIndex[electrode]] = time;

                        // save if needed
                        spikeIndex[electrode]++;
                        if (spikeIndex[electrode] == projectionsFile.writeBufferSize) {
                            try { // BEU
                                projectionsFile.saveData(
                                    electrode, spikeIndex[electrode], spikeTimes,
                                    spikeProjections);
                            } catch (IOException e) {
                                e.printStackTrace();
                            }
                            spikeIndex[electrode] = 0;
                        }
                    }
                } catch (CannotEvaluateException e) {
                    e.printStackTrace();
                }

                synchronized (emptyList) {
                    emptyList.addLast(spike);
                    emptyList.notifyAll();
                }
            }
        }

    }


    public void setParameters(HashMap<String, String> parameters) {
        rawDataSource = parameters.get("Raw_Data_File");
        modelFileName = parameters.get("Model File");
        String datasetFolder = parameters.get("Dataset Folder");
        nDimensions = Integer.parseInt(parameters.get("PCA Dimensions"));
        nThreads = Integer.parseInt(parameters.get("nThreads"));
        
        String datasetName = new File(datasetFolder).getName();
        spikeFileName = datasetFolder + File.separator + datasetName + ".spikes";
        covFileName = datasetFolder + File.separator + datasetName + ".wcov";
        File covF = new File(covFileName);
        //use standard covariance file if whitened covariance file was not calculated.
        if (!covF.exists()) {
            covFileName = datasetFolder + File.separator + datasetName + ".cov";
        }

        // shlens addition; in case user wants to specify the projections file
        projectionsFileName = parameters.get("Projections File");

        
        // use default projections file name if none provided
        if (projectionsFileName == null ||
            projectionsFileName.trim().length() == 0) {
            projectionsFileName = datasetFolder + File.separator + datasetName + VisionParams.PROJ_FILE_EXTENSION;
        }
    }

}
