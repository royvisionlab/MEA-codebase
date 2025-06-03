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
 * Calculation to map the spikes in a given dataset to the neurons identified in another
 * dataset in the same retina. Needs the master's dataset .model file to work. Produces
 * a .neurons file.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class PCANeuronMapping
        extends AbstractCalculation implements SampleListener {
    
    private String datasetFolder, rawFileName, prjPath, masterFolder, datasetName;
    private int nlPoints, nrPoints, nPoints;
    NeuronFile masterNeuronFile;
    private SpikeFile spikes;
    ElectrodeMap map;
    RawDataHeader header;
    VisionHeader h;
    boolean projectionsOnly;
    int nElectrodes;
    int nSamples;
    int nDimensions;
    ElectrodeUsage electrodeUsage;
    long sampleIndex;
    LinkedList<short[]> sampleList = new LinkedList<short[]>();
    int[][] adjacent;
    int[][] spikeTimes;
    short[] data;
    double[][][] eigenvectors;
    float[][][] projectedSpikes;
    double[] projectedSpike;
    ProjectionsFile pf;
    int[] spikeCount;
    SpikeAligner spikeAligner;
    double[][][] means, sigmas;
    int[][] neuronID;
    int[] nGaussians;
    double[][] pJ;
    IntegerList[] spikeLists;
    Vision app;
    long t1;
    float[] mean;
    SpikeIterator spikeIterator;
    Spike currentSpike;
    int nSpikes;
    int[] idList;
    double f;


    public void startCalculation() throws Exception {
        app = Vision.getInstance();
        app.sendMessage("Preparing...");

        // get the raw data file
        MultipleCompressedSampleInputStream sis = new MultipleCompressedSampleInputStream(rawFileName);
        header = sis.getHeader();
        this.nElectrodes = header.getNumberOfElectrodes();
        this.nSamples = header.getNumberOfSamples();
        
        datasetName = new File(datasetFolder).getName();
        String masterDatasetName = new File(masterFolder).getName();
        
        // Get the spikes and neurons
        spikes = new SpikeFile(datasetFolder + File.separator + datasetName + VisionParams.SPIKES_FILE_EXTENSION);
        map = ElectrodeMapFactory.getElectrodeMap(spikes.getArrayID());
        masterNeuronFile = new NeuronFile(masterFolder + File.separator + masterDatasetName + VisionParams.NEURON_FILE_EXTENSION);
        idList = masterNeuronFile.getIDList();
        final int maxID = MathUtil.max(masterNeuronFile.getFullIDList()) * 10;
        spikeLists = new IntegerList[maxID + 1];
        for (int i = 0; i < idList.length; i++)
            spikeLists[idList[i]] = new IntegerList(10000, 2);
        
        spikeIterator = spikes.iterator();
        currentSpike = spikeIterator.next();
        nSpikes = spikes.getSpikesCount();
        n++;
        
        mean = new float[nElectrodes];
        
        // load the model
        String modelPath = masterFolder + File.separator + masterDatasetName + VisionParams.MODEL_FILE_EXTENSION;
        loadModel(modelPath);
        // figure out the name of the output projections file
        prjPath = datasetFolder + File.separator + datasetName + VisionParams.PROJ_FILE_EXTENSION;
        // done here because of speed reasons
        f = Math.pow(2 * Math.PI, nDimensions / 2.0);

        for (int i = 0; i < nPoints; i++) sampleList.add(new short[nElectrodes]);

        // create the adjacency table
        adjacent = new int[nElectrodes][];
        if (electrodeUsage == ElectrodeUsage.ONE_ELECTRODE) {
            for (int electrode = 0; electrode < nElectrodes; electrode++)
                adjacent[electrode] = new int[] {electrode};
        } else {
            for (int electrode = 0; electrode < nElectrodes; electrode++)
                adjacent[electrode] = map.getAdjacentsTo(electrode);
        }

        
        //TOBEDONE
        //        spikeAlligner = new SpikeAlligner(nPoints, nlPoints, adjacent,
                //                                          masterNeuronFile.getHeader().minimizationError);
        spikeAligner = new SpikeAligner(nPoints, nlPoints, adjacent, 0.001);

        data = new short[nPoints * 7];
        spikeCount = new int[nElectrodes];
        projectedSpike = new double[nDimensions];


        // make a new projections file if necessary
        if (projectionsOnly) {
            int[] nSpikesElectrode = new int[nElectrodes];
            int maxSpikes = 0;
            for (int i=0; i<nElectrodes; i++) {
                spikeCount[i] = 0;
                nSpikesElectrode[i] = spikes.getSpikesCount(i);
                if (nSpikesElectrode[i] > maxSpikes)
                    maxSpikes = nSpikesElectrode[i];
            }
            
            // build a projections file using the header from the master model file
            pf = new ProjectionsFile(prjPath, h, nSpikesElectrode, spikes.getTTLTimes());
            projectedSpikes = new float[nElectrodes][ProjectionsFile.writeBufferSize][nDimensions];
            spikeTimes = new int[nElectrodes][ProjectionsFile.writeBufferSize];
        }

        // start
        app.sendMessage("Mapping...");
        System.out.println("Mapping...");
        app.startProgressBar();

        t1 = System.currentTimeMillis();
        sis.addSampleListener(this);
        sis.start();
    }


    private final double timeConstant = 0.010; // in seconds
    private final float alpha = (float) (1.0 / (timeConstant * 20000.0));
    long n = 0;
    int percentage, oldPercentage = -1;
    int oldTime = -1;


    public void processSample(short[] sample) {
        short[] s = (short[]) sampleList.removeFirst();
        System.arraycopy(sample, 0, s, 0, nElectrodes);
        for (int i = 0; i < nElectrodes; i++) {
            mean[i] += alpha * (s[i] - mean[i]);
            s[i] -= mean[i];
        }
        sampleList.addLast(s);

        while (sampleIndex - currentSpike.time == nrPoints) {
            Iterator<short[]> iter = sampleList.iterator();
            for (int i = 0; iter.hasNext(); i++) {
                s = (short[]) iter.next();
                for (int adj = 0; adj < adjacent[currentSpike.electrode].length; adj++) {
                    data[adj * nPoints + i] = s[adjacent[currentSpike.electrode][adj]];
                }
            }

            // we have a complete spike for this electrode, process it
            if (sampleIndex >= nPoints) {
                try {
                    processSpike(currentSpike.electrode, currentSpike.time, data);
                } catch (CannotEvaluateException e) {
                    e.printStackTrace();
                }
            }
            if (!spikeIterator.hasNext()) break;
            nextSpike();
        }

        sampleIndex++;
    }


    private void nextSpike() {
        // get the next spike, but ignore TTL signals
        do {
            currentSpike = spikeIterator.next();
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


    private final void processSpike(final int electrode, int time, short[] spike) throws CannotEvaluateException {
        if (eigenvectors[electrode] == null) return;

        // align the spike time
        double[] alignedSpike = spikeAligner.align(electrode, spike, true);

        // get the projection on the PCA eigenvectors
        for (int d = 0; d < nDimensions; d++) {
            projectedSpike[d] = MathUtil.project(eigenvectors[electrode][d], alignedSpike);
            if (projectionsOnly) {
                projectedSpikes[electrode][spikeCount[electrode]][d] = (float) projectedSpike[d];
                spikeTimes[electrode][spikeCount[electrode]] = time;
            }
        }

        // save projections if necessary
        if (projectionsOnly) {
            try {
                spikeCount[electrode]++;

                // save the projections out (in batches less than BUFFER SIZE)
                if (spikeCount[electrode] == ProjectionsFile.writeBufferSize) {
                    pf.saveData(electrode, spikeCount[electrode], spikeTimes, projectedSpikes);
                    spikeCount[electrode] = 0;
                }

            } catch(IOException e) { Vision.reportException(e); }
        } else {
            // get the closest gaussian
            double denom = 0;
            for (int j = 0; j < nGaussians[electrode]; j++)
                denom += pJ[electrode][j] * ExpectationMaximization._pXJ(projectedSpike, means[electrode][j], sigmas[electrode][j], f);

            double maxP = Double.NEGATIVE_INFINITY;
            int closestGaussian = -1;
            for (int j = 0; j < nGaussians[electrode]; j++) {
                double p = pJ[electrode][j] * ExpectationMaximization._pXJ(
                        projectedSpike, means[electrode][j], sigmas[electrode][j], f) / denom;
                if (p > maxP) {
                    maxP = p;
                    closestGaussian = j;
                }
            }

            if (closestGaussian != -1) {
                int id = neuronID[electrode][closestGaussian];
                if (spikeLists[id] != null) spikeLists[id].add(time); // only good neurons are saved
            }
        }
    }


    public void finishSampleProcessing() throws IOException {
        app.endProgressBar();

        String name = datasetFolder + File.separator + datasetName + ".neurons";
        VisionHeader h = masterNeuronFile.getHeader();

        if (!projectionsOnly) {
            NeuronFile reconstructedNeuronFile = new NeuronFile(name, h, VisionParams.NEURONS_HEADER_CAPACITY, spikes.getTTLTimes());

            for (int i = 0; i < idList.length; i++) {
                int neuronID = idList[i];
                reconstructedNeuronFile.addNeuron(
                    masterNeuronFile.getElectrode(neuronID), neuronID,
                    spikeLists[neuronID], spikeLists[neuronID].size());
            }

            reconstructedNeuronFile.close();
        } else {
            try {

                // save the projections out (in batches less than BUFFER SIZE)
                for (int electrode = 0; electrode < nElectrodes; electrode++) {
                    if ( (eigenvectors[electrode] != null) && (spikeCount[electrode] != 0))
                        pf.saveData(electrode, spikeCount[electrode], spikeTimes, projectedSpikes);
                }
                pf.close();

            } catch(IOException e) { Vision.reportException(e); }
        }

        long t2 = System.currentTimeMillis();
        app.sendMessage("Done in: " + (t2 - t1) / 1000. + " s.");
        Vision.getInstance().getCalculationManager().calculationDone();
    }


    private void loadModel(String modelFileName) throws IOException {
        ClusteringModelFile modelFile = new ClusteringModelFile(modelFileName);

        // read and discard nElectrodes
        h = modelFile.getUserHeader();
        this.nDimensions = h.nDimensions;
        this.nlPoints = h.nlPoints;
        this.nrPoints = h.nrPoints;
        this.nPoints = nlPoints + nrPoints + 1;
        this.electrodeUsage = h.electrodeUsage;

        means = new double[nElectrodes][][];
        sigmas = new double[nElectrodes][][];
        neuronID = new int[nElectrodes][];
        nGaussians = new int[nElectrodes];
        pJ = new double[nElectrodes][];

        for (int electrode = 1; electrode < nElectrodes; electrode++) {
            ClusteringModelFile.EMModel m =
                (ClusteringModelFile.EMModel) modelFile.getNeuronExtraction(electrode);

            if (m != null) {
                nGaussians[electrode] = m.nClusters; //??? or m.nGaussians
                        neuronID[electrode] = m.neuronID;
                pJ[electrode] = m.probability;
                means[electrode] = m.means;
                sigmas[electrode] = m.covariances;
            }
        }

        eigenvectors = new double[nElectrodes][][];
        int[] elec = modelFile.getExtractionIds();
        for (int i = 0; i < elec.length; i++) {
            int electrode = elec[i];
            ClusteringModelFile.Model m = modelFile.getNeuronExtraction(electrode);

            eigenvectors[electrode] = new double[nDimensions][];
            for (int d = 0; d < nDimensions; d++)
                eigenvectors[electrode][d] = m.eigenvectors[d];
        }

        modelFile.close();
    }


    public void setParameters(HashMap<String,String> parameters) {
        masterFolder  = new File(parameters.get("Master Dataset Folder")).getAbsolutePath();
        datasetFolder = new File(parameters.get("Dataset Folder")).getAbsolutePath();
        rawFileName   = new File(parameters.get("Raw Data File")).getAbsolutePath();
        projectionsOnly = Boolean.parseBoolean(parameters.get("Projections Only"));
    }
}