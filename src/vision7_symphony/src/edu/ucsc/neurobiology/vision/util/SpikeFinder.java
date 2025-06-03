package edu.ucsc.neurobiology.vision.util;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.electrodemap.*;
import edu.ucsc.neurobiology.vision.io.*;


/**
 * A spike finder which finds spikes according to a fixed threshold and using
 * an updating mean. The spikes generated are NOT GUARANTEED TO BE IN CHRONOLOGICAL ORDER.
 * Use a SpikeBuffer to sort the spikes.
 *
 * @see SpikeBuffer
 * @author Dumitru Petrusca, University of California, Santa Cruz <br>
 *         Charles A. Loomis, University of California, Santa Cruz
 */
public class SpikeFinder
    implements SampleListener, SampleInitializationListener {

    private final int nElectrodes;
    private final boolean[] buildingSpike;
    private final float[] maxAmplitude;
    private final int[] maxTime, previousSpikeTime;
    private final int[] startTime;
    private final float[] spikeThreshold;
    private final float ttlThreshold;
    private final boolean[] disconnected;
    private final static float delta = 50e-6f;
    private final float alpha;
    private float[] means;
    private int electrode;
    private int currentSample = -1; // important to be -1
    private float amplitude;
    private Spike spike = new Spike(0, 0, 0);
    private SpikeListener spikeListener1, spikeListener2, spikeListener3, spikeListener4;
    public int nSpikes;
    private boolean initialized;

    double previousTTL = -1;
    double average = 0;
    int nTTL = 0;

    private static final int minTimeSeparation = (int) (0.25 * 20); // .25 ms
    public static final int maxSpikeWidth = (int) (2.5 * 20); // 2.5 ms


    public SpikeFinder(
        ElectrodeMap electrodeMap, float[] spikeThreshold, float ttlThreshold, float timeConstant) {

        this.spikeThreshold = spikeThreshold;
        this.ttlThreshold = ttlThreshold;
        this.nElectrodes = electrodeMap.getNumberOfElectrodes();
        this.disconnected = electrodeMap.getDisconnectedElectrodesList();

        buildingSpike = new boolean[nElectrodes];
        maxTime = new int[nElectrodes];
        previousSpikeTime = new int[nElectrodes];
        Arrays.fill(previousSpikeTime, -1000);
        maxAmplitude = new float[nElectrodes];
        startTime = new int[nElectrodes];
        means = new float[nElectrodes];

        this.alpha = delta / timeConstant;
    }


    public void addSpikeListener(SpikeListener listener) {
        if (spikeListener1 == null) {
            spikeListener1 = listener;
        } else if (spikeListener2 == null) {
            spikeListener2 = listener;
        } else if (spikeListener3 == null) {
            spikeListener3 = listener;
        } else if (spikeListener4 == null) {
            spikeListener4 = listener;
        } else {
            throw new IllegalStateException("Only 4 listeners allowed");
        }
    }


    public synchronized void removeSpikeListener(SpikeListener listener) {
        if (spikeListener1 == listener) {
            spikeListener1 = null;
        } else if (spikeListener2 == listener) {
            spikeListener2 = null;
        } else if (spikeListener3 == listener) {
            spikeListener3 = null;
        } else if (spikeListener4 == listener) {
            spikeListener4 = null;
        }
    }


    //    DoubleHistogram h = new DoubleHistogram("", 0, 1000, 1);

    public final void processSample(short[] sample) {
        if (!initialized) {
            throw new Error("The initialize method was not called");
        }

        currentSample++;

        // Treat the TTL signal specially.
        if (sample[0] < -ttlThreshold) {
            // Make a new spike the first time the threshold goes below ttlThreshold counts.
            if (!buildingSpike[0]) {
                buildingSpike[0] = true;

                spike.electrode = 0;
                spike.time = currentSample;
                spike.amplitude = 1500;

                if (spikeListener1 != null) {
                    spikeListener1.processSpike(spike);
                }
                if (spikeListener2 != null) {
                    spikeListener2.processSpike(spike);
                }
                if (spikeListener3 != null) {
                    spikeListener3.processSpike(spike);
                }
                if (spikeListener4 != null) {
                    spikeListener4.processSpike(spike);
                }

                if (previousTTL > 0) {
                    average += currentSample - previousTTL;
                    nTTL++;
                }

                previousTTL = currentSample;
                nSpikes++;
            }
        } else {
            buildingSpike[0] = false;
        }

        // If this isn't the TTL signal, process the spike normally.
        for (electrode = 1; electrode < nElectrodes; electrode++) {
            if (disconnected[electrode]) {
                continue;
            }

            // data not filtered
            amplitude = sample[electrode];
            means[electrode] += alpha * (amplitude - means[electrode]);
            amplitude = means[electrode] - amplitude;

            // data filtered
//            amplitude = -sample[electrode];

            if (amplitude > spikeThreshold[electrode]) {
                if (!buildingSpike[electrode]) {
                    // Create a new spike
                    maxTime[electrode] = currentSample;
                    startTime[electrode] = currentSample;
                    maxAmplitude[electrode] = amplitude;
                    buildingSpike[electrode] = true;
                } else {
                    // Add information to existing spike
                    if (amplitude > maxAmplitude[electrode]) {
                        maxTime[electrode] = currentSample;
                        maxAmplitude[electrode] = amplitude;
                    }
                }
            } else if (buildingSpike[electrode]) {
                // The voltage dropped bellow threshold, the spike is finished
                // Save the information about the spike and prepare for the next one

//                h.fill(currentSample - startTime[electrode], 1);

                if (maxTime[electrode] - previousSpikeTime[electrode] > minTimeSeparation
                    && currentSample - startTime[electrode] <= maxSpikeWidth) {

                    spike.electrode = electrode;
                    spike.time = maxTime[electrode];
                    spike.amplitude = maxAmplitude[electrode];

                    if (spikeListener1 != null) {
                        spikeListener1.processSpike(spike);
                    }
                    if (spikeListener2 != null) {
                        spikeListener2.processSpike(spike);
                    }
                    if (spikeListener3 != null) {
                        spikeListener3.processSpike(spike);
                    }
                    if (spikeListener4 != null) {
                        spikeListener4.processSpike(spike);
                    }
                    nSpikes++;
                }

                previousSpikeTime[electrode] = maxTime[electrode];
                buildingSpike[electrode] = false;
            }
        } // for (electrode
    }


    public void initialize(short[][] sampleBuffer) {
        if (initialized) {
            throw new Error("Multiple initialization, surely a bug !");
        }

        if (currentSample != -1) {
            throw new Error("processSample() might have been called before initialize");
        }

        for (int i = 0; i < sampleBuffer.length; i++) {
            for (int elec = 0; elec < nElectrodes; elec++) {
                means[elec] += sampleBuffer[i][elec];
            }
        }

        for (int elec = 0; elec < nElectrodes; elec++) {
            means[elec] /= sampleBuffer.length;
//            System.err.println(elec + ": " + means[elec]);
        }

        initialized = true;
    }
    
    public void initialize() {
        initialize(new short[0][]);
    }


    public double getRefreshRate() {
        return average / nTTL;
    }


    /**
     * The spike finder is not doing any buffering, therefore it just passes the
     * "finish" signal directly to listeners.
     */
    public void finishSampleProcessing() {
        if (spikeListener1 != null) {
            try {
                spikeListener1.finishSpikeProcessing();
            } catch (IOException e) {
                Vision.reportException("Could not finish processing for " +
                                       spikeListener1.getClass().getName(), e);
            }
        }
        if (spikeListener2 != null) {
            try {
                spikeListener2.finishSpikeProcessing();
            } catch (IOException e) {
                Vision.reportException("Could not finish processing for " +
                                       spikeListener2.getClass().getName(), e);
            }
        }
        if (spikeListener3 != null) {
            try {
                spikeListener3.finishSpikeProcessing();
            } catch (IOException e) {
                Vision.reportException("Could not finish processing for " +
                                       spikeListener2.getClass().getName(), e);
            }
        }
        if (spikeListener4 != null) {
            try {
                spikeListener4.finishSpikeProcessing();
            } catch (IOException e) {
                Vision.reportException("Could not finish processing for " +
                                       spikeListener4.getClass().getName(), e);
            }
        }
    }


    public final static int findSpikes(
        short[] data, double threshold, double[] spikeTimeList,
        float[] spikeAmplitudeList, final int ignoreAtBeginning, final int ignoreAtEnd) {

        int maxTime = 0, previousSpikeTime = -1000;
        double amplitude;
        int maxAmplitude = 0, startTime = 0;
        int sIndex = 0;
        boolean buildingSpike = false;

        // calculate the initial mean from the 1st second of data
        double mean = 0;
        int n = 20000;
        for (int i = ignoreAtBeginning; i < ignoreAtBeginning + n; i++) {
            mean += data[i];
        }
        mean /= n;

        double timeConstant = 0.01;
        double alpha = delta / timeConstant;

        // ignore the filtering artifact
        for (int currentSample = ignoreAtBeginning;
                                 currentSample < data.length - ignoreAtEnd; currentSample++) {

            // data not filtered
            amplitude = data[currentSample];
            mean += alpha * (amplitude - mean);
            amplitude = mean - amplitude;

            // data filtered
//            amplitude = -data[currentSample];

            if (amplitude > threshold) {
                if (!buildingSpike) {
                    // Create a new spike
                    buildingSpike = true;
                    maxTime = currentSample;
                    startTime = currentSample;
                    maxAmplitude = (int) amplitude;
                } else {
                    // Add information to existing spike
                    if (amplitude > maxAmplitude) {
                        maxTime = currentSample;
                        maxAmplitude = (int) amplitude;
                    }
                }
            } else if (buildingSpike) {
                // The spike is finished
                // Save the information about the spike and prepare for the next one
                if (maxTime - previousSpikeTime > minTimeSeparation
                    && currentSample - startTime <= maxSpikeWidth) {

                    spikeTimeList[sIndex] = maxTime;
                    spikeAmplitudeList[sIndex] = maxAmplitude;
                    sIndex++;
                }

                previousSpikeTime = maxTime;
                buildingSpike = false;
            }
        }

        return sIndex;
    }

    public int numSamplesProcessed() {
        return currentSample + 1;
    }
    
    public int numSpikesFound() {
        return nSpikes;
    }
}
