package edu.ucsc.neurobiology.vision.util;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.io.*;


/**
 * This class can be registered as a listener to the SpikeFinder and it will buffer
 * and sort the spikes chronologically before issuing them to it's own listeners.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class SpikeBuffer
    implements SpikeListener, SampleListener {

    private final SpikeListener[] spikeListeners;
    private int nListeners;
    public LinkedList<Spike> spikeList;
    private int currentSample = 0;
    private int minTime = Integer.MAX_VALUE;
    private Spike minTimeSpike = null;
    private final int timeWindow;


    public SpikeBuffer(int timeWindow) {
        this.spikeList = new LinkedList<Spike> ();
        this.spikeListeners = new SpikeListener[5];
        this.nListeners = 0;

        this.timeWindow = timeWindow;
    }


    public void addSpikeListener(SpikeListener listener) {
        if (nListeners == spikeListeners.length) {
            throw new IllegalArgumentException(
                                               "Only " + spikeListeners.length + " listenrs possible");
        }
        spikeListeners[nListeners] = listener;
        nListeners++;
    }


    public final void processSpike(final Spike _spike) {
        Spike spike = (Spike) _spike.clone();

        // add the spike to the list
        spikeList.add(spike);

        // update the min information
        if (spike.time < minTime) {
            minTime = spike.time;
            minTimeSpike = spike;
        }
    }


    public void processSample(short[] sample) {
        currentSample++;

        while (currentSample - minTime >= timeWindow && minTimeSpike != null) {
            // relase the min time spike and remove it from list
            for (int i = 0; i < nListeners; i++) {
                spikeListeners[i].processSpike(minTimeSpike);
            }
            spikeList.remove(minTimeSpike);

            // find the minimum spike time
            minTime = Integer.MAX_VALUE;
            minTimeSpike = null;
            final int n = spikeList.size();

            for (int i = 0; i < n; i++) {
                Spike spike = spikeList.get(i);
                if (spike.time < minTime) {
                    minTime = spike.time;
                    minTimeSpike = spike;
                }
            }
        }
    }


    public void finishSpikeProcessing() throws IOException {

    }


    public void finishSampleProcessing() throws IOException {
        // Patch for interrupted streaming failure //P
        // Remaining spikes in buffer when streaming falls
        // Throw an OutOfBound in the CovarianceCalculator
        // as right samples just don't exist for those spikes
        // See RunScript.neuronFinding()
                
        int failedSpikes = 0;
        
        // release all spikes
        while (minTimeSpike != null) {
            // relase the min time spike and remove it from list
            try { //P
                for (int i = 0; i < nListeners; i++) {
                    spikeListeners[i].processSpike(minTimeSpike);
                }
            } catch (IndexOutOfBoundsException e) { //P
                failedSpikes++; //P
            } //P
                
            spikeList.remove(minTimeSpike);

            // find the minimum spike time
            minTime = Integer.MAX_VALUE;
            minTimeSpike = null;
            final int n = spikeList.size();

            for (int i = 0; i < n; i++) {
                Spike spike = spikeList.get(i);
                if (spike.time < minTime) {
                    minTime = spike.time;
                    minTimeSpike = spike;
                }
            }
        }
        
        if (failedSpikes > 0) { //P
            System.out.println("Warning: " + failedSpikes + " spikes couldn't be processed at EOF."); //P
        } //P
        
        // tell the listeners we are done
        for (int i = 0; i < nListeners; i++) {
            spikeListeners[i].finishSpikeProcessing();
        }
    }

}
