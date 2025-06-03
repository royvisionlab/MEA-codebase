package edu.ucsc.neurobiology.vision.util;

import java.io.IOException;

import edu.ucsc.neurobiology.vision.io.SampleListener;
import edu.ucsc.neurobiology.vision.io.Spike;
import edu.ucsc.neurobiology.vision.io.SpikeListener;
import edu.ucsc.neurobiology.vision.math.*;

public class NoiseFinder implements SpikeListener, SampleListener {

    private final SpikeListener[] noiseListeners;
    private int nListeners;
    private int nElectrodes;
    private int minimumNoiseEvents;
    private double eventProbabilityPerSample;
    private int voidDistance; //minimum distance between noise "spike" and real spike, or between two noise "spikes"
    private int[] readySample; // time till next possible noise "spike"
    private int currentSample;
    private int[] eventsGenerated; //count of how many noise "spikes" have been generated.
    RandomMersenne rand = new RandomMersenne(4357);
    final static double EVENTS_MULTIPLIER = 2; //because we can't be sure of the correct probability to get the desired number of events
                                            //we multiply by a constant, to get extra.
    
    public NoiseFinder(int nElectrodes, int nSamples, int minimumNoiseEvents, int voidDistance, int spikeBufferSize ) {
        this.nElectrodes = nElectrodes;
        this.minimumNoiseEvents = minimumNoiseEvents;
        this.eventProbabilityPerSample = EVENTS_MULTIPLIER* (double) minimumNoiseEvents/ (double) nSamples;
        this.voidDistance = voidDistance;
        this.noiseListeners = new SpikeListener[5];
        this.nListeners = 0;
        this.readySample = new int[nElectrodes];   
        this.eventsGenerated = new int[nElectrodes];
        currentSample = -spikeBufferSize + 1; //spike buffer causes a delay between spikes and samples.  This corrects it.		
    }

    
    public void addNoiseListener(SpikeListener listener) {
        if (nListeners == noiseListeners.length){
            throw new IllegalArgumentException(
                    "Only " + noiseListeners.length + " listeners possible.");
        }
        noiseListeners[nListeners] = listener;
        nListeners++;
    }

    
    public void processSpike(Spike spike) {
        readySample[spike.electrode] = spike.time + voidDistance;
//		System.out.println("spike.time:" + spike.time);
    }


    public void processSample(short[] sample) {
//		System.out.println("currentSample: " + currentSample);
        for (int electrode = 1; electrode < nElectrodes; electrode++) {
            if (currentSample >= readySample[electrode] + voidDistance) {
                if (rand.nextDouble() < eventProbabilityPerSample) {
                    Spike noise = new Spike(currentSample - voidDistance, electrode, 0);
                    for (int listener = 0; listener < nListeners; listener++) {
                        noiseListeners[listener].processSpike(noise);
                    }
                    eventsGenerated[electrode]++;
                    readySample[electrode] = currentSample + voidDistance;
                }
            }
        }
        currentSample++;
    }

    
    public void finishSpikeProcessing() throws IOException {
    }

    
    public void finishSampleProcessing() throws IOException {
        for (int electrode = 1; electrode < nElectrodes; electrode++) {
            if (eventsGenerated[electrode] < minimumNoiseEvents) {
                System.out.println();
                System.out.println("Warning:  Only " + eventsGenerated[electrode] + " noise events generated on electrode " + electrode + ".");
            }
        }
        
        for (int listener = 0; listener < nListeners; listener++) {
            noiseListeners[listener].finishSpikeProcessing();
        }
    }


}
