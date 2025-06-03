package edu.ucsc.neurobiology.vision.io;

import java.io.*;


/**
 * This class can receive spikes from the SpikeFinder or the SpikeBuffer and saves them
 * to the spike file. Only the time and the amplitude of each spike are saved.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class SpikeSaver
    implements SpikeListener {

    private SpikeFile spikeFile;
    private final int nElectrodes;
    private int spikesFound;
    private int[] nSpikes;


    public SpikeSaver(SpikeFile spikeFile) throws IOException {
        this.nElectrodes = spikeFile.getNumberOfElectrodes();
        this.nSpikes = new int[nElectrodes];
        this.spikeFile = spikeFile;
    }


    public int getSpikesFound() {
        return spikesFound;
    }


    int oldTime = -1, oldEl = -10;

    public void processSpike(Spike spike) {
        // check order
        if (spike.time < oldTime) {
            System.out.println(
                "sorting error, old: " + oldTime + "/" + oldEl + ", current: " +
                spike.time + "/" + spike.electrode + " (" + spikesFound + ")");
        }
        oldTime = spike.time;
        oldEl = spike.electrode;

        try { //BEU
            spikeFile.addSpike( (short) spike.electrode, spike.time);
        } catch (IOException e) {
            e.printStackTrace();
        }

        nSpikes[spike.electrode]++;
        spikesFound++;
    }


    boolean first = true;


    synchronized public void finishSpikeProcessing() throws IOException {
        spikeFile.close();
    }

}
