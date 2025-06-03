package edu.ucsc.neurobiology.vision.io;

import java.io.*;


/**
 * This interface has to be implemented by any classes that wants to receive spikes.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public interface SpikeListener {

    /**
     * Called by the spike producer whenever a new spike is available.
     */
    public void processSpike(Spike spike);


    /**
     * Called by the spike producer to announce the end of the spike stream.
     */
    public void finishSpikeProcessing() throws IOException;

}
