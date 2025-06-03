package edu.ucsc.neurobiology.vision.io;

import java.io.*;


/**
 * This interface has to be implemented by any classes that wants to work with raw data.
 * The <tt>MultipleCompressedSampleInputStream</tt> class (and others) registers a set
 * of <tt>SampleListener</tt> classes and sends them the next data sample when it is read in.

 <br>
 <br>
 <b>The listeners should not modify the short[] containing the sample !!!</b>
 <br>
 <br>
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public interface SampleListener {

    /**
     * Called by raw data reading class to announce the listener about the new sample.
     */
    public void processSample(short[] sample);


    /**
     * Called by raw data reading class to announce the end of data.
     */
    public void finishSampleProcessing() throws IOException;
}
