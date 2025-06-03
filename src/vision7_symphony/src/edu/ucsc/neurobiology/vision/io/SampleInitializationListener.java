package edu.ucsc.neurobiology.vision.io;


/**
 * This interface must be implemnented to allow the implementing class to
 * receive a chuck of data from the beginning of the stream before any data gets
 * sent through the SampleListener interface.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public interface SampleInitializationListener {

    public void initialize(short[][] buffer);

}
