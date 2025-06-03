package edu.ucsc.neurobiology.vision.io;


/**
 * An iterator that iterates through spikes. An instance can be obtained form the
 * SpikeFile class.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public interface SpikeIterator {

    public Spike next();

    public boolean hasNext();
    
    public Spike current();
}