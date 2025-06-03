package edu.ucsc.neurobiology.vision.math;

/**
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public abstract class RandomAbstract {

    // The starting seed for this sequence of random numbers.
    protected long seed;

    public RandomAbstract() {}

    public abstract short nextShort();

    public abstract char nextChar();

    public abstract int nextBits(int bits);

    /**
     * Set up the random number generator given a initial seed.
     * Resets generator to initial conditions.
     */
    public abstract void initialize(long seed);

    /**
     * Set the seed value for this random number generator.
     */
    public abstract void setSeed(long seed);

    /**
     * Returns the current seed of the generator.
     */
    public final long getSeed() {
        return seed;
    }

    /**
     * Simply advance the seed to the next value the given number of times
     */
    public abstract void advanceSeed(int numToAdvance);
}