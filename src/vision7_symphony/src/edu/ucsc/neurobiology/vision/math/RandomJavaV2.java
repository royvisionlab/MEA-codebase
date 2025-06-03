package edu.ucsc.neurobiology.vision.math;

/**
 This class was copied from standard java source, java.lang.Object,
 version 1.4.2, to insure future compatibility.  Small changes were
 made to allow natural integration with vision.
 
 The only change in V2 is to make the nextBits function take the most significant
 bit instead of some middle bits from the seed.  Otherwise, the frames are
 highly correlated to each other.
 *
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class RandomJavaV2 extends RandomAbstract {
//    private boolean haveNextNextGaussian = false;


    public RandomJavaV2() {
        this(System.currentTimeMillis());
    }


    public RandomJavaV2(long seed) {
        initialize(seed);
    }


    synchronized public void initialize(long seed) {
        super.seed = (seed ^ 0x5DEECE66DL) & ( (1L << 48) - 1);
//        haveNextNextGaussian = false;
    }


    synchronized public void setSeed(long seed) {
        super.seed = seed;
//        haveNextNextGaussian = false;
    }
    
    
    synchronized protected int next(int bits) {
        super.seed = (seed * 0x5DEECE66DL + 0xBL) & ( (1L << 48) - 1);
        return (int) (seed >>> (48 - bits));
    }
    
    
    public int nextInt() {
        return next(32);
    }
    
    
    public int nextInt(int n) {
        if (n <= 0) {
            throw new IllegalArgumentException("n must be positive");
        }
        
        if ( (n & -n) == n) { // i.e., n is a power of 2
            return (int) ( (n * (long) next(31)) >> 31);
        }
        
        int bits, val;
        do {
            bits = next(31);
            val = bits % n;
        } while (bits - val + (n - 1) < 0);
        return val;
    }


    public long nextLong() {
        return ( (long) next(32) << 32) + next(32);
    }


    public boolean nextBoolean() {
        return next(1) != 0;
    }


    public float nextFloat() {
        return next(24) / ( (float) (1 << 24));
    }


    @Override
    public short nextShort() {
        return (short) (this.nextInt() >> 16);
    }

    @Override
    public char nextChar() {
        return (char) (this.nextInt() >> 16);
    }

    /**
     *
     * Return a random number, that uses bit shifting to give the correct number of
     * bits.
     */
    @Override
    public int nextBits(int bits) {
        return next(bits);
    }

    
    /**
     * Just advances the seed to its next value the given number of times
     */
    @Override
    synchronized public void advanceSeed(int numToAdvance) {
        long seed = super.seed;
        long mask = ( (1L << 48) - 1);
        for (int i = 0; i < numToAdvance; i++) {
            seed = (seed * 0x5DEECE66DL + 0xBL) & mask;
        }
        super.seed = seed;
    }
}
