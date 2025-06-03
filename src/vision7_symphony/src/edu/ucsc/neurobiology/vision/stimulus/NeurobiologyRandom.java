package edu.ucsc.neurobiology.vision.stimulus;


/**
 * This replicates the original Macintosh Toolbox Random() routine.  The
 * functon is defined here explicitly for portability and independence from
 * changes in the MacOS.  EJC 1999-12-22
 * Ported to Java. CAL 2001-04-18
 *
 * @author Charles A. Loomis, University of California, Santa Cruz
 */
public class NeurobiologyRandom {
    private final static short NEG_ZERO = (short) 0x8000;

    // The starting seed for this sequence of random numbers.
    private int seed;


    /**
     * Return the next signed short value in the sequence.
     */
    public final short nextShort() {
        int temp1, temp2, temp3;

        temp1 = (seed & 0xFFFF) * 0x41A7;
        temp2 = (seed >> 16) * 0x41A7 + (temp1 >> 16);
        temp3 = temp2 * 2 >> 16;
        temp1 = (temp1 & 0xFFFF) - 0x7FFFFFFF;
        temp2 &= 0x7FFF;
        temp1 += (temp2 << 16) | (temp2 >> 16) + temp3;
        if (temp1 < 0) {
            temp1 += 0x7FFFFFFF;
        }
        seed = temp1;

        short result = (short) (temp1 & 0xFFFF);
        return (result == NEG_ZERO) ? 0 : result;
    }


    /**
     * Set the seed value for this random number generator.
     */
    public void setSeed(int seed) {
        this.seed = seed;
    }


    /**
     * Returns the current seed of the generator.
     */
    public final int getSeed() {
        return seed;
    }

}
