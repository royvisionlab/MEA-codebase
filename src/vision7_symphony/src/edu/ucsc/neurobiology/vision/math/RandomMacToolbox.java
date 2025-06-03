package edu.ucsc.neurobiology.vision.math;


/**
 * This replicates the original Macintosh Toolbox Random() routine.  The
 * functon is defined here explicitly for portability and independence from
 * changes in the MacOS.  EJC 1999-12-22
 * Ported to Java. CAL 2001-04-18
 *
 * @author Charles A. Loomis, University of California, Santa Cruz
 */
public class RandomMacToolbox
    extends RandomAbstract {
    private final static short NEG_ZERO = (short) 0x8000;


    public RandomMacToolbox() {
        super();
    }


    public final void initialize(long seed) {
        super.seed = seed;
    }


    public final void setSeed(long seed) {
        super.seed = seed;
    }


    /**
     * Return the next signed short value in the sequence.
     */
    @Override
    public final short nextShort() {
        int temp1, temp2, temp3;
        temp1 = ( (int) (seed & 0xFFFFL)) * 0x41A7;
        temp2 = ( (int) (seed >> 16L)) * 0x41A7 + (temp1 >> 16);
        temp3 = temp2 * 2 >> 16;
        temp1 = (temp1 & 0xFFFF) - 0x7FFFFFFF;
        temp2 &= 0x7FFF;
        temp1 += (temp2 << 16) | (temp2 >> 16) + temp3;
        if (temp1 < 0) {
            temp1 += 0x7FFFFFFF;
        }
        seed = (long) temp1;
        short result = (short) (temp1 & 0xFFFF);
        return (result == NEG_ZERO) ? 0 : result;
    }
    
    @Override
    public final char nextChar() {
        return (char) nextShort();
    }


    /**
     * Return a random number, that uses masking to give the correct number of
     * bits.
     */
    public int nextBits(int bits) {
        return (int) nextShort() & ( (1 << bits) - 1);
    }

    
    /**
     * Just advances the seed to its next value
     */
    public void advanceSeed() {
        // FIXME: could be optimized
        int temp1, temp2, temp3;
        temp1 = ( (int) (seed & 0xFFFFL)) * 0x41A7;
        temp2 = ( (int) (seed >> 16L)) * 0x41A7 + (temp1 >> 16);
        temp3 = temp2 * 2 >> 16;
        temp1 = (temp1 & 0xFFFF) - 0x7FFFFFFF;
        temp2 &= 0x7FFF;
        temp1 += (temp2 << 16) | (temp2 >> 16) + temp3;
        if (temp1 < 0) {
            temp1 += 0x7FFFFFFF;
        }
        seed = (long) temp1;

    }
    
    /**
     * Advance the seed the given number of times
     */
    public void advanceSeed(int numToAdvance) {
        // FIXME: could be optimized
        for (int i = 0; i < numToAdvance; i++) {
            advanceSeed();
        }
    }
}
