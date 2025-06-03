package edu.ucsc.neurobiology.vision.stimulus;

import java.io.*;

import edu.ucsc.neurobiology.vision.math.RandomAbstract;


/**
 * This replicates the original Macintosh Toolbox Random() routine.  The
 * function is defined here explicitly for portability and independence from
 * changes in the MacOS.  EJC 1999-12-22
 * 
 * Ported to Java. CAL 2001-04-18
 *
 * @author Charles A. Loomis, University of California, Santa Cruz
 *         Dumitru Petrusca, University of California, Santa Cruz
 */
public class BinaryFrameGenerator extends AbstractConstantSeedFrameGenerator {
    public static final long BITS_PER_ELEMENT = 64 - 1;
    public static final long BITS_PER_COLOR = 7;
    public static final long PIXELS_PER_ELEMENT = BITS_PER_ELEMENT / (BITS_PER_COLOR*3);

    // The color mask is used to mask the lower bits of a long which correspond to a pixels color.
    public static final long COLOR_MASK = (long) Math.pow(2, BITS_PER_COLOR) - 1;
    public static final long LEFT_SHIFT = BITS_PER_COLOR * (PIXELS_PER_ELEMENT * 3 - 1);

    public static final long GREY_LONG = 1;

    // This is the default color lookup table.
    private final static int[] defaultColorLookup = {
         1, 1, 1,    1, 1, -1,    1, -1, 1,    1, -1, -1,
        -1, 1, 1,   -1, 1, -1,   -1, -1, 1,   -1, -1, -1
    };

    // must be defined as long[]
    private final static long[] defaultColorLookupInt = {
        2, 2, 2,   2, 2, 0,   2, 0, 2,   2, 0, 0,
        0, 2, 2,   0, 2, 0,   0, 0, 2,   0, 0, 0
    };

    // This is the default color lookup table for binary Independent (RGB) Separated stimulus.
    private final static int[] defaultColorLookupSeparated = {
        1, 0, 0,   -1, 0, 0,   0, 1, 0,   0, -1, 0,   0, 0, 1,   0, 0, -1
    };


    // must be declared as long[]
    private final static long[] defaultColorLookupSeparatedInt = {
        2, 1, 1,   0, 1, 1,   1, 2, 1,   1, 0, 1,   1, 1, 2,   1, 1, 0
    };
    
    
    // The color lookup table which is actually used to generate a new random
    // color.  The table contains positive or negative ones indicating whether
    // the given color component should be increased or decreased by delta
    // from a neutral value of 0.5.  Consequently, this table should have a
    // length which is a factor of three.
    //
    // For efficiency we mask off the necessary number of bits to get the index
    // of the RGB triplet.  Therefore, the total number of colors should be a
    // factor of two.
    //
    // Overall, the size should be 3*2^n where n is the number of RGB triplets.
    private float[] colorLookup;
    private long[] colorLookupInt;

    // The number of bits required for color random generated numbers.
    private int colorBits;

    // The size of the Lookup Table
    private int lutSize;

    // The size of the lookup table used
    private int lookupsize;

    // The delta by which to offset the color components from 0.5.
    private float colorContrast;
    public final int packedFrameLength;
        
    /**
     * Create a random number generator with the given seed.
     */
    public BinaryFrameGenerator(int width, int height, float colorContrast,
                                RandomNumberGenerator rng, long seed, ColorType colorType) {
        super(width, height, rng, seed, colorType);
        packedFrameLength = (int) Math.ceil(width * height / ((double) PIXELS_PER_ELEMENT));
        init(colorContrast);
    }
    
    public BinaryFrameGenerator(int width, int height, float colorContrast, RandomAbstract random, ColorType colorType) {
        super(width, height, random, colorType);
        packedFrameLength = (int) Math.ceil(width * height / ((double) PIXELS_PER_ELEMENT));
        init(colorContrast);
    }
    
    private void init(float colorContrast) {
        this.colorContrast = colorContrast;
        initColorLookup();
        initPixelGenerator();
    }
    
    private void initColorLookup() {
        lookupsize = defaultColorLookup.length;
        if (colorType == ColorType.SEPARATED) lookupsize = defaultColorLookupSeparated.length;

        float[] lookup    = new float[lookupsize];
        long[]  lookupInt = new long[lookupsize];

        for (int i = 0; i < lookup.length; i++) {
            if (colorType != ColorType.SEPARATED) {
                lookup[i] = FLOAT_GREY + colorContrast * defaultColorLookup[i];
                lookupInt[i] = defaultColorLookupInt[i];
            } else {
                lookup[i] = FLOAT_GREY + colorContrast * defaultColorLookupSeparated[i];
                lookupInt[i] = defaultColorLookupSeparatedInt[i];
            }
        }
        setColorLookup(lookup, lookupInt);
    }
    
    
    @Override
    public long[] createFrameBuffer() {
        return new long[packedFrameLength];
    }


    /**
     * A utility routine which sets the color lookup table to use.  Note that
     * the table is not copied, so changes to the table from elsewhere can
     * have unintended side effects.
     */
    private void setColorLookup(float[] lookup, long[] lookupInt) {
        // First, check that size is a multiple of three.
        int size = lookup.length;
        if (size % 3 != 0) 
            throw new IllegalArgumentException("Lookup size must be multiple of 3.");
 
        // Check that size is a factor of two.
        int nbits = 0;
        int temp = size / 3;
        for (int i = 0; i < 32; i++) {
            if ( (temp & 0x1) != 0) {
                nbits++;
            }
            temp >>>= 1;
        }

        // Everything is OK so find the number of bits required to calculate one color stixel.
        if (nbits != 1) {
            colorBits = 16;
        } else {
            colorBits = 0;
            temp = size / 3 - 1;
            while (temp > 0) {
                temp = temp >>> 1;
                colorBits++;
            }
        }
//        System.out.println("colorBits:" + colorBits);

        lutSize = size / 3;

        colorLookup = lookup;
        colorLookupInt = lookupInt;
    }
    
    
    private void initPixelGenerator() {
        switch (colorType) {
            case DEPENDENT:
                pixelGenerator = new PixelGenerator() {
                    public void setPixel(int pixelIndex, long[] frame) {
                        int longIndex = (int) (pixelIndex / PIXELS_PER_ELEMENT);
                        long intensity = 1 - (random.nextBits(1));
                        intensity *= 2;
                        for (int c = 0; c < 3; c++) {
                            frame[longIndex] >>= BITS_PER_COLOR;
                            frame[longIndex] += (intensity << LEFT_SHIFT);
                        }
                    }

                    public void setPixel(int pixelIndex, float[] frame) {
                        int dataindex = pixelIndex * 3;
                        int v = 1 - ( (random.nextBits(1)) << 1);
                        frame[dataindex++] = frame[dataindex++] = frame[dataindex++] = FLOAT_GREY + colorContrast * v;
                    }
                };
                break;
                
            case INDEPENDENT:
                pixelGenerator = new PixelGenerator() {
                    public void setPixel(int pixelIndex, long[] frame) {
                        int longIndex = (int) (pixelIndex / PIXELS_PER_ELEMENT);
                        int position = 3 * (random.nextBits(colorBits) % lutSize);
                        for (int c = 0; c < 3; c++) {
                            frame[longIndex] >>= BITS_PER_COLOR;
                            frame[longIndex] +=	(colorLookupInt[position++] << LEFT_SHIFT);
                        }
                    }

                    public void setPixel(int pixelIndex, float[] frame) {
                        int dataindex = pixelIndex * 3;
                        int lutindex = 3 * (random.nextBits(colorBits) % lutSize);
                        frame[dataindex++] = colorLookup[lutindex++];
                        frame[dataindex++] = colorLookup[lutindex++];
                        frame[dataindex]   = colorLookup[lutindex];
                    }
                    
                };
                break;
        }
    }
    
    
    @Override
    public void nextFrame(Object frame) throws IOException {
        if (frame instanceof float[]) {
            nextFrame((float[]) frame);
        } else if (frame instanceof long[]) {
            nextFrame((long[]) frame);
        }
    }


    @Override
    public void nextFrame(long seed, Object frame) throws IOException {
        random.setSeed(seed);
        nextFrame(frame);
    }
    
    
    public char[] nextMGLFrame(long seed) {
        random.setSeed(seed);
        return nextMGLFrame();
    }

    /**
     * Get frame in MGL format; RGBAxMxN uint8
     * @return
     */
    public char[] nextMGLFrame() {
        char[] frame = new char[4*nPixels];
        
        float[] dummy = new float[3];
        for (int i = 0; i < nPixels; i++) {
            pixelGenerator.setPixel(0, dummy);
            int mglIndex = 4*i;
            frame[mglIndex++] = (char) Math.round(dummy[0] * 255);
            frame[mglIndex++] = (char) Math.round(dummy[1] * 255);
            frame[mglIndex++] = (char) Math.round(dummy[2] * 255);
            frame[mglIndex] = 255;
        }
        
        return frame;
    }
    
    
    private void nextFrame(long[] frame) {
        for (int pixelIndex = 0; pixelIndex < nPixels; pixelIndex++)
            pixelGenerator.setPixel(pixelIndex, frame);
        finishFrame(frame);
    }
    
    @Override
    /**
     * If the nPixels in the frame is not a multiple of PIXELS_PER_ELEMENT, we have to shift the values over
     * to get them in the right final position.
     * @param frame
     */
    void finishFrame(long[] frame) {
        frame[frame.length - 1] >>= (frame.length * PIXELS_PER_ELEMENT - nPixels) * 3 * BITS_PER_COLOR;    	
    }
    
    private final void nextFrame(float[] frame) {
        for (int i = 0; i < nPixels; i++)
            pixelGenerator.setPixel(i, frame);
    }


    // @Override only okay on Java >= 6
    public void setGreyPixel(int pixelIndex, long[] frame) {
        int longIndex = (int) (pixelIndex / PIXELS_PER_ELEMENT);
        for (int c = 0; c < 3; c++) {
            frame[longIndex] >>= BITS_PER_COLOR;
        frame[longIndex] += (GREY_LONG << LEFT_SHIFT);
        }
    }


    @Override
    public int seedsPerPixel(Object frame) {
        return 1;
    }
    
}