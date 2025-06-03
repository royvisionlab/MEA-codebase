package edu.ucsc.neurobiology.vision.stimulus;

import java.io.*;

import edu.ucsc.neurobiology.vision.math.RandomAbstract;


/**
 * Most things (spatial and temporal parameters, etc) are the same as for binary.
 * The only significant difference is how the intensity values are drawn
 * for a given pixel in a given frame.
 * I generate Gaussian distributed intensities using the standard
 * method: use a random number generator that has a (approximately)
 * uniform distribution, and use this number to look up a value in a
 * one-dimensional array.  The array is constructed so that its values,
 * when plotted against the index of the array, take the form of the
 * inverse of the cumulative normal (or, rather, an approximation of
 * this).  Call this function F.  The approximation to the F is given
 * with C code below.  I'm sorry but I don't remember where I got this
 * code -- it was a fairly standard source (I don't think it was NR, but
 * it may have been a book by Abramowitz and Stegun).
 * The inverse cumulative normal table I use is of size 256.
 * Specifically, for each index i from 0 to 255, I compute evenly spaced
 * floating point numbers between 0 and 1 that are bin-centers, without
 * the endpoints, that is x = (i+1)/257.  The value y in the table at
 * the ith location is then given by F(x).
 * The random number generator I use to index into the table is the same
 * one that you already have.  The low order 8 bits obtained from the
 * random number are used to index into the table.
 * The only remaining issue is the standard deviation.  You will recall
 * we specify intensities as fractions of maximum gun intensity, as
 * deviations from a mean (typically, the mean is 0.5).  Given that the
 * gun intensities are limited between 0 and 1 (actually, a little
 * higher than zero, because the darkest you can go is not zero
 * intensity), one has to live with an approximation to a Gaussian
 * distribution.  We often use a "Gaussian" distribution with a standard
 * deviation of 0.16 (in our stimulus specification, this would be RGB =
 * [0.16, 0.16, 0.16]).  That is, we create the normal table as I
 * described above, then multiply the values in that table by 0.16.
 * This procedure yields a table with lowest value of -0.4258835 and
 * highest value 0.4258835 (use these values to check your code).  Add
 * this to 0.5 (the mean), and you get the intensity that we place on
 * the monitor.  Infrequently, when we need more contrast and are
 * willing to make a poorer approximation to the Gaussian, we use a
 * contrast of 0.24.  In this case the smallest number in the table is
 * -0.638825 and the largest is 0.638825.  These are not physically
 * realizable of course, and the distribution is clipped when the
 * stimulus is actually presented.  When we compute the STAs, we use the
 * values in the table, rather than the clipped displayed values.  This
 * is perhaps not ideal but should not matter much if the tables have
 * few values out of range.
 * As with binary, we either draw the guns independently of one another,
 * or yoked.  As with binary, in the yoked case, only one random number
 * is drawn rather than three, and it is used to index into the red,
 * green and blue gun tables.
 *
 * @author Sasha Sher, University of California, Santa Cruz
 *         Dumitru Petrusca, University of California, Santa Cruz
 */
public class GaussFrameGenerator
    extends AbstractConstantSeedFrameGenerator {

    public static final long BITS_PER_COLOR = 21;
    public static final long TWICE_BITS_PER_COLOR = BITS_PER_COLOR * 2;
    public static final long COLOR_MASK = (long) Math.pow(2, BITS_PER_COLOR) - 1;
    
    public static final int INT_GREY = Math.round(FLOAT_GREY * 255);
    public static final long LONG_GREY = (INT_GREY << TWICE_BITS_PER_COLOR) + 
                                         (INT_GREY << BITS_PER_COLOR) +
                                          INT_GREY;

    public static final long LEFT_SHIFT = BITS_PER_COLOR * 2;


    private float[][] lookup = new float[3][256];
    // it is very important for this table to be declared as long[][] !!!
    // Otherwise there are problems with the binary shifts
    private long[][] intLookup = new long[3][256];
    private float[] sigma;
    public final int packedFrameLength;

    /**
     * Create a gaussian random number generator with the given seed.
     * Sets up three lookup tables with values from the inverse cumulative
     * normal function (F): double lookup[256] = F
     * double lookup[256][i=0..2]=lookup[256]*sigma[i]+0.5
     * long intlookup[256][i=0..2]=lookup[256][i]*255
     * the last one is "clipped" at 0 and 255.
     * That is, all the values <0 or >255 become 0 or 255 respectively.
     */
    public GaussFrameGenerator(int width, int height, float[] sigma,
                               RandomNumberGenerator rng, long seed, ColorType colorType) {
        super(width, height, rng, seed, colorType);
        packedFrameLength = width * height;
        init(sigma);
    }
    
    public GaussFrameGenerator(int width, int height, float[] sigma, RandomAbstract random, ColorType colorType) {
        super(width, height, random, colorType);
        packedFrameLength = width * height;
        init(sigma);
    }
    
    private void init(float[] sigma) {
        this.sigma = sigma;
        
        float p;
        for (int i = 0; i < lookup[0].length; i++) {
            p = inverseCumulativeNormal( (i + 1.0) / 257.0);
            for (int c = 0; c < 3; c++) {
                lookup[c][i] = p * sigma[c] + FLOAT_GREY;
                intLookup[c][i] = (int) Math.round(lookup[c][i] * 255);
            }
        }
        
        initPixelGenerator();
    }
    
    private void initPixelGenerator() {
        switch (colorType) {
            case DEPENDENT:
                pixelGenerator = new PixelGenerator() {
                    public void setPixel(int pixelIndex, long[] frame) {
                        // It is very important for this variable to be declared as long, otherwise there are problems with the binary shifts
                        int intensity = random.nextBits(8);
                        long r = intLookup[0][intensity];
                        long g = intLookup[1][intensity];
                        long b = intLookup[2][intensity];
                        frame[pixelIndex] = (b << TWICE_BITS_PER_COLOR) + (g << BITS_PER_COLOR) + r;
                    }

                    public void setPixel(int pixelIndex, float[] frame) {
                        int dataindex = pixelIndex * 3;
                        int position = random.nextBits(8);
                        for (int c = 0; c < 3; c++)
                            frame[dataindex++] = lookup[0][position];
                    }
                };
                break;

            case INDEPENDENT:
                pixelGenerator = new PixelGenerator() {
                    public void setPixel(int pixelIndex, long[] frame) {
                        long r = intLookup[0][random.nextBits(8)];
                        long g = intLookup[1][random.nextBits(8)];
                        long b = intLookup[2][random.nextBits(8)];
                        frame[pixelIndex] = (b << TWICE_BITS_PER_COLOR) + (g << BITS_PER_COLOR) + r;
                    }

                    public void setPixel(int pixelIndex, float[] frame) {
                        int dataindex = pixelIndex * 3;
                        for (int c = 0; c < 3; c++) {
                            int position = random.nextBits(8);
                            frame[dataindex++] = lookup[c][position];
                        }
                    }    				
                };
                break;    	        
        }
    }


    public Object createFrameBuffer() {
        return new long[packedFrameLength];
    }


    @Override
    public void nextFrame(long seed, Object frame) throws IOException {
        random.setSeed(seed);
        nextFrame(frame);
    }

    public void nextFrame(Object frame) {
        if (frame instanceof float[]) {
            nextFrame((float[]) frame);
        } else if (frame instanceof long[]) {
            nextFrame((long[]) frame);
        }
    }


    public void nextFrame(float[] frame) {
        for (int i = 0; i < nPixels; i++) 
            pixelGenerator.setPixel(i, frame);
    }
    
    public void nextFrame(long[] frame) {
        for (int i = 0; i < nPixels; i++) 
            pixelGenerator.setPixel(i, frame);
    }
    
        
    /**
     * This currently doesn't work because GaussFrameGenerator.LONG_GREY does not give a true
     * unbiased grey.  This is because the values in the Gaussian int lookup table straddle the
     * true grey value.  Getting true unbiased middle grey values in Gaussian packed encoding 
     * will require changing the packed encoding or adding an additional structure to transmit
     * information about grey pixels between the GaussFrameGenerator and the GaussUpdatingFrame.
     * :(
     * 
     * If you want to get some STAs to look at, you can uncomment the line below and comment out
     * the error.  This should give you reasonable looking STAs, but they will be biased due to 
     * the grey value bias, so they cannot be used for higher order calculations. 
     * 
     * @Override only works on Java >=6
     *
     * @param pixelIndex
     * @param frame
     */
    public void setGreyPixel(int pixelIndex, long[] frame) {
//		frame[pixelIndex] = LONG_GREY;
        throw new Error("Sparse gaussian movies cannot currently get unbiased grey values.");
    }
    
    
    /**
     * Calculates the Inverse Of The Cumulative Normal Distribution Function
     * as per EJ's C code.
     */
    private float inverseCumulativeNormal(double p) {
        double tmp;

        if (p < 0 || p > 1) {
            tmp = -1;
        } else if (p > 0.5) {
            tmp = -inverseCumulativeNormal(1 - p);
        } else {
            double tmp_sqr = ( -2 * Math.log(p));
            tmp = (float) Math.sqrt(tmp_sqr);

            tmp = (2.515517 + 0.802853 * tmp + 0.010328 * tmp_sqr) /
                  (1.0 + 1.432788 * tmp + 0.189269 * tmp_sqr +
                   0.001308 * tmp * tmp_sqr) - tmp;
        }

        return (float) tmp;
    }

    
    public int seedsPerPixel(Object frame) {
        if (colorType == ColorType.INDEPENDENT) return 3;
        return 1;
    }

    
    /**
     * Doesn't need to do anything here; BinaryFrameGenerator does something.
     */
    @Override
    void finishFrame(long[] frame) {}
}
