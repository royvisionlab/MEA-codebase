package edu.ucsc.neurobiology.vision.stimulus;

import java.io.*;

import edu.ucsc.neurobiology.vision.math.*;


/**
 * Defines the shared functionality of a FrameGenerator. Such generators are used to
 * generate movie frames aw well as in STA calculation.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public abstract class FrameGenerator {
    public static enum ColorType {
        INDEPENDENT, DEPENDENT, SEPARATED
    }

    public static enum RandomNumberGenerator {
        MAC_RANDOM, JAVA_RANDOM, JAVA_RANDOM_V2,
    }


    protected RandomAbstract random;
    protected final int width, height;
    protected final int nPixels, nColors;


    public FrameGenerator(int width, int height, RandomNumberGenerator rng, long seed) {
        this.width = width;
        this.height = height;
        this.nPixels = width * height;
        this.nColors = nPixels * 3;
        
        if (rng == RandomNumberGenerator.MAC_RANDOM) {
            random = new RandomMacToolbox();
        } else if (rng == RandomNumberGenerator.JAVA_RANDOM) {
            random = new RandomJava();
        } else {
            random = new RandomJavaV2();
        }
        random.initialize(seed);
    }
    
    public FrameGenerator(int width, int height, RandomAbstract random) {
        this.width = width;
        this.height = height;
        this.nPixels = width * height;
        this.nColors = nPixels * 3;

        this.random = random;
    }

   
    public long getSeed() {
        return random.getSeed();
    }


    public int getNumberOfFrames() {
        return -1;
    }


    /**
     * This method is called one for every frame.
     *
     * @param frame Object
     * @throws IOException
     */
    public abstract void nextFrame(Object frame) throws IOException;


    /**
     * Must create a frame buffer (not a frame!) compatible with this frame.
     * @return Object
     */
    public abstract Object createFrameBuffer();
    
    public float[] nextFloatFrame() throws IOException {
        float[] frame = new float[nColors];
        nextFrame(frame);
        return frame;
    }
}
