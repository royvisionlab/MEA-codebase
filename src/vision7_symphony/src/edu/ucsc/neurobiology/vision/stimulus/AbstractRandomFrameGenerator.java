package edu.ucsc.neurobiology.vision.stimulus;

import java.io.IOException;

import edu.ucsc.neurobiology.vision.math.RandomAbstract;

/**
 * Defines shared functionality for FrameGenerators that are based on AbstractRandom seeded 
 * random number generators.  Actually, significant sections of FrameGenerator should be 
 * moved here eventually; most of the subclasses of FrameGenerator no longer use the 
 * AbstractRandom generator at all.
 * 
 * @author Peter H. Li, The Salk Institute
 */
public abstract class AbstractRandomFrameGenerator extends FrameGenerator {
    final ColorType colorType;
    
    public AbstractRandomFrameGenerator(int width, int height, RandomNumberGenerator rng, long seed, ColorType colorType) {
        super(width, height, rng, seed);
        this.colorType = colorType;
    }
    
    public AbstractRandomFrameGenerator(int width, int height, RandomAbstract random, ColorType colorType) {
        super(width, height, random);
        this.colorType = colorType;
    }
    
    /**
     * This advances the seed the given number of pixels.
     * @param numOffset
     */
    public abstract void offsetSeed(int numOffset, Object frame);
    
    public abstract void advanceFrameSeed(Object frame);
    
    /**
     * This version of nextFrame is used in random number frame generators that need
     * to know a seed and can be called out of order.
     *
     * @param seed long
     * @param frame Object
     * @throws IOException
     */
    public abstract void nextFrame(long seed, Object frame) throws IOException;
}
