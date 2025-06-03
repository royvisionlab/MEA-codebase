package edu.ucsc.neurobiology.vision.stimulus;

import edu.ucsc.neurobiology.vision.math.RandomAbstract;

public abstract class AbstractConstantSeedFrameGenerator 
        extends AbstractRandomFrameGenerator implements SparseFrameGenerator.PackedGreyPixelSetter {
    public static final float FLOAT_GREY = 0.5f;
    protected PixelGenerator pixelGenerator;
    
    public AbstractConstantSeedFrameGenerator(int width, int height, 
            RandomNumberGenerator rng, long seed, ColorType colorType) {
        super(width, height, rng, seed, colorType);
    }
    
    public AbstractConstantSeedFrameGenerator(int width, int height, RandomAbstract random, ColorType colorType) {
        super(width, height, random, colorType);
    }
    
    /**
     * This advances the seed the given number of pixels.
     * @param numOffset
     */
    @Override
    public void offsetSeed(int numOffset, Object frame) {
        random.advanceSeed(numOffset * seedsPerPixel(frame));
    }
    
    @Override
    public void advanceFrameSeed(Object frame) {
        random.advanceSeed(nPixels * seedsPerPixel(frame));
    }
    
    public abstract int seedsPerPixel(Object frame);
    
    void finishFrame(long[] frame) {}
    
    protected interface PixelGenerator {
        void setPixel(int pixelIndex, long[] frame);
        void setPixel(int pixelIndex, float[] frame);
    }
    
    protected interface MGLPixelGenerator {
        void setPixel(int pixelIndex, char[] frame);
    }
}
