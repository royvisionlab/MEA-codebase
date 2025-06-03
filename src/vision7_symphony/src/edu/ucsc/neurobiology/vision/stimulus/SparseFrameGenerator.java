package edu.ucsc.neurobiology.vision.stimulus;

import java.io.IOException;

public class SparseFrameGenerator extends AbstractRandomFrameGenerator {
    private static long charRange = Short.MAX_VALUE - Short.MIN_VALUE;
    
    private float probability;
//	private char threshold;
  
    // Somewhat ugly kludge.  We want the pixel generation code to be DRY; we shouldn't have a copied and pasted 
    // algorithm for getting Binary vs. Gaussian pixels.  But the underlying FrameGenerators (Binary/Gaussian) are 
    // not abstracted/modular enough to do this nicely.  So for now, we create an underlying Binary/Gaussian 
    // frameGenerator to handle most of the low level work.  The underlying generator shares our random number 
    // generator, so when we advance the seed here it will reflect this.
    private AbstractConstantSeedFrameGenerator underlyingFrame;
        
    /**
     * Create a new SparseFrameGenerator based on an underlyingFrame, which we keep around to help generate pixel 
     * values.  We share the same seeded random number generator object with the underlyingFrame.
     * 
     * @param probability
     * @param underlyingFrame
     */
    public SparseFrameGenerator(float probability, AbstractConstantSeedFrameGenerator underlyingFrame) {
        super(underlyingFrame.width, underlyingFrame.height, underlyingFrame.random, underlyingFrame.colorType);
        this.underlyingFrame = underlyingFrame;
        setProbability(probability);
    }
    
    
    private void setProbability(float probability) {
        this.probability = probability;
//		long scaledCharRange = (long) Math.round(this.probability * charRange);
//		this.threshold = (char) (scaledCharRange + Character.MIN_VALUE);
    }
    

    /**
     * For some reason it seems I wrote this to return true if the next pixel is NOT grey?
     * @return
     */
    private boolean nextPixelGrey() {
//		return threshold > random.nextChar();
        return probability > (float) (random.nextChar() - Character.MIN_VALUE) / (float) (charRange);
    }
    
    
    @Override
    public void offsetSeed(int numOffset, Object frame) {
        for (int i = 0; i < numOffset; i++) {
            if (nextPixelGrey())
                random.advanceSeed(underlyingFrame.seedsPerPixel(frame));
        }
    }
    
    
    @Override
    public void advanceFrameSeed(Object frame) {
        offsetSeed(nPixels, frame);
    }

    
    @Override
    public void nextFrame(long seed, Object frame) throws IOException {
        random.setSeed(seed);
        nextFrame(frame);
    }
    
    @Override
    public void nextFrame(Object frame) {
        if (frame instanceof float[]) {
            nextFrame((float[]) frame);
        } else if (frame instanceof long[]) {
            nextFrame((long[]) frame);
        }
    }
    
    
    private void nextFrame(float[] frame) {
        for (int pixelIndex = 0; pixelIndex < nPixels; pixelIndex++) {
            if (nextPixelGrey())
                underlyingFrame.pixelGenerator.setPixel(pixelIndex, frame);
            else setGreyPixel(pixelIndex, frame);
        }
    }

    private void setGreyPixel(int pixelIndex, float[] frame) {
        int dataindex = pixelIndex * 3;
        frame[dataindex++] = frame[dataindex++] = frame[dataindex++] = AbstractConstantSeedFrameGenerator.FLOAT_GREY;
    }

    private void nextFrame(long[] frame) {
        for (int pixelIndex = 0; pixelIndex < nPixels; pixelIndex++) {
            if (nextPixelGrey()) {
                underlyingFrame.pixelGenerator.setPixel(pixelIndex, frame);
            } else {
                underlyingFrame.setGreyPixel(pixelIndex, frame);
            }	
        }
        underlyingFrame.finishFrame(frame);
    }

    
    @Override
    public Object createFrameBuffer() {
        return underlyingFrame.createFrameBuffer();
    }
    
    protected interface PackedGreyPixelSetter {
        void setGreyPixel(int pixelIndex, long[] frame);
    }
}
