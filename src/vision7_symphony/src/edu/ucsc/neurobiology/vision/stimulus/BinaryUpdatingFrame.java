package edu.ucsc.neurobiology.vision.stimulus;

import java.util.*;

import edu.ucsc.neurobiology.vision.stimulus.FrameGenerator.ColorType;


/**
 * This class is a single frame of a binary RGB checkerboard that can be updated with
 * another frame. The frame uses PACKED_ENCODONG.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class BinaryUpdatingFrame
    extends UpdatingFrame {

    private int[] nPlus;
    private long[] packedNPlus;
    private int nCombinations, oldNCombinations;
    private final int packedFrameLength;
    private long countLimit;
    private final ColorType colorType;

    public BinaryUpdatingFrame(int width, int height, ColorType colorType) {
        super(width, height);
        this.colorType = colorType;

        nPlus = new int[width * height * 3];

        packedFrameLength = (int) Math.ceil(
            width * height / ( (double) BinaryFrameGenerator.PIXELS_PER_ELEMENT));
        packedNPlus = new long[packedFrameLength];

        // Set the maximum number of additions before the frame array is copied and cleaned
        countLimit = (BinaryFrameGenerator.COLOR_MASK / 2) - 1;
    }


    public final void update(Object frame) {
        long[] newPackedFrame = (long[]) frame;

        for (int i = 0; i < packedFrameLength; i++)
            packedNPlus[i] += newPackedFrame[i];

        nCombinations++;
        if (nCombinations - oldNCombinations >= countLimit) copy();
    }


    public void finish() {
        copy();
    }


    private final void copy() {
        if (nCombinations == oldNCombinations) return;
        
        int nPixels = width * height, index = 0, longIndex = -1;
        long value = 0;

        for (int pixelIndex = 0; pixelIndex < nPixels; pixelIndex++) {
            if (pixelIndex % BinaryFrameGenerator.PIXELS_PER_ELEMENT == 0) {
                longIndex++;
                value = packedNPlus[longIndex];
            }

            // red
            nPlus[index++] += (int) (value & BinaryFrameGenerator.COLOR_MASK);
            value >>= BinaryFrameGenerator.BITS_PER_COLOR;

            // green
            nPlus[index++] += (int) (value & BinaryFrameGenerator.COLOR_MASK);
            value >>= BinaryFrameGenerator.BITS_PER_COLOR;

            // blue
            nPlus[index++] += (int) (value & BinaryFrameGenerator.COLOR_MASK);
            value >>= BinaryFrameGenerator.BITS_PER_COLOR;
        }

        oldNCombinations = nCombinations;
        Arrays.fill(packedNPlus, 0);
    }


    public float[] getColorBuffer() {
        float[] buffer = new float[width * height * 3];
        float norm = 1.f / ( (float) nCombinations);
        
        int index = 0;
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                for (int c = 0; c < 3; c++) {
                    double p = nPlus[index] * norm / 2.0;
                    
                    // Normalize, assuming hardcoded stimulus range from 0.02 to 0.98.
                    buffer[index] = (float) (p * 0.98 + (1 - p) * 0.02 - 0.5);
                    
                    index++;
                }
            }
        }
        
        return buffer;
    }


    public float[] getErrorBuffer() {
        float[] buffer = new float[width * height * 3];
        float norm = 1.f / ( (float) nCombinations);
        
        int index = 0;
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                for (int c = 0; c < 3; c++) {
                    double p = nPlus[index] * norm / 2.0;
                    
                    // Assuming binomial distribution, calculate mean variance for the observed value
                    buffer[index] = (float) (0.96 * Math.sqrt(p * (1 - p) / nCombinations));
                    if (isSeparated()) buffer[index] = buffer[index] * 0.666f;
                    
                    index++;
                }
            }
        }
        
        return buffer;
    }

    
    public boolean isSeparated() {
        return colorType == FrameGenerator.ColorType.SEPARATED;
    }
    
    //Variance buffer is meaningless for binary movies.
    public float[] getCovarianceBuffer() {
        return new float[width * height * 3];
    }

}
