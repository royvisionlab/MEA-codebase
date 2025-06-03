package edu.ucsc.neurobiology.vision.stimulus;

import java.util.*;


/**
 * This class is a single frame of a Gaussian RGB checkerboard that can be updated with
 * another frame. The frame uses PACKED_ENCODING.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class GaussUpdatingFrame extends UpdatingFrame {
    private long[][] packedBuffers; // 0 is mean, 1 is var
    private int[][] buffers; 		// 0 is mean, 1 is var
    private int nCombinations, oldNCombinations;
    private final int packedFrameLength;
    private final boolean calculateSTV;
    private final int combinationRate;

    public GaussUpdatingFrame(int width, int height, boolean calculateSTV) {    	
        super(width, height);
        
        this.calculateSTV = calculateSTV;
        this.combinationRate = calculateSTV ? 31 : 8191; // 2^13 = 8192
        buffers = new int[2][width * height * 3];
        packedFrameLength = width * height;
        packedBuffers = new long[2][packedFrameLength];
    }
    
    
    public final void update(Object frame) {
        long[] newPackedFrame = (long[]) frame;
        
        for (int i = 0; i < packedFrameLength; i++)
            packedBuffers[0][i] += newPackedFrame[i];
        
        if (calculateSTV) {
            for (int i = 0; i < packedFrameLength; i++) {
                long value = newPackedFrame[i];
                long r = value & GaussFrameGenerator.COLOR_MASK;
                value >>= GaussFrameGenerator.BITS_PER_COLOR;
                long g = value & GaussFrameGenerator.COLOR_MASK;
                value >>= GaussFrameGenerator.BITS_PER_COLOR;
                long b = value & GaussFrameGenerator.COLOR_MASK;
                
                packedBuffers[1][i] += ((b*b) << GaussFrameGenerator.TWICE_BITS_PER_COLOR) +
                                       ((g*g) << GaussFrameGenerator.BITS_PER_COLOR) +
                                        (r*r);
            }
        }
        
        nCombinations++;
        if (nCombinations - oldNCombinations >= combinationRate) copy();
    }
    
    
    public void finish() {
        copy();
    }
    
    
    private final void copy() {
        int nPixels = width * height;
        long value;
        int nBuffers = calculateSTV ? 2 : 1;
        for (int i = 0; i < nBuffers; i++) {
            for (int pixelIndex = 0, index = 0; pixelIndex < nPixels; pixelIndex++) {
                value = packedBuffers[i][pixelIndex];
                
                // red
                buffers[i][index++] += (int) (value & GaussFrameGenerator.COLOR_MASK);
                value >>= GaussFrameGenerator.BITS_PER_COLOR;
                
                // green
                buffers[i][index++] += (int) (value & GaussFrameGenerator.COLOR_MASK);
                value >>= GaussFrameGenerator.BITS_PER_COLOR;
                
                // blue
                buffers[i][index++] += (int) (value & GaussFrameGenerator.COLOR_MASK);
                value >>= GaussFrameGenerator.BITS_PER_COLOR;
            }
            
            Arrays.fill(packedBuffers[i], 0);
        }
        oldNCombinations = nCombinations;
    }
    
    
    public float[] getColorBuffer() {
        float[] buffer = new float[width * height * 3];
        final float norm = 1.f / ( (float) nCombinations);
        int index = 0;
        
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                for (int c = 0; c < 3; c++) {
                    // Normalize, assuming hardcoded stimulus range from 0.02 to 0.98.
                    buffer[index] = buffers[0][index] * norm / 255.0f - 0.5f;

                    index++;
                }
            }
        }
        
        return buffer;
    }
    
    
    public float[] getErrorBuffer() {
        float[] buffer = new float[width * height * 3];
        double N = nCombinations;
        
        for (int i = 0, k = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                for (int c = 0; c < 3; c++) {
//                    double m = buffers[k] / N;
//                    buffer[k] = (float) Math.sqrt( (sigma[k] - N * m * m) / (N - 1))
//                                / 255.f / 255.f;
                    
                    buffer[k] = (float) (0.16 / Math.sqrt(N));
                    
//                    double m = buffers[k] / N;
//                    buffer[k] = (float) Math.sqrt( (sigma[k] - N * m * m) / (N - 1))
//                                / (255.f * (float)Math.sqrt(N));
                    
                    k++;
                }
            }
        }
        
        return buffer;
    }


    public float[] getCovarianceBuffer() {
        float[] buffer = new float[width * height * 3];
        float N = nCombinations;

        for (int i = 0, k = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                for (int c = 0; c < 3; c++) {
                    // (X + 1) / 257 is the correct conversion between short
                    // and float (from ej)
                    // This leads to the values being symmetric about .5
                    float meanFloat = (buffers[0][k] + N) / 257f;

                    // SUM((x+1)(x+1)/257/257) = SUM((x^2 + 2x + 1) / 257 / 257)
                    float sigmaFloat = (buffers[1][k] + 2 * buffers[0][k] + N) / 257f / 257f;

                    //calculate spike triggered variance, given SUM(x^2) and SUM(x)
                    //The buffers used is the buffers from the stimulus: .5, not the buffers of
                    //the response
                    //buffer[k] = (1/(N-1)) SUM (x_i - .5)^2
                    buffer[k] = (float)
                                (sigmaFloat - meanFloat + 0.5f * 0.5f * N) / (N - 1f);
                    k++;
                }
            }
        }

        return buffer;
    }

}
