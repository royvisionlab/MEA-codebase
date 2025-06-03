package edu.ucsc.neurobiology.vision.stimulus;

import java.util.*;


/**
 * This class defines how to create each frame of a pink movie, version 1.
 * The code currently does not call this function.  It is kept on the outside
 * chance that someone will want to analyze the early pink noise data.
 *
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class PinkFrameGeneratorV1
    extends FrameGenerator {

    private float[][] lookup = new float[3][256]; //lookup table for random number generator (float)
    private int[][] intLookup = new int[3][256]; //lookup table for random number generator (int)
    private float[] sigma; //standard deviation of final stimulus.
    private int spatialFreqDepth; //number of spatial frequencies
    private int temporalFreqDepth; //number of temporal frequencies
    private float[][] currentSpatialFrames; //spatial frames that are averaged together to create actual frame.
    private int[] framesTillSwitch; //keeps track of when each spatial frame switches again.
    private int[] frameSwitchFreq; //record of how frequently each spatial frame switches.
    private float normalizer; //used to normalize the frames to the correct standard deviation
    private float[] currentColors; //current frame
    private ColorType colorType;
    private MovieType noiseType; //RGB vs GRAY and BINARY_MOVIE vs GAUSSIAN_MOVIE
//    public static final int RGB = 0, GRAY = 2;

    static int bwColor = 0;
    int frameNumber = 0;


    public PinkFrameGeneratorV1(int width, int height, int spatialFreqDepth,
                                int temporalFreqDepth, float[] sigma, int startSeed,
                                ColorType colorType, MovieType noiseType) {

        super(width, height, RandomNumberGenerator.MAC_RANDOM, 0);

        this.sigma = sigma;
        this.colorType = colorType;
        this.noiseType = noiseType;

        /**
         * Create a gaussian random number generator with the given seed.
         * Sets up three lookup tables with values from the inverse cumulative
         * normal function (F): double lookup[256] = F
         * double lookup[256][i=0..2]=lookup[256]*sigma[i]+0.5
         * long intlookup[256][i=0..2]=lookup[256][i]*255
         * the last one is "clipped" at 0 and 255.
         * That is, all the values <0 or >255 become 0 or 255 respectively.
         */

        random.setSeed(startSeed);
        float p;
        for (int i = 0; i < lookup[0].length; i++) {
            //i+1/257 gives an even set of 256 numbers.  These map from 0 to 255 in the saved file.
            //This is ej's standard method.  It makes it so that all the bins are symmetric.
            p = inverseCumulativeNormal( (i + 1.0) / 257.0);
            for (int c = 0; c < 3; c++) {
//                lookup[c][i] = Math.abs(p * sigma[c]);
                lookup[c][i] = p * sigma[c];
                //             intLookup[c][i] = (int) Math.round(lookup[c][i] * 255);
            }
        }

        this.spatialFreqDepth = spatialFreqDepth;
        this.temporalFreqDepth = temporalFreqDepth;

        currentSpatialFrames = new float[temporalFreqDepth][width * height * 3];
        framesTillSwitch = new int[temporalFreqDepth];
        frameSwitchFreq = new int[temporalFreqDepth];
        normalizer = (float) Math.sqrt(spatialFreqDepth * temporalFreqDepth);
        currentColors = new float[width * height * 3];

//        Temporal Frames switch in the following way:
//         XXXXXXXXXXXXXXX
//         X X X X X X X X
//          X   X   X   X
//            X       X
//                X
//       //Create first frame, which is a sum of temporalFreqDepth spatial frames.
        // Set up times for each temporal frame to switch.
        int powerOfTwo = 1;
        for (int i = 0; i < temporalFreqDepth; i++) {

            currentSpatialFrames[i] = createSpatialFrame(currentSpatialFrames[i]);
            for (int j = 0; j < currentColors.length; j++) {
                currentColors[j] += currentSpatialFrames[i][j];
            }
            frameSwitchFreq[i] = powerOfTwo;
            if (i == 0) {
                framesTillSwitch[i] = 2;
            } else if (i == 1) {
                framesTillSwitch[i] = 3;
            } else {
                framesTillSwitch[i] = powerOfTwo / 2;
            }
            powerOfTwo *= 2;
        }

    }


    public Object createFrameBuffer() {
        return new int[width * height * 3];
    }


    public void nextFrame(Object frame) {
        if (frame instanceof float[]) {
            nextFrame( (float[]) frame);
        } else if (frame instanceof int[]) {
            throw new IllegalStateException(
                "Frames are not created as int[] with pink noise.");
        }
    }


    private final void nextFrame(float[] color) {
        //for each temporal frequency replace the spatial frame, if necessary
        for (int j = 0; j < temporalFreqDepth; j++) {
            framesTillSwitch[j]--;
            if (framesTillSwitch[j] == 0) {
                framesTillSwitch[j] = frameSwitchFreq[j];
                //subtract old spatial frame from mix
                for (int k = 0; k < currentColors.length; k++) {
                    currentColors[k] -= currentSpatialFrames[j][k];
                }

                currentSpatialFrames[j] = createSpatialFrame(currentSpatialFrames[j]);

                //add new spatial frame to mix
                for (int i = 0; i < currentColors.length; i++) {
                    currentColors[i] += currentSpatialFrames[j][i];
                }
            }
        }

        switch (colorType) {
            case INDEPENDENT:

                //put colors in correct form for display.
                for (int i = 0; i < currentColors.length; i++) {
//                        color[i] = 2*(currentColors[i]/normalizer)*(currentColors[i]/normalizer )+.5f;
                    color[i] = (currentColors[i] / normalizer) + .5f;
                    if (noiseType == MovieType.GAUSSIAN_MOVIE) {
                        if (color[i] < 1.0 / 257.0) {
                            color[i] = (float) (1.0 / 257.0);
                        } else if (color[i] > 256.0 / 257.0) {
                            color[i] = (float) (256.0 / 257.0);
                        }
                    } else {
                        if (color[i] < .5) {
                            color[i] = (float) .5 - sigma[i % 3];
                        } else {
                            color[i] = (float) .5 + sigma[i % 3];
                        }
                    }

                }
                break;

            case DEPENDENT:
                for (int i = 0; i < currentColors.length; i += 3) {
                    color[i] = (currentColors[i] / normalizer) + .5f;
//                           color[i] = 2*(currentColors[i]/normalizer)*(currentColors[i]/normalizer )+.5f;
                    if (noiseType == MovieType.GAUSSIAN_MOVIE) {
                        if (color[i] < 1.0 / 257.0) {
                            color[i] = (float) (1.0 / 257.0);
                        } else if (color[i] > 256.0 / 257.0) {
                            color[i] = (float) (256.0 / 257.0);

                        }
                    } else {

                        if (color[i] < .5) {
                            color[i] = (float) .5 - sigma[i % 3];
                        } else {
                            color[i] = (float) .5 + sigma[i % 3];
//                               System.out.println(bwColor);
//                                  if (bwColor == 0)
//                                      color[i] = (float) .5 - sigma[i % 3];
//                                  else
//                                      color[i] = (float) .5 + sigma[i % 3];
//

//                               if (i /3  < frameNumber %(width*height))
//                                   color[i] = (float) .5 - sigma[i % 3];
//                               else
//                                   color[i] = (float) .5 + sigma[i % 3];
                        }
                    }

                    color[i + 1] = color[i];
                    color[i + 2] = color[i];

                }

//                       bwColor = (bwColor + 1) % 2;

//                       frameNumber++;
                break;
        }

    }


    //create a spatial frame by adding together each spatial frequency.
    private float[] createSpatialFrame(float[] colors) {
        Arrays.fill(colors, 0);
        for (int i = 0; i < spatialFreqDepth; i++) {
            colors = addOneSpatialFreq(colors, i);
        }
        return colors;
    }


    //add each spatial frequency to spatial frame.
    private float[] addOneSpatialFreq(float[] colors, int depth) {
        int position;
        int pixelSize;
        int currentWidth, currentHeight, currentNPixels;
        if (depth > 0) {
            //pixelSize = 2^depth;
            pixelSize = 1 << depth;
            //calculate number of pixels in spatial frequency.
            //round up width/pixelSize in case there is a fractional stixel
            //add one so that pixels on the left will not wrap around to pixels on the right.
            currentWidth = (int) Math.ceil( (float) width / (float) pixelSize) + 1;
            currentHeight = (int) Math.ceil( (float) height / (float) pixelSize) + 1;
        } else {
            //depth == 0 will have any of the issues for depth != 0
            pixelSize = 1;
            currentWidth = width;
            currentHeight = height;
        }
        currentNPixels = currentWidth * currentHeight;
        //generate current frequency.
        float[] tempColors = new float[currentNPixels * 3];
        for (int i = 0, index = 0; i < currentNPixels; i++) {
            if (colorType == ColorType.INDEPENDENT) {
                for (int c = 0; c < 3; c++) {
                    position = random.nextShort() & 0xff;
                    tempColors[index++] = lookup[c][position];
                }
            } else {
                position = random.nextShort() & 0xff;
                for (int c = 0; c < 3; c++) {
                    tempColors[index++] = lookup[0][position];
                }

            }
        }

        int horPhase = 0;
        int vertPhase = 0;

        //Calculate the offset of each frequency.
        //Phase is not necessary for depth == 0;
        if (depth != 0) {
//CHECK ARE THESE REALLY RANDOM, will be if random over whole range of 16 bit short.?
            horPhase = (Math.abs(random.nextShort())) % pixelSize;
            vertPhase = (Math.abs(random.nextShort())) % pixelSize;
            //these if statements may only be true if pixelSize > width or height
            //In this case, phase = 0 means than no transitions occur in the visible
            //part of the stimulus.
            if (horPhase > width) {
                horPhase = 0;
            }
            if (vertPhase > height) {
                vertPhase = 0;
            }
        }

        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                for (int c = 0; c < 3; c++) {
                    int offsetJ;
                    int offsetI;
                    //If you are past the phase starting point, get values from
                    //the beginning of the coarse frame.
                    //If you are before the phase starting point, get values from
                    //the end of the coarse frame.
                    if (j >= vertPhase) {
                        offsetJ = (j - vertPhase) / pixelSize;
                    } else {
                        offsetJ = currentHeight - 1 - (vertPhase - j - 1) / pixelSize;
                    }
                    if (i >= horPhase) {
                        offsetI = (i - horPhase) / pixelSize;
                    } else {
                        offsetI = currentWidth - 1 - (horPhase - i - 1) / pixelSize;
                    }

                    //Add the fine noise to the coarse noise.
                    //color(row, column, color) = color(row, column, color) + coarseColor(coarseRow, coarseColumn, color)
//                    if (depth == depth) {
                        colors[ (j * width + i) * 3 + c] =
                            (colors[ (j * width + i) * 3 + c] +
                             tempColors[ (offsetJ * currentWidth + offsetI) * 3 + c]);
//                    }
                }
            }
        }

        return colors;
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

}
