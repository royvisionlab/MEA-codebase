package edu.ucsc.neurobiology.vision.stimulus;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.math.*;


/**
 * This class defines how to create each frame of a object-motion sensor movie.
 *
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class OMSFrameGenerator
    extends FrameGenerator {


    private float[] sigma; //standard deviation of final stimulus.

    private float normalizer = 1.0f; //used to normalize the frames to the correct standard deviation
    private float[] currentColors; //current frame

    int nFrames;

    FileInputStream[] fis;
    BufferedInputStream[] fileReader;


    RandomMersenne mersenne;

    int objectCount;
    double objectRadius;
    double objectSpacing;
    double stepSize;
    double[][] centers;
    ArrayList centerList;
    double[][] gratingCenters;
    double[] backgroundCenter = {60.0, 0.0};
    int barWidth;
    int runLength;
    int frameNumber = 0;
    boolean showBackground, randomizeMotion;


    public OMSFrameGenerator(int width, int height, float[] sigma, int startSeed,
                             int nFrames, double objectRadius,
                             double objectSpacing, double stepSize, int barWidth,
                             int runLength,
                             boolean showBackground, boolean randomizeMotion) {

        super(width, height, RandomNumberGenerator.JAVA_RANDOM, 0);
        this.sigma = sigma;
        this.nFrames = nFrames;
        this.objectRadius = objectRadius;
        this.objectSpacing = objectSpacing;
        this.stepSize = stepSize;
        this.barWidth = barWidth;
        this.runLength = runLength;
        this.showBackground = showBackground;
        this.randomizeMotion = randomizeMotion;

        currentColors = new float[width * height * 3];
        mersenne = new RandomMersenne(startSeed);

        initializeObjects();

    }


    private void initializeObjects() {
        // FIXME

        // Remainder of 32-bit signed random integer
        // This code generates a random signed integer and then computes the remainder
        // of that value modulo another value. Since the random number can be negative,
        // the result of the remainder operation can also be negative. Be sure this is
        // intended, and strongly consider using the Random.nextInt(int) method instead.

        // leave it like this for compatibility
        int xPosInitial = mersenne.nextInt() % ( (int) objectSpacing / 2);
        int yPosInitial = mersenne.nextInt() % ( (int) objectSpacing / 4);
        int xPos = xPosInitial;
        int yPos = yPosInitial;
        int row = 0;
        ArrayList centerList = new ArrayList();
        while (yPos < height + objectRadius) {
            while (xPos < width + objectRadius) {

                centerList.add(new double[] {xPos, yPos});
                xPos += objectSpacing;
            }
            row++;
            xPos = xPosInitial - ( ( (int) objectSpacing / 2) * (row % 2));
            yPos += objectSpacing / 2;

        }
        objectCount = centerList.size();
        centers = new double[objectCount][2];
        gratingCenters = new double[objectCount][2];
        for (int i = 0; i < centers.length; i++) {

            centers[i][0] = ( (double[]) centerList.get(i))[0]; //(double) width /2;
            centers[i][1] = ( (double[]) centerList.get(i))[1]; //(double) height /2;
            gratingCenters[i][0] = 0;
            gratingCenters[i][1] = 0;

        }

    }


    //Use this function to get new run centers without writing frames.
    //Used for analysis puproses.
    private double[][] getOneRunCenters() {
        for (int i = 0; i < runLength; i++) {
            if (showBackground) {
                mersenne.nextInt();
            }
            if (randomizeMotion) {
                for (int j = 0; j < gratingCenters.length; j++) {
                    mersenne.nextInt();
                }
            }
        }
        initializeObjects();
        return centers;
    }


    public double[][][] getRunCenters(int nRuns) {
        double[][][] allCenters = new double[nRuns][][];
        allCenters[0] = (double[][]) centers.clone();
        for (int i = 1; i < nRuns; i++) {
            allCenters[i] = (double[][]) getOneRunCenters().clone();

        }
        return allCenters;
    }


    public Object createFrameBuffer() {
        return new int[width * height * 3];
    }


    public void nextFrame(Object frame) {
        if (frame instanceof float[]) {
            nextFrame( (float[]) frame);
        } else if (frame instanceof int[]) {
            throw new IllegalStateException(
                "Frames are not created as int[] with power law movies.");
        }

        frameNumber = (frameNumber + 1) % runLength;
        if (frameNumber == 0) {
            initializeObjects();
        }
    }


    private final void nextFrame(float[] color) {

        currentColors = createFrame(currentColors);

        for (int i = 0; i < currentColors.length; i += 3) {
//                    System.out.println(currentColors[i]);
            color[i] = currentColors[i];
//                    color[i] = ((currentColors[i]-.5f)*sigma[0]/.5f +.5f) ;
            if (color[i] == 0) {
                color[i] = (float) .5 - sigma[i % 3];
            } else if (color[i] == 1) {
                color[i] = (float) .5 + sigma[i % 3];
            } else {
                color[i] = .5f;
            }

            color[i + 1] = color[i];
            color[i + 2] = color[i];

        }

    }


    private float[] createFrame(float[] colors) {
        Arrays.fill(colors, 0);

        int direction = ( (mersenne.nextInt() & 0x01) == 0) && randomizeMotion ? -1 : 1;

        backgroundCenter[0] = backgroundCenter[0] + stepSize * direction;

        for (int i = 0; i < gratingCenters.length; i++) {
            direction = ( (mersenne.nextInt() & 0x01) == 0) && randomizeMotion ? -1 : 1;
            gratingCenters[i][0] = gratingCenters[i][0] + stepSize * direction;

        }
        int intCent = (int) Math.floor(backgroundCenter[0]);
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {

                for (int c = 0; c < 3; c++) {

                    //paint background
                    if (x - intCent >= 0 && showBackground) {

                        colors[3 * (y * width + x) +
                            c] = (float) ( (x - intCent) /
                                          barWidth) % 2;
                    } else if (showBackground) {
                        colors[3 * (y * width + x) +
                            c] = (float) ( (barWidth - 1 + intCent - x) /
                                          barWidth) % 2;

                    } else {
                        colors[3 * (y * width + x) + c] = .5f;
                    }
                    //paint objects
                    for (int i = 0; i < centers.length; i++) {
                        int objCentX = (int) Math.floor(centers[i][0]);
                        int objCentY = (int) Math.floor(centers[i][1]);
                        int gratCentX = (int) Math.floor(gratingCenters[i][0]);
                        if ( ( (objCentX - x) *
                              (objCentX - x) +
                              (objCentY - y) *
                              (objCentY - y)) <
                            objectRadius * objectRadius) {
                            if (x - gratCentX >= 0) {

                                colors[3 * (y * width + x) +
                                    c] = (float) ( (x - gratCentX) /
                                                  barWidth) % 2;
                            } else {
                                colors[3 * (y * width + x) +
                                    c] = (float) ( (barWidth - 1 + gratCentX - x) /
                                                  barWidth) % 2;
                            }
                        }
                    }
                }
            }
        }

        return colors;
    }

}
