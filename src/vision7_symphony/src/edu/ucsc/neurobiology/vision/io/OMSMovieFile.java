package edu.ucsc.neurobiology.vision.io;


import java.io.*;


/**
 * This class adds the necessary information to the RawMovieFile header for power law movies.
 *
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class OMSMovieFile
    extends RawMovieFile {

//    public int seed, barWidth, runLength;
//    public float sigma;
//    public double objectRadius, objectSpacing, stepSize;



    public OMSMovieFile(String fileName, int width, int height, int framesGenerated,
                        String algorithm, int seed,
                        float sigma, double objectRadius,
                        double objectSpacing, double stepSize, int barWidth,
                        int runLength, boolean showBackground, boolean randomizeMotion) throws
        IOException {

        super(fileName, width, height, framesGenerated, algorithm);
//        this.seed = seed;
//
//        this.seed = seed;
//        this.sigma = sigma;
//        this.objectRadius = objectRadius;
//        this.objectSpacing = objectSpacing;
//        this.stepSize = stepSize;
//        this.barWidth = barWidth;
//        this.runLength =runLength;


        addToHeader("seed", seed);
        addToHeader("standard-deviation", sigma);
        addToHeader("object-radius", objectRadius);
        addToHeader("object-spacing", objectSpacing);
        addToHeader("step-size", stepSize);
        addToHeader("bar-width", barWidth);
        addToHeader("run-length", runLength);
        addToHeader("show-background", showBackground);
        addToHeader("randomize-motion", randomizeMotion);

        writeHeader();
        //    System.out.println("headerLength      " + headerLength);
    }


    public int getSeed() {
        return (new Integer(getValues("seed")[0])).intValue();
    }


    public float getSigma() {
        return (new Float(getValues("standard-deviation")[0])).floatValue();
    }


    public double getObjectRadius() {
        return (new Double(getValues("object-radius")[0])).doubleValue();
    }


    public double getObjectSpacing() {
        return (new Double(getValues("object-spacing")[0])).doubleValue();
    }


    public double getStepSize() {
        return (new Double(getValues("step-size")[0])).doubleValue();
    }


    public int getBarWidth() {
        return (new Integer(getValues("bar-width")[0])).intValue();
    }


    public int getRunLength() {
        return (new Integer(getValues("run-length")[0])).intValue();
    }


    public boolean getShowBackground() {
        return (new Boolean(getValues("show-background")[0])).booleanValue();
    }


    public boolean getRandomizeMotion() {
        return (new Boolean(getValues("randomize-motion")[0])).booleanValue();
    }


    public OMSMovieFile(String fileName) throws IOException {
        super(fileName);
    }
}
