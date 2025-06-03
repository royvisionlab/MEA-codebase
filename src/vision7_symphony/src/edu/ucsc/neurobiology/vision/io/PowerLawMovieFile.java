package edu.ucsc.neurobiology.vision.io;


import java.io.*;

import edu.ucsc.neurobiology.vision.stimulus.*;


/**
 * This class adds the necessary information to the RawMovieFile header for power law movies.
 *
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class PowerLawMovieFile
    extends RawMovieFile {

    public int seed, colorType, spatialPeriodCutoff, temporalPeriodCutoff;
    public MovieType noiseType;
    public float sigma;
    public double spatialPower, temporalPower;

    public PowerLawMovieFile(String fileName, int width, int height, int framesGenerated,
                             String algorithm, int seed, int colorType,
                             MovieType noiseType,
                             float sigma, int spatialPeriodCutoff,
                             int temporalPeriodCutoff, double spatialPower,
                             double temporalPower) throws IOException {

        super(fileName, width, height, framesGenerated, algorithm);
        this.seed = seed;
        this.colorType = colorType;
        this.noiseType = noiseType;
        this.seed = seed;
        this.sigma = sigma;
        this.spatialPeriodCutoff = spatialPeriodCutoff;
        this.temporalPeriodCutoff = temporalPeriodCutoff;
        this.spatialPower = spatialPower;
        this.temporalPower = temporalPower;

        addToHeader("seed", seed);
        addToHeader("independent", colorType == PowerLawFrameGenerator.RGB);
        addToHeader("distribution",
                    noiseType == MovieType.GAUSSIAN_MOVIE ? "GAUSSIAN" : "BINARY");
        addToHeader("standard-deviation", sigma);
        addToHeader("white-spatial-period-cutoff", spatialPeriodCutoff);
        addToHeader("white-temporal-period-cutoff", temporalPeriodCutoff);
        addToHeader("spatial-power", spatialPower);
        addToHeader("temporal-power", temporalPower);
        writeHeader();
        //    System.out.println("headerLength      " + headerLength);
    }


    public PowerLawMovieFile(String fileName) throws IOException {
        super(fileName);

        seed = (new Integer(getValues("seed")[0])).intValue();
        colorType = (new Boolean(getValues("independent")[0])).booleanValue() ?
                    PowerLawFrameGenerator.RGB : PowerLawFrameGenerator.GRAY;

        noiseType = (getValues("distribution")[0].compareTo("GAUSSIAN")) == 0 ?
                    MovieType.GAUSSIAN_MOVIE : MovieType.BINARY_MOVIE;
        sigma = (new Float(getValues("standard-deviation")[0])).floatValue();
        spatialPeriodCutoff = (new Integer(getValues("white-spatial-period-cutoff")[0])).
                              intValue();
        temporalPeriodCutoff = (new Integer(getValues("white-temporal-period-cutoff")[0])).
                               intValue();
        spatialPower = (new Double(getValues("spatial-power")[0])).doubleValue();
        temporalPower = (new Double(getValues("temporal-power")[0])).doubleValue();
        /*
                System.out.println("seed: " + seed);
         System.out.println("colorType: " + (colorType==PinkFrameGenerator.RGB ? "RGB" : "GRAY"));
                System.out.println("noiseType: " + (noiseType==PinkFrameGenerator.GAUSSIAN_MOVIE ? "GAUSSIAN" : "BINARY"));
                System.out.println("standardDeviation: " + standardDeviation);
                System.out.println("spatialFreqCount: " + spatialFreqCount);
                System.out.println("tempFreqCount: " + tempFreqCount);
         */
    }
}
