package edu.ucsc.neurobiology.vision.io;

import java.io.*;

import edu.ucsc.neurobiology.vision.stimulus.*;


/**
 * This class adds the necessary information to the RawMovieFile header for pink movies.
 *
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class PinkMovieFile
    extends RawMovieFile {

    public int seed, colorType, spatialFreqCount, tempFreqCount;
    private MovieType noiseType;
    public float standardDeviation;

    public PinkMovieFile(String fileName, int width, int height, int framesGenerated,
                         String algorithm, int seed, int colorType, MovieType noiseType,
                         float standardDeviation, int spatialFreqCount,
                         int tempFreqCount) throws IOException {

        super(fileName, width, height, framesGenerated, algorithm);
        this.seed = seed;
        this.colorType = colorType;
        this.noiseType = noiseType;
        this.spatialFreqCount = spatialFreqCount;
        this.tempFreqCount = tempFreqCount;
        this.standardDeviation = standardDeviation;
        addToHeader("seed", seed);
        addToHeader("independent", colorType == PinkFrameGeneratorV2.RGB);
        addToHeader("distribution",
                    noiseType == MovieType.GAUSSIAN_MOVIE ? "GAUSSIAN" :
                    "BINARY");
        addToHeader("standard-deviation", standardDeviation);
        addToHeader("N-spatial-frequencies", spatialFreqCount);
        addToHeader("N-temporal-frequencies", tempFreqCount);
        writeHeader();
        //    System.out.println("headerLength      " + headerLength);
    }


    public PinkMovieFile(String fileName) throws IOException {
        super(fileName);

        seed = (new Integer(getValues("seed")[0])).intValue();
        colorType = (new Boolean(getValues("independent")[0])).booleanValue() ?
                    PinkFrameGeneratorV2.RGB : PinkFrameGeneratorV2.GRAY;

        noiseType = (getValues("distribution")[0].compareTo("GAUSSIAN")) == 0 ?
                    MovieType.GAUSSIAN_MOVIE : MovieType.BINARY_MOVIE;
        standardDeviation = (new Float(getValues("standard-deviation")[0])).floatValue();
        spatialFreqCount = (new Integer(getValues("N-spatial-frequencies")[0])).
                           intValue();
        tempFreqCount = (new Integer(getValues("N-temporal-frequencies")[0])).intValue();
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
