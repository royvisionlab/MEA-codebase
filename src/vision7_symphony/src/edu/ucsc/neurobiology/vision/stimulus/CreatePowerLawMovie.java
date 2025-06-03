package edu.ucsc.neurobiology.vision.stimulus;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.calculations.*;
import edu.ucsc.neurobiology.vision.io.*;


/**
 * Create a power law noise movie, from the action menu
 *
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class CreatePowerLawMovie
    extends AbstractCalculation {

    String tempMovieFileName, movieFileName;
    int width;
    int height;
    int spatialPeriodCutoff;
    int temporalPeriodCutoff;
    double spatialPower;
    double temporalPower;
    int nFrames;
    int seed;
    int colorType;
    MovieType noiseType;
    float sigma;
    private Vision app;


    public void setParameters(HashMap<String, String> p) {
        String workingDirectory = (String) p.get("File_Path");
        String datasetName = new File(workingDirectory).getName();
        movieFileName = workingDirectory + File.separator + datasetName +
                        ".rawMovie";
        tempMovieFileName = workingDirectory + File.separator + datasetName;


        width = Integer.parseInt( (String) p.get("Width"));
        height = Integer.parseInt( (String) p.get("Height"));
        spatialPeriodCutoff = Integer.parseInt( (String) p.get(
            "White Spatial Period Cutoff"));
        temporalPeriodCutoff = Integer.parseInt( (String) p.get(
            "White Temporal Period Cutoff"));
        spatialPower = Double.parseDouble( (String) p.get("Spatial Power"));
        temporalPower = Double.parseDouble( (String) p.get("Temporal Power"));
        nFrames = Integer.parseInt( (String) p.get("Number of Frames"));
        seed = Integer.parseInt( (String) p.get("Seed"));
        colorType = (int) Double.parseDouble( (String) p.get("ColorType"));
        noiseType = MovieType.values()[ (int) Double.parseDouble( (String) p.get(
            "NoiseType"))];
        sigma = (float) Double.parseDouble( (String) p.get("Sigma"));
    }


    public void startCalculation() throws IOException {
        app = Vision.getInstance();

        app.sendMessage("Calculating Natural Power Movie...");

        try {
            PowerLawFrameGenerator powerLawFrameGenerator = new PowerLawFrameGenerator(
                width,
                height, spatialPeriodCutoff, temporalPeriodCutoff, spatialPower,
                temporalPower, new float[] {sigma, sigma, sigma}
                , seed, colorType, noiseType, nFrames, tempMovieFileName);

            PowerLawMovieFile powerFile = new PowerLawMovieFile(
                movieFileName, width, height, nFrames,
                "PowerLaw-v1", seed, colorType,
                noiseType,
                sigma, spatialPeriodCutoff,
                temporalPeriodCutoff, spatialPower,
                temporalPower);
            
            float[] colorFrame = new float[width * height * 3];

            for (int i = 0; i < nFrames; i++) {
                powerLawFrameGenerator.nextFrame(colorFrame);
                powerFile.writeFrame( (float[]) colorFrame);
            }
            
            powerFile.close();


            powerLawFrameGenerator.deleteTempFile();
        } catch (IOException e) {
            Vision.reportException(e);
        }
        app.sendMessage("Create Natural Power Movie: Done.");
        Vision.getInstance().getCalculationManager().calculationDone();
    }

}
