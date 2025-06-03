package edu.ucsc.neurobiology.vision.stimulus;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.calculations.*;
import edu.ucsc.neurobiology.vision.io.*;


/**
 * Create a power law noise movie, from the action menu
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class CreateOMSMovie
    extends AbstractCalculation {

    private int width;
    private int height;
    private int nFrames;
    private int seed;
    private float sigma;
    private double objectRadius;
    private double objectSpacing;
    private double stepSize;
    private int barWidth;
    private int runLength;
    private boolean showBackground;
    private boolean randomizeMotion;
    private String movieFileName;


    public void startCalculation() throws IOException {
        OMSFrameGenerator frameGenerator = new OMSFrameGenerator(width, height,
            new float[] { (float) sigma, (float) sigma, (float) sigma}
            , seed, nFrames, objectRadius, objectSpacing, stepSize, barWidth,
            runLength, showBackground, randomizeMotion);
        OMSMovieFile file = null;
        file = new OMSMovieFile(
            movieFileName, width, height, nFrames, "OMS-v1", seed, sigma, objectRadius,
            objectSpacing, stepSize, barWidth, runLength, showBackground, randomizeMotion);
        
        
        float[] colorFrame = new float[width * height * 3];

        for (int i = 0; i < nFrames; i++) {
            frameGenerator.nextFrame(colorFrame);
            file.writeFrame( (float[]) colorFrame);
        }
        
        file.close();

        Vision.getInstance().getCalculationManager().calculationDone();
    }


    public void setParameters(HashMap<String, String> p) {
        String workingDirectory = (String) p.get("File_Path");
        String datasetName = new File(workingDirectory).getName();
        movieFileName = workingDirectory + File.separator + datasetName +
                        ".rawMovie";
        width = Integer.parseInt( (String) p.get("Width"));
        height = Integer.parseInt( (String) p.get("Height"));
        nFrames = Integer.parseInt( (String) p.get("Number of Frames"));
        seed = Integer.parseInt( (String) p.get("Seed"));
        sigma = Float.parseFloat( (String) p.get("Sigma"));
        objectRadius = Double.parseDouble( (String) p.get("Object Radius"));
        objectSpacing = Double.parseDouble( (String) p.get("Object Spacing"));
        stepSize = Double.parseDouble( (String) p.get("Step Size"));
        barWidth = Integer.parseInt( (String) p.get("Bar Width"));
        runLength = Integer.parseInt( (String) p.get("Run Length"));
        showBackground = Boolean.valueOf( (String) p.get("Show Background")).booleanValue();
        randomizeMotion = Boolean.valueOf( (String) p.get("Randomize Motion")).
                          booleanValue();
    }

}
