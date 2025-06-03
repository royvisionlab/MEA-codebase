package edu.ucsc.neurobiology.vision.stimulus;

import java.io.File;
import java.util.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.calculations.*;
import edu.ucsc.neurobiology.vision.io.chunk.*;
import edu.ucsc.neurobiology.vision.stimulus.FrameGenerator.*;


/**
 * A Vision calculation that creates a white noise movie file.  Only parameters are created.
 * Frames are generated as needed.
 *
 *
 * @author Matthew Grivich, The Salk Institute
 */
public class CreateWhiteNoiseMovie extends AbstractCalculation {

    String workingDirectory;
    int width;
    int height;
    int seed;
    ColorType colorType;
    MovieType noiseType;
    RandomNumberGenerator rng;
    double contrast;
    boolean sparse = false;
    float probability;


    public void startCalculation() throws Exception {
        
        String datasetName = new File(workingDirectory).getName();	
        String movieFileName = workingDirectory + File.separator + datasetName +
        ".movie";

        
        WhiteNoiseMovieFile file = new WhiteNoiseMovieFile(movieFileName, ChunkFile.WRITE);
        file.setWhiteMovieParams(width, height, seed, colorType, noiseType, rng, contrast, sparse, probability);
        file.close();
        Vision.getInstance().getCalculationManager().calculationDone();
    }


    public void setParameters(HashMap<String, String> p) {
        workingDirectory = p.get("filePath");
        noiseType = MovieType.values()[ (int) Double.parseDouble(p.get("NoiseType"))];
        colorType = ColorType.values()[ (int) Double.parseDouble(p.get("ColorType"))];
        rng = RandomNumberGenerator.values()[ (int) Double.parseDouble(p.get("RandomNumberGenerator"))];
        seed = Integer.parseInt(p.get("Seed"));
        width = Integer.parseInt(p.get("Width"));
        height = Integer.parseInt(p.get("Height"));
        contrast = Double.parseDouble(p.get("ContrastSigma"));
        if (p.containsKey("Sparse")) {
            sparse = Boolean.parseBoolean(p.get("Sparse"));
            if (sparse) probability = Float.parseFloat(p.get("Sparse.Probability"));
        }
    }

}
