package edu.ucsc.neurobiology.vision.actions;

import java.io.*;

import java.awt.event.*;
import javax.swing.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.parameters.*;
import edu.ucsc.neurobiology.vision.stimulus.*;


/**
 *
 * @author Matthew Grivich, University of California, Santa Cruz
 */

/*Create a pink noise movie, from the action menu*/


public class CreatePinkNoiseMovie
    extends AbstractAction {


    public CreatePinkNoiseMovie() {
        super("Create Pink-Noise Movie");
    }


    public void actionPerformed(ActionEvent event) {
        Vision app = Vision.getInstance();

        ParametersTable table = app.getConfig().showDialog(
            "CreatePinkNoiseMovie", "Create Pink Noise Movie", app.getMainFrame());
        if (table == null) {
            return;
        }

        String workingDirectory = ( (FileParameter) table.getParameter("File_Path")).
                                  getValue();
        String datasetName = new File(workingDirectory).getName();
        String movieFileName = workingDirectory + File.separator + datasetName +
                               ".rawMovie";

        int width = ( (IntegerParameter) table.getParameter("Width")).getValue();
        int height = ( (IntegerParameter) table.getParameter("Height")).getValue();
        int spatialFreqCount = ( (IntegerParameter) table.getParameter(
            "Spatial Frequency Count")).getValue();
        int tempFreqCount = ( (IntegerParameter) table.getParameter(
            "Temporal Frequency Count")).getValue();
        int nFrames = ( (IntegerParameter) table.getParameter("Number of Frames")).
                      getValue();
        int seed = ( (IntegerParameter) table.getParameter("Seed")).getValue();
        int colorType = (int) ( (EnumeratorParameter) table.getParameter("ColorType")).
                        getValue();
        MovieType noiseType =
            MovieType.values()[ (int) ( (EnumeratorParameter) table.getParameter(
                "NoiseType")).
            getValue()];

        float sigma = (float) ( (DoubleParameter) table.getParameter("ContrastSigma")).
                      getValue();

        PinkFrameGeneratorV2 pinkFrameGenerator = new PinkFrameGeneratorV2(
            width, height,
            spatialFreqCount,
            tempFreqCount, new float[] {sigma, sigma, sigma}
            , seed, colorType, noiseType);
        PinkMovieFile pinkFile = null;
        try {
            pinkFile = new PinkMovieFile(movieFileName, width, height, nFrames,
                                         "Voss-McCartney-v1", seed, colorType,
                                         noiseType,
                                         sigma, spatialFreqCount,
                                         tempFreqCount);
            
            float[] colorFrame = new float[width * height * 3];

            for (int i = 0; i < nFrames; i++) {
                pinkFrameGenerator.nextFrame(colorFrame);
                pinkFile.writeFrame( (float[]) colorFrame);
            }
            
            pinkFile.close();
            
        } catch (IOException e) {
            Vision.reportException(e);
        }
    }

}
