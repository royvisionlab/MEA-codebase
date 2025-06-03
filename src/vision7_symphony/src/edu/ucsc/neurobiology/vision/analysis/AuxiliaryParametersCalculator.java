package edu.ucsc.neurobiology.vision.analysis;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.calculations.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.io.chunk.*;
import edu.ucsc.neurobiology.vision.io.chunk.GlobalsFile.*;
import edu.ucsc.neurobiology.vision.stimulus.WhiteNoiseMovie;
import edu.ucsc.neurobiology.vision.electrodemap.*;
import edu.ucsc.neurobiology.vision.util.VisionParams;


public class AuxiliaryParametersCalculator extends AbstractCalculation {
    String workingDirectory;
    
    //default values used when setMovie is false
    int xOffset = 0, yOffset = 80; //pixels, from bottom left
    int refreshInterval = 1; //frames
    int pixelsPerStixelX = 1, pixelsPerStixelY = 1;
    int[] xCorners = null,yCorners = null; //pixels

    //default values used when setElectrodes is false
    private int arrayID = 504, arrayPart = 1, arrayNParts = 1;
    private boolean flipX = false, flipY = false;
    
    private double monitorFrequency;
    private int framesPerTTL;
    
    private boolean setMovie, overrideElectrodes;
    private Vision app;


    public void startCalculation() throws IOException {
        app = Vision.getInstance();

        app.sendMessage("Calculating Auxiliary Parameters...");
        String datasetName = new File(workingDirectory).getName();
        
        String globalsFileName = workingDirectory + File.separator + datasetName + ".globals";
        GlobalsFile oldGlobals = new GlobalsFile(globalsFileName, ChunkFile.READ); // Reads in whole file and closes
        GlobalsFile updatedGlobals = new GlobalsFile(globalsFileName, ChunkFile.WRITE);
        
        // Copy non-updated fields from oldGlobals to updatedGlobals
        Set<Integer> copyTags = oldGlobals.getTags();
        copyTags.remove(GlobalsFile.VERSION_TAG); // This is set automatically when updatedGlobals is created
        copyTags.remove(GlobalsFile.ICP_TAG);
        if (setMovie) copyTags.remove(GlobalsFile.RTMP_TAG);
        oldGlobals.copyTagsTo(updatedGlobals, copyTags);
        
        
        // If overrideElectrodes, this stuff gets set in setParameters
        if (!overrideElectrodes) {
            ImageCalibrationParams imgParams =  oldGlobals.getImageCalibrationParams();
            arrayID = imgParams.arrayID;
            arrayPart = imgParams.arrayPart;
            arrayNParts = imgParams.arrayNParts;
            flipX = imgParams.flipX;
            flipY = imgParams.flipY;
        }
        
        ElectrodeMap map = null;
        int width = 0, height = 0;
        int[] droppedFrameNumbers = null;
        double refreshPeriod = 0.0;
        int nFramesRequired = 0;
        map = ElectrodeMapFactory.getElectrodeMap(arrayID + (arrayPart << 16) + (arrayNParts << 24));
        
        if (setMovie) {
            try {
                map.SetAlignment(xCorners, yCorners, flipX, flipY);
            } catch(Exception e) {
                Vision.reportException(e);
            }

            String neuronFileName = workingDirectory + File.separator + datasetName + ".neurons";
            String movieFileName = workingDirectory + File.separator + datasetName + ".movie";

            if (new File(movieFileName).exists()) {
                WhiteNoiseMovieFile movieFile = new WhiteNoiseMovieFile(movieFileName, ChunkFile.READ);
                WhiteNoiseMovieFile.WhiteMovieParams params = movieFile.getWhiteMovieParams();
                width = params.width;
                height = params.height;
            } else {
                movieFileName = workingDirectory + File.separator + datasetName + ".rawMovie";
                if (new File(movieFileName).exists()) {
                    RawMovieFile rawFile = new RawMovieFile(movieFileName);
                    width = rawFile.getWidth();
                    height = rawFile.getHeight();
                } else {
                    Vision.reportException("Movie file not found.", new IllegalStateException());
                }
            }

            NeuronFile neuronFile = new NeuronFile(neuronFileName);

            int[] ttlTimes = neuronFile.getTTLTimes();
            ArrayList<Integer> droppedFramesList = new ArrayList<Integer> ();
            refreshPeriod = WhiteNoiseMovie.calculateRefreshPeriod(droppedFramesList,
                    ttlTimes, refreshInterval, monitorFrequency, framesPerTTL);

            //The frame number of each dropped frame
            droppedFrameNumbers = new int[droppedFramesList.size()];

            for (int i = 0; i < droppedFramesList.size(); i++) {
                droppedFrameNumbers[i] = (int) (droppedFramesList.get(i));
            }

            //Number of frames required for STA calculation to be possible.
            nFramesRequired = ( ttlTimes.length * 100) / refreshInterval + 1;

            // Copy into new globals
            updatedGlobals.setRunTimeMovieParams(pixelsPerStixelX,  pixelsPerStixelY, width, 
                    height,  pixelsPerStixelX*map.micronsPerPixelX, pixelsPerStixelY*map.micronsPerPixelY,
                    xOffset*map.micronsPerPixelX, yOffset*map.micronsPerPixelY, refreshInterval, monitorFrequency, framesPerTTL,
                    refreshPeriod, nFramesRequired, droppedFrameNumbers);
        }

        updatedGlobals.setImageCalibrationParams(map.micronsPerPixelX, map.micronsPerPixelY, 
            map.centerX, map.centerY, flipX, flipY, map.angle, arrayID, arrayPart, arrayNParts);
        
        
        updatedGlobals.close();
        app.sendMessage("AuxiliaryParametersCalculator: Done.");
        app.getCalculationManager().calculationDone();
    }


    public void setParameters(HashMap<String, String> parameters) {
        workingDirectory = parameters.get("File_Path");
        
        setMovie = Boolean.valueOf(parameters.get("Set Movie"));
        if (setMovie) {
            String cornersString = parameters.get("Set Movie.corners (pixels)");
            String[] cornersStrings = cornersString.split(";");
            xCorners = new int[cornersStrings.length/2];
            yCorners = new int[cornersStrings.length/2];
            try {
                for (int i = 0; i < cornersStrings.length; i += 2) {
                    xCorners[i/2] =Integer.parseInt(cornersStrings[i]);
                    yCorners[i/2] =Integer.parseInt(cornersStrings[i+1]);
                }
            } catch (NumberFormatException e) {
                System.out.println("Electrode map string not provided.  Electrode positions will be modeled from the data, and should not be trusted.");
                xCorners = null;
                yCorners = null;

            }

            pixelsPerStixelX = Integer.parseInt(parameters.get("Set Movie.pixelsPerStixelX")); 
            pixelsPerStixelY = Integer.parseInt(parameters.get("Set Movie.pixelsPerStixelY"));
            xOffset = Integer.parseInt(         parameters.get("Set Movie.xOffset")); //pixels
            yOffset = Integer.parseInt(         parameters.get("Set Movie.yOffset")); //pixels
            refreshInterval = Integer.parseInt( parameters.get("Set Movie.refreshInterval")); //frames
            
            monitorFrequency = VisionParams.DEFAULT_MONITOR_FREQUENCY; // default
            if (parameters.containsKey("Set Movie.Monitor Frequency")) 
                monitorFrequency = Double.valueOf(parameters.get("Set Movie.Monitor Frequency"));

            framesPerTTL = VisionParams.DEFAULT_FRAMES_PER_TTL; // default
            if (parameters.containsKey("Set Movie.Frames Per TTL"))
                framesPerTTL = Integer.valueOf(parameters.get("Set Movie.Frames Per TTL"));
        }
        
        overrideElectrodes = Boolean.valueOf(parameters.get("Override Electrodes From Spike Finding"));
        if (overrideElectrodes) {
            arrayID = Integer.valueOf(    parameters.get("Override Electrodes From Spike Finding.arrayID"));
            arrayPart = Integer.valueOf(  parameters.get("Override Electrodes From Spike Finding.arrayPart"));
            arrayNParts = Integer.valueOf(parameters.get("Override Electrodes From Spike Finding.arrayNParts"));
            flipX = Boolean.valueOf(      parameters.get("Override Electrodes From Spike Finding.flipX")).booleanValue();
            flipY = Boolean.valueOf(      parameters.get("Override Electrodes From Spike Finding.flipY")).booleanValue();
        }
    }
}
