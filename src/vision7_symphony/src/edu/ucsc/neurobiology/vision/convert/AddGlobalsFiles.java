package edu.ucsc.neurobiology.vision.convert;

import java.io.*;
import java.util.HashMap;

import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.UnflaggedOption;


import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.io.chunk.*;
import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.util.StringUtil;
import edu.ucsc.neurobiology.vision.util.VisionJSAP;
import edu.ucsc.neurobiology.vision.util.VisionParams;
import edu.ucsc.neurobiology.vision.calculations.*;

public class AddGlobalsFiles extends AbstractCalculation {
    
    private String inputPath;
    private boolean fullTree;
    
    public void startCalculation() throws IOException {
        topLevel(inputPath, fullTree);
        Vision.getInstance().getCalculationManager().calculationDone();
    }
    
    
    public void topLevel(String rootDirectoryName, boolean fullTree) throws IllegalStateException {
        File root = new File(rootDirectoryName);
        if (!root.isDirectory()) throw new IllegalStateException("Requested root directory does not exist.");
        
        boolean error = true;
        while (error) {
            error = false;
            try {
                process(root, fullTree);
            } catch (IOException e) {
                e.printStackTrace();
                error = true;
                try {
                    Thread.sleep(1000);
                } catch (InterruptedException ex) {
                    ex.printStackTrace();
                }
                
            }
        }
    }
    
    
    private void process(File dir, boolean fullTree) throws IOException {
        System.out.println("Generating globals file: " + dir.getAbsolutePath());
        
        File[] files = dir.listFiles(); //includes directories
        
        // find the name of the neurons file in dir
        String name = findName(files, VisionParams.NEURON_FILE_EXTENSION);
        if (name != null) {
            
            // check if globals file exists, if it does then don't make a new one
            String globalsName = name + VisionParams.GLOBALS_FILE_EXTENSION;
            String globalsPath = dir.getAbsolutePath() + "//" + globalsName;
            if (new File(globalsPath).exists()) {
                System.err.println("A globals file already exists in this directory. Skipping...");
                return;
            }
            
            // otherwise make a new one
            System.err.println("Creating a new globals file...");
            GlobalsFile gf = null;
            String staName     = name + VisionParams.STA_FILE_EXTENSION;
            String neuronsName = name + VisionParams.NEURON_FILE_EXTENSION;
            String eiName      = name + VisionParams.EI_FILE_EXTENSION;
            
            if (contains(files, neuronsName)) {
                int rawArrayID = 0;
                int arrayID = 0;
                int part = 0;
                int nParts = 1;
                if (contains(files, eiName)) {
                    PhysiologicalImagingFile eif = new PhysiologicalImagingFile(dir.getAbsolutePath() + "//" + eiName);
                    rawArrayID = eif.getArrayID();
                    arrayID = rawArrayID & 0xFFFF;
                    part = (rawArrayID >> 16) & 0xFF;
                    nParts = rawArrayID >> 24;
                }

                gf = new GlobalsFile(globalsPath, GlobalsFile.WRITE);
                gf.setImageCalibrationParams(1, 1, 0, 0, false, false, 0, arrayID, part, nParts);
            } else System.err.println("Skipping directory, no neurons file found!");

            if (gf != null && contains(files, staName)) {
                RandomAccessFile sf = new RandomAccessFile(dir.getAbsolutePath() + "//" + staName, "r");
                sf.readInt();
                sf.readInt();
                int width = sf.readInt();
                int height = sf.readInt();
                sf.readInt();
                int pixelsPerStixel = (int) sf.readDouble(); //reads deprecated pixel size parameter.
                double refreshTime = sf.readDouble();
                int interval = (int) Math.round(refreshTime/8.33);
                refreshTime /= interval;

                sf.close();

                double monitorFrequency = VisionParams.DEFAULT_MONITOR_FREQUENCY;
                int framesPerTTL = VisionParams.DEFAULT_FRAMES_PER_TTL;
                System.out.println("Using default monitorFrequency " + monitorFrequency + " and default framesPerTTL " + framesPerTTL);
                gf.setRunTimeMovieParams(pixelsPerStixel, pixelsPerStixel, width, height, 
                        5.8*pixelsPerStixel, 5.8*pixelsPerStixel, 0.0, 0.0, interval, 
                        monitorFrequency, framesPerTTL,
                        refreshTime, -1, new int[0]);
            }

            if (gf != null) {
                gf.writeCreatedBy(GlobalsFile.UPDATER);
                gf.close();
            }
        }
        
        for (int i = 0; i < files.length; i++) {
            if (files[i].isDirectory() && fullTree) {
                boolean error = true;
                while (error) {
                    error = false;
                    try {
                        process(files[i], fullTree);
                    } catch (IOException e) {
                        e.printStackTrace();
                        error = true;
                        try {
                            Thread.sleep(1000);
                        } catch (InterruptedException ex) {
                            ex.printStackTrace();
                        }

                    }
                }
                //process(files[i], fullTree);
            }
        }
    }

    private String findName(File[] files, String extension) {
        //store information about current file
        String name, subName;
        int dotIndex;

        //store name of the neurons file
        String fileName = null;

        //return null if there is more than one file with extension
        for(int i=0; i<files.length; i++) {
            name = files[i].getName();
            dotIndex = name.lastIndexOf(".");
            
            //if there is no dot (or if the file is hidden), then this isn't a valid .neurons file
            if(dotIndex == -1 || name.charAt(0) == '.') continue;
            
            //get the file extension
            subName = name.substring(dotIndex);
            
            if(extension.matches("\\Q" + subName  + "\\E") && files[i].isFile()) {
                if(fileName != null) {
                    System.err.println("There are MULTIPLE neurons files in this directory. Skipping...");
                    return null;
                }
                fileName = name.substring(0,dotIndex);
            }
        }

        if (fileName == null) {
            System.err.println("There are NO neurons files in this directory. Skipping...");
        }
        
        return fileName;

    }

    private boolean contains(File[] files, String name) {

        for(int i=0; i<files.length; i++) {
            if(name.matches("\\Q" + files[i].getName() + "\\E") && files[i].isFile()) return true;
        }
        return false;
    }

    public void setParameters(HashMap<String, String> parameters) {
        inputPath = ( (String) parameters.get("rootPath")).trim();
        fullTree = Boolean.valueOf( (String) parameters.get("fullTree")).booleanValue();
    }
    
    public static void main(String[] args) throws Exception {
        long t1 = System.currentTimeMillis();

        VisionJSAP jsap = new VisionJSAP( 
                AddGlobalsFiles.class.getName(), 
                new com.martiansoftware.jsap.Parameter[] {
                        new UnflaggedOption("rootPath", JSAP.STRING_PARSER, JSAP.REQUIRED, "Root analysis folder."),
                        new UnflaggedOption("fullTree", JSAP.BOOLEAN_PARSER, JSAP.REQUIRED, "Recurse subdirectories?"),
                }
        );

        JSAPResult parsedArgs = jsap.parse(args);

//        HashMap<String,String> parameters;
        boolean fullTree = parsedArgs.getBoolean("fullTree");
        String rootPath = parsedArgs.getString("rootPath");
        File dir = new File(rootPath);
        
        if (!dir.exists()) {
            System.out.println("Root directory is invalid!");
            System.exit(1);
        }
        
        AddGlobalsFiles adder = new AddGlobalsFiles();
        adder.topLevel(rootPath,fullTree);

        long t2 = System.currentTimeMillis();
        System.out.print("\nThe whole analysis took: ");
        double t = (t2 - t1) / 1000.0;
        if (t < 60) {
            System.out.println(t + " sec.");
        } else if (t < 3600) {
            System.out.println(StringUtil.format(t / 60.0, 1) + " min.");
        } else {
            System.out.println(StringUtil.format(t / 3600.0, 1) + " hours.");
        }
        
        System.exit(0);
    }

}
