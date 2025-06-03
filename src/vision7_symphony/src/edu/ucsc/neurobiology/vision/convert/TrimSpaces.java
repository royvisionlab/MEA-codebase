package edu.ucsc.neurobiology.vision.convert;

import java.io.*;
import java.util.HashMap;


import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.calculations.*;

/*
 * Utility to strip outlying spaces on file names for a full directory.
 * 
 * @author Matthew Grivich
 */

public class TrimSpaces extends AbstractCalculation {
    
    private String inputPath;
    private boolean fullTree;
    
    public void startCalculation() throws IOException {
        topLevel(inputPath, fullTree);
        Vision.getInstance().getCalculationManager().calculationDone();
    }

    
    public void topLevel(String rootDirectoryName, boolean fullTree) throws IllegalStateException {
        File root = new File(rootDirectoryName);
        if(!root.isDirectory()) {
            throw new IllegalStateException("Requested root directory does not exist.");
            
        } else {
            try {
            process(root, fullTree);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        
    }
    
    private void process(File dir, boolean fullTree) throws IOException {
        System.out.println("Stripping spaces: " + dir.getAbsolutePath());
        File[] files = dir.listFiles(); //includes directories
        
        for(int i=0; i<files.length; i++) {
            
            String name = (files[i].getAbsolutePath()).trim();
            if(!name.matches("\\Q" + files[i].getName() + "\\E")) {
                files[i].renameTo(new File(name));
            }
        }
        files = dir.listFiles();
        
        
        for(int i=0; i<files.length; i++) {
            if(files[i].isDirectory() && fullTree) {
                process(files[i], fullTree);
            }
        }
    }
    

    
    public void setParameters(HashMap<String, String> parameters) {
        inputPath = ( (String) parameters.get("rootPath")).trim();

        fullTree = Boolean.valueOf( (String) parameters.get("fullTree")).
                   booleanValue();

    }
    
    public static void main(String args[]){
        if(args.length != 2) {
            System.out.println("Syntax: java -cp Vision.jar edu.ucsc.neurobiololgy.vision.convert.TrimSpaces inputPath fullTree(true or false)");
            System.exit(1);
        }
        TrimSpaces strip = new TrimSpaces();
        strip.inputPath = args[0];
        strip.fullTree = Boolean.valueOf( args[1]).booleanValue();
        try{
        strip.topLevel(strip.inputPath, strip.fullTree);
        } catch(IllegalStateException e) {
            e.printStackTrace();
        }
        
    }



}
