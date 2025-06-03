package edu.ucsc.neurobiology.vision.analysis;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.calculations.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * Calculation that joins the results of farmed out STA calculations
 * Roughly derived from anf.Join
 * 
 * @author Peter H. Li, The Salk Institute
 */
public class STAJoin extends AbstractCalculation {
    String paths;

    public void startCalculation() throws Exception {
        // For now, assume just one directory is given
        File analysisDir = new File(paths);
        final String dataSetName = new File(paths).getName();
        String[] staFilenames = analysisDir.list(new FilenameFilter() {
            public boolean accept(File dir, String name) {
                return name.matches(dataSetName + ".sta-[0-9]+");
            }
        });
        
        // Open STA subfiles
        List<STAFile> staFiles = new ArrayList<STAFile>();
        for (String staFilename : staFilenames)
            staFiles.add(new STAFile(analysisDir + File.separator + staFilename));
        
        // Get total number of neurons
        int nSTAs = 0;
        for (STAFile staFile : staFiles) nSTAs += staFile.getSTACount();
        // TODO: Could check this against the neurons file in the folder...

        String msg = "Found " + StringUtil.listString(staFilenames) + ", " + nSTAs + " total";
        Vision.getInstance().sendMessage(msg);
        System.out.println(msg);
        
        // Create output STA file
        STAFile outSTAFile = new STAFile(
                paths + File.separator + dataSetName + VisionParams.STA_FILE_EXTENSION,
                2 * nSTAs,
                staFiles.get(0).getWidth(),
                staFiles.get(0).getHeight(),
                staFiles.get(0).getSTADepth(),
                staFiles.get(0).getSTAOffset(),
                staFiles.get(0).getStixelWidth(),
                staFiles.get(0).getStixelHeight(),
                staFiles.get(0).getRefreshTime()
        );

        // Write STAs
        int stasSaved = 0;
        Vision.getInstance().startProgressBar();
        for (STAFile staFile : staFiles) {
            int[] idList = staFile.getIDList();
            for (int i = 0; i < idList.length; i++) {
                int id = idList[i];
                outSTAFile.addSTA(id, staFile.getSTA(id));

                Vision.getInstance().setProgress(stasSaved * 100 / nSTAs);
                stasSaved++;
            }        	
        }
        Vision.getInstance().endProgressBar();
        
        // Close files
        for (STAFile staFile : staFiles) staFile.close();
        outSTAFile.close();

        Vision.getInstance().getCalculationManager().calculationDone();
    }


    public void setParameters(HashMap<String, String> parameters) {
        paths = parameters.get("paths").trim();
    }

}
