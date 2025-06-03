package edu.ucsc.neurobiology.vision.convert;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Set;

import edu.ucsc.neurobiology.vision.Vision;
import edu.ucsc.neurobiology.vision.calculations.AbstractCalculation;
import edu.ucsc.neurobiology.vision.io.RawDataFile;
import edu.ucsc.neurobiology.vision.io.RawDataHeader512;
import edu.ucsc.neurobiology.vision.io.chunk.GlobalsFile;

public class CopyRawDataHeader512ToGlobals extends AbstractCalculation {
    private String rawdatapath;
    private String globalspath;
    
    @Override
    public void setParameters(HashMap<String, String> parameters) {
        rawdatapath = new File(parameters.get("RawDataPath")).getAbsolutePath();
        globalspath = new File(parameters.get("GlobalsPath")).getAbsolutePath();
    }
    
    @Override
    public void startCalculation() throws IOException {
        GlobalsFile oldGlobals = new GlobalsFile(globalspath, GlobalsFile.READ); // Reads whole file and closes
        if (oldGlobals.rawDataHeader512Exists()) {
            System.err.println("GlobalsFile " + globalspath + " already contains RawDataHeader512");
            finish();
            return;
        }
        
        // Get list of chunks to copy from oldGlobals into updatedGLobals
        Set<Integer> copyTags = oldGlobals.getTags();
        copyTags.remove(GlobalsFile.VERSION_TAG);
        
        // Get RawDataHeader512
        RawDataFile rawdata = new RawDataFile(rawdatapath);
        RawDataHeader512 rdh512;
        try {
            rdh512 = rawdata.getHeader();
        } finally { rawdata.close(); }
        
        // Write RDH512
        GlobalsFile updatedGlobals = new GlobalsFile(globalspath, GlobalsFile.WRITE);
        try {
            oldGlobals.copyTagsTo(updatedGlobals, copyTags);
            updatedGlobals.writeRDH512(rdh512);
        } finally { updatedGlobals.close();	}
        
        finish();
    }
    
    private void finish() {
        Vision.getInstance().sendMessage("CopyRawDataHeader512ToGlobals: Done.");
        Vision.getInstance().getCalculationManager().calculationDone();
    }

}
