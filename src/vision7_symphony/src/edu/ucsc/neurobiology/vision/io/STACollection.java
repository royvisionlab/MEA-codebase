package edu.ucsc.neurobiology.vision.io;

import java.io.File;
import java.io.IOException;

import edu.ucsc.neurobiology.vision.stimulus.STA;

/**
 * Abstracted out some elements from STAFile to allow NeuronViewer to work with STAFileSplit
 * @author peterli
 *
 */
public interface STACollection {
    public int[] getIDList();
    public STA getSTA(int neuronID) throws IOException;
    public double getRefreshTime();
    public int getHeight();
    public int getWidth();
    public double getStixelWidth();
    public double getStixelHeight();
    public int getSTADepth();
    public int getSTAOffset();
    
    /**
     * This should be designed so that it is safe to run it more than once.
     */
    public void close() throws IOException;
    
    public class Factory {
        public static STACollection fromFilename(String path) throws IOException {
            File f = new File(path);
            
            STACollection r;
            if (f.exists() && f.canRead()) r = new STAFile(path);	
            else r = SplitSTAFile.fromBaseFilename(path);
            
            return r;
        }
    }
    
}
