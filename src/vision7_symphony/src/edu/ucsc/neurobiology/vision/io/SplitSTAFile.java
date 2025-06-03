package edu.ucsc.neurobiology.vision.io;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Set;

import edu.ucsc.neurobiology.vision.stimulus.STA;

public class SplitSTAFile implements STACollection {
    private STAFile[] staFiles;
    private Hashtable<Integer,Integer> neuronIDMap;

    /**
     * Use a basefile name, e.g. data000.sta to find the split files, e.g. data000.sta-0, data000.sta-1, etc
     * @param baseFilename
     * @throws IOException 
     */
    public static SplitSTAFile fromBaseFilename(String baseFilename) throws IOException {
        File baseFile = new File(baseFilename);

        // If this is a dir, use it, otherwise get parent
        File baseDir;
        if (baseFile.isDirectory()) baseDir = baseFile;
        else baseDir = baseFile.getParentFile();

        // If we still haven't found a dir, give up
        if (!baseDir.isDirectory()) return null;

        // Get all files in parent dir matching ".sta-"
        File[] staFiles = baseDir.listFiles(new FilenameFilter() {
            public boolean accept(File dir, String filename) {
                if (filename.contains(".sta-")) return true;
                return false;
            }
        });
        
        if (staFiles.length == 0) return null;
        return new SplitSTAFile(staFiles);
    }
    
    
    public SplitSTAFile(File[] files) throws IOException {
        staFiles = new STAFile[files.length];
        for (int i = 0; i < staFiles.length; i++) staFiles[i] = new STAFile(files[i]);
        
        // Check that STAFiles are consistent
        for (int i = 1; i < staFiles.length; i++) {
            if (!staFiles[i].isConsistentWith(staFiles[0])) {
                String errString = "STAFiles not consistent: " + files[0].getName() + ", " + files[i].getName();
                throw new IOException(errString);
            }
        }
        
        // Build neuronID map
        neuronIDMap = new Hashtable<Integer,Integer>();
        for (int filenum = 0; filenum < staFiles.length; filenum++) {
            int[] neuronIDs = staFiles[filenum].getIDList();
            for (int id : neuronIDs) neuronIDMap.put(id, filenum);
        }
    }
    
    public STA getSTA(int neuronID) throws IOException {
        Integer filenum = neuronIDMap.get(neuronID);
        if (filenum == null) return null;
        return staFiles[filenum].getSTA(neuronID);
    }

    public int getWidth()           { return staFiles[0].getWidth(); }
    public int getHeight()          { return staFiles[0].getHeight(); }
    public double getStixelWidth()  { return staFiles[0].getStixelWidth(); }
    public double getStixelHeight() { return staFiles[0].getStixelHeight();	}
    public int getSTADepth()        { return staFiles[0].getSTADepth(); }
    public int getSTAOffset()       { return staFiles[0].getSTAOffset(); }
    public double getRefreshTime()  { return staFiles[0].getRefreshTime(); }	
    
    public void close() throws IOException {
        for (STAFile staFile : staFiles) staFile.close();
    }


    public int[] getIDList() {
        Set<Integer> ids = neuronIDMap.keySet();
        Iterator<Integer> idIterator = ids.iterator();
        int[] r = new int[ids.size()];
        for (int i = 0; i < r.length; i++) r[i] = idIterator.next();
        return r;
    }
}
