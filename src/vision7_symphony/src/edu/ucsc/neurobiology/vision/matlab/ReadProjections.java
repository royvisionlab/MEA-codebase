package edu.ucsc.neurobiology.vision.matlab;

import java.io.IOException;
import edu.ucsc.neurobiology.vision.io.*;


/**
 * Wrapper class for vision ProjectionsFile class
 * This is necessary because Matlab cannot pass/return
 * variables by reference. Consequently, this class is needed
 * to interface with the projections file io code.
 * 
 * Added to Vision 3/09 by tamachado at request of shlens
 * 
 * @author tamachado@salk.edu (1/23/08)
 *
 */

public class ReadProjections {
    private ProjectionsFile projectionsFile;
    private VisionHeader header;
    private int electrode; //most recently read electrode ("current electrode")
    private float[][] data;//projections from the current electrode
    private int[]   times; //spike times from the current electrode
    private int nSpikes;   //number of spikes on the current electrode

    /**
     * Create or open a neurons file to write to.
     * 
     * @param path     path to projections file
     */
    public ReadProjections (String path) {
        
        try {
           projectionsFile = new ProjectionsFile(path);
        } catch (IOException e) {
           e.printStackTrace();
        }
        
        header = projectionsFile.getHeader();
        int dims = header.nDimensions;
        int max  = projectionsFile.maxSpikesPerElectrode;
        
        data  = new float[dims][max];
        times = new int[max];
        nSpikes = 0;
        electrode = 0;
    }
    
    /**
     * The member variables are updated after this method is called
     * so after calling this method, use the getter methods to examine
     * the data stored on electrode.
     * 
     * 
     * @param electrode   read from electrode
     * @return boolean    true upon success
     */
    public boolean readElectrode(int electrode) {
        try {
            
           //read the projections in from the prj file 
           nSpikes = projectionsFile.readProjections(electrode, data, times);
           
           //if there are no spikes, this is a bad electrode
           if (nSpikes <= 0) {
               return false;
           }
           
           //update current electrode being viewed
           this.electrode = electrode;
           return true;
           
        } catch (IOException e) {
           e.printStackTrace();
           return false;
        }
    }
    
    
    //getter functions
    public int[] getTTLPulses() {
        return projectionsFile.getTTLTimes();
    }
    
    public ProjectionsFile getProjectionsFile() {
        return projectionsFile;
    }
    
    public int[] getSpikeTimes() {
        return times;
    }
    
    public float[][] getProjections() {
        return data;
    }
    
    public int getSpikeCount() {
        return nSpikes;
    }

    public int getElectrode() {
        return electrode;
    }
}
