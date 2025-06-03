package edu.ucsc.neurobiology.vision.matlab;

import java.io.File;
import java.io.IOException;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.util.VisionParams;

/**
 * Wrapper class for vision NeuronFile class
 *  This is necessary because Matlab cannot pass/return
 *  variables by reference. Consequently, this class is needed
 *  to interface with the neurons file io code in MATLAB.
 * 
 *  Added to Vision 3/09 by tamachado at request of shlens
 * 
 * @author tamachado@salk.edu (2/4/08)
 *
 */
public class ReadNeurons {
    
    // Parameters for file header
    public final static int HEADER_SIZE = VisionParams.NEURONS_HEADER_CAPACITY;
    public final static int MIN_SPIKES  = 100;
    public final static int MAX_CONTAM  = 1;
    
    // Neurons file
    private NeuronFile nFile;
    
    
    /**
     * Create/open a neurons file to write to.
     * 
     * @param path     path to existing or new neurons file
     * @param prjFile  projections file object
     */
    public ReadNeurons(String path, ProjectionsFile prjFile) {
        
        // Check if a neurons file already exists at path
        File f = new File(path);
        
        // If it exists, then open it 
        if (f.isFile()) {
            try {
                this.nFile = new NeuronFile(path);
            } catch (IOException e) {
                e.printStackTrace();
            }  
       // Otherwise, create a new file     
        } else {
            // Update the given header to work in the neurons file
            VisionHeader header = prjFile.getHeader();
            header.version = NeuronFile.INT_VERSION;
            header.minNeuronSpikes = MIN_SPIKES;
            header.maxContamination = 1;

            // Construct the new neurons file
            try {
                this.nFile = new NeuronFile(path, header, HEADER_SIZE, prjFile.getTTLTimes());
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }
    
    /**
     * This works the same as the above constructor except gets the vision header
     * from another neurons file instead of from a projections file
     * 
     * @param path     path to existing or new neurons file
     * @param nFile    neurons file object
     */
   public ReadNeurons(String path, NeuronFile nFile) {
        
        // Check if a neurons file already exists at path
        File f = new File(path);
        
        // If it exists, then open it 
        if (f.isFile()) {
            try {
                this.nFile = new NeuronFile(path);
            } catch (IOException e) {
                e.printStackTrace();
            }  
       // Otherwise, create a new file     
        } else {
            // Update the given header to work in the neurons file
            VisionHeader header = nFile.getHeader();
            header.version = NeuronFile.INT_VERSION;
            header.minNeuronSpikes = MIN_SPIKES;
            header.maxContamination = 1;

            // Construct the new neurons file
            try {
                this.nFile = new NeuronFile(path, header, HEADER_SIZE, nFile.getTTLTimes());
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }
    
    
    /**
     * Write an individual cluster to the neurons file.
     * 
     * @param electrode    electrode where cluster was found
     * @param clusterIndex nth cluster on electrode we are writing
     * @param spikeTimes   list of spike times for this cluster
     * @param nSpikes      total number of spikes for this cluster
     * @return             boolean on success/failure
     */
    public boolean addNeuron(int electrode, int clusterIndex, int[] spikeTimes, int nSpikes) {
        
        // Get Neuron ID
        int neuronID = NeuronFile.getNeuronID(electrode, clusterIndex);
        
        //Add this neuron to the file
        return this.addNeuron(electrode, spikeTimes, neuronID, nSpikes);
    }
    
    
    /**
     * Write an individual cluster to the neurons file.
     * 
     * @param electrode    electrode where cluster was found
     * @param spikeTimes   list of spike times for this cluster
     * @param neuronID     id of neuron that we are writing
     * @param nSpikes      total number of spikes for this cluster
     * @return             boolean on success/failure
     */
    public boolean addNeuron(int electrode, int[] spikeTimes, int neuronID, int nSpikes) {
        // Add this neuron to the file
        try {
            nFile.addNeuron(electrode, neuronID, spikeTimes, nSpikes);
            return true;
        } catch (IOException e) {
            e.printStackTrace();
            return false;
        }
    }
    
    
    /**
     * Close the IO stream to the neurons file
     * 
     * The Finalize method for FileOutputStream ensures that close has been
     * called before the object is disposed.
     *
     */
    public void closeNeurons() {
        try {
           nFile.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
