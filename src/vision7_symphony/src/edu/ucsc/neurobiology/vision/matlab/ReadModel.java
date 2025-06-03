package edu.ucsc.neurobiology.vision.matlab;

import java.io.File;
import java.io.IOException;

import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.electrodemap.*;

/**
 * These methods can be used by Matlab to write out a 
 * clustering model file after it has been modified externally.
 * This allows output from an arbitrary clustering algorithm
 * to be examined using Cell-Finder, and also makes it 
 * possible to use the standard mapping scripts on the new
 * model file.
 * 
 * Much of this code is based on code from Cell-Finder
 * that was written by Jon Shlens.
 * 
 * Added to Vision 3/09 by tamachado at request of shlens
 * 
 * @author tamachado@salk.edu (4/2/08)
 * 
 *
 */
public class ReadModel {
    
    // Store the actual model file
    private ClusteringModelFile modelFile;
    private VisionHeader modelHeader;
    private ClusteringModelFile.EMModel[] extractions;
    private boolean error = false;
    
    private int nElectrodes;
    
    // Maximum number of slots (clusters) within a model file
    public final static int maxModelSlots = 5000;

    /**
     * Create/open a model file to write to.
     * 
     * @param path      path to existing or new model file
     * @param existPath path to existing model file containing eigenvectors
     */
    public ReadModel(String path, String existPath) {
        
        // See if model file already exists
        File f = new File(path);

        // If it does, notify the user and stop
        if (f.exists()) {
            System.err.println("Model file already exists. Stopping...");
            error = true;
            return;
        }
            
        
        // If we don't have a (partial) model file containing eigenvectors, exit with an error
        // because we can't get them from anywhere else.
        ClusteringModelFile partialModel;
        VisionHeader h;
        
        try {
            partialModel = new ClusteringModelFile(existPath);

            // Get the header from the old model file
            h = partialModel.getUserHeader();
            
            // Save number of electrodes: we must add all the electrodes at once
            // writing to a partial model file is not currently supported
            ElectrodeMap map = ElectrodeMapFactory.getElectrodeMap(h.arrayID);
            nElectrodes = map.getNumberOfElectrodes() + 1;

            // Obtain necessary fields from the old model file:
            // model.electrodes[]     --- neighbors of this electrode
            // model.threshold        --- sigma used in neuron finding
            // model.eigenvectors[][] --- from PCA
            //
            // MUST be present in the old model file
            extractions = new ClusteringModelFile.EMModel[nElectrodes];
            for (int i = 1; i < nElectrodes; i++) {
                extractions[i] = new ClusteringModelFile.EMModel();
            }
            
            // Save models (skip TTL pulses)
            for (int electrode = 1; electrode < nElectrodes; electrode++) {
                try {
                    extractions[electrode] = (ClusteringModelFile.EMModel) partialModel.getNeuronExtraction(electrode);
                } catch (IOException ex) {
                    ex.printStackTrace();
                }
            }
            
            modelHeader = new VisionHeader(h, 1);

            // Create unused, filler quantities (these values are never used)     
            modelHeader.nlPointsEI = 0;
            modelHeader.nrPointsEI = 0;
            modelHeader.minimizationError = 0;
            modelHeader.maxElectrodePatternSize = 0;
            modelHeader.minThreshold = 0;

            // This indicates that MATLAB was used to create this model file.
            modelHeader.unused3 = ("MATLAB").hashCode();

            partialModel.close();

        } catch (IOException e) {
            System.err.println("Model file provided is invalid! Check path to OLD model file.");
            e.printStackTrace();
            error = true;
            return;
        }

        // Create the model file at the path provided
        try {
            modelFile = new ClusteringModelFile(path, maxModelSlots, modelHeader);
        } catch (IOException e) {
            System.err.println("The model file (.model) file cannot be created.");
        }
    }
    

    /**
     * Create a clusterModel object with information about the clusters on electrode e.
     * Then add that clusterModel to the model file. This function is be called
     * from MATLAB to add data to the model file.
     * 
     * @param cProb  Mixture weights for all gaussians on electrode e
     * @param cMeans Mean value of all gaussians in nDimensional space
     * @param cCov   Diagonal of the covariance matrix for each gaussian
     * @param cellId     First neuron id to be used on this electrode
     * @param e      Electrode where the cluster model was obtained from
     * 
     * @return boolean indicating success or failure
     */    
    public boolean addElectrode(double[] cProb, double[][] cMeans, double[][] cCov, int e, int nClusters) {

        // Check if we are trying to add a bad electrode
        if(extractions[e] == null) {
            System.err.println("Skipping electrode " + e + " because it is bad.");
            error = true;
            return false;
        }
        
        // We use "EMModel" because the we have a mixture of gaussians
        // the only other option is "ManualModel" which is used for manual boxing.
        ClusteringModelFile.EMModel model = extractions[e];

        model.extractionID = e;
        
        // Populate the model
        model.cleaningLevel = 0;
        
        model.nClusters = nClusters;
        model.nGaussians = nClusters;
        
        model.neuronIndex = new int[model.nClusters];
        model.neuronID = new int[model.nClusters];
        
        //First neuron id on this electrode
        int id = NeuronFile.getNeuronID(e, 0);
        
        for (int i = 0; i < model.nClusters; i++) {
            model.neuronIndex[i] = i;
            model.neuronID[i] = id++;
        }
        
        // Since nDimensions might be less than what is stored in the prj file
        // we need to compute it here as opposed to just assuming that it is the same
        // as the projections file. Matlab doesn't allow "ragged" arrays: multidimensional
        // arrays are a fixed number of rows x columns, a given row cannot be a random
        // length. Consequently, cMeans[1].length gives us nDimensions.
        model.nDimensions = cMeans[1].length;
        
        //System.err.println("clusters " + model.nClusters + " nDims " + model.nDimensions);
        
        //Weight for each cluster in the mixture
        model.probability = cProb;
        
        //Mean of each cluster, in each dimension
        model.means = cMeans;
        
        //Covariance of each cluster, in each dimension
        model.covariances = cCov;

        // Add this electrode to the model file
        writeElectrode(model, e);
        
        return true;
    }
    
    
    /**
     * Add all clusters from an electrode to the model file contained in this object.
     * 
     * @param clusterModel Contains model information about all clusters on electrode
     * @param e    Electrode where the cluster model was obtained from
     */
    private void writeElectrode(ClusteringModelFile.EMModel clusterModel, int e) {
        try {
            if (clusterModel != null) {
                modelFile.addExtraction(clusterModel);
            } else {
                System.err.println("Cluster model for electrode " + e + " is null.");
            }
        } catch (IOException ex) {
            System.err.println("Error writing to .model file.");
        }
    }
    
    /**
     * Close the IO stream to the model file
     * 
     * The Finalize method for FileOutputStream ensures that close has been
     * called before the object is disposed.
     *
     */
    public void closeModel() {
        try {
           modelFile.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Return error state
     *
     * @return boolean error state, true indicates an error has occurred
     */
    public boolean getError(){
        return error;
    }
    
    
}
