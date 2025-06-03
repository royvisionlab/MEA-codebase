package edu.ucsc.neurobiology.vision.tasks;

import java.awt.event.ActionEvent;
import java.util.Arrays;
import java.util.HashMap;
import java.io.File;
import java.lang.Exception;

import javax.swing.AbstractAction;

import com.martiansoftware.jsap.FlaggedOption;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Switch;
import com.martiansoftware.jsap.UnflaggedOption;

import edu.ucsc.neurobiology.vision.util.VisionJSAP;
import edu.ucsc.neurobiology.vision.util.VisionParams;
import edu.ucsc.neurobiology.vision.Config;
import edu.ucsc.neurobiology.vision.Vision;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.io.ClusteringModelFile.*;
import edu.ucsc.neurobiology.vision.math.ExpectationMaximization1;
import edu.ucsc.neurobiology.vision.parameters.ParametersTable;

/**
 * This class is used for performing mapping on datasets analyzed using
 * noise whitened covariance matrices and works in the following way:
 * 
 * mapping(dataset000, dataset001, ..., datasetN)
 *  concatenate all input datasets and save output .prj, .model, .neurons
 *  for d = dataset000 to datasetN:
 *     project d using eigenvectors from concatenated model file (save d.prj)
 *     apply clusters from concatenated model file to d
 *     re-fit clusters (save d.model, d.neurons)
 *  end
 * end
 * 
 * @see edu.salk.snl.cellfinder.convert.MappingAnalysis
 * @see edu.salk.neurobiology.vision.tasks.MappingAnalysis
 * @author tamachado@salk.edu
 * @since 2009-06-25
 */


public class ConcatenatedMappingAnalysis extends AbstractAction {
    private static final boolean DEBUG = false;
    private static final long serialVersionUID = 2L;
    static HashMap<String, String> params;


    public ConcatenatedMappingAnalysis() {
        super("Concatenated Mapping Analysis");
    }

    public void actionPerformed(ActionEvent event) {
        Vision app = Vision.getInstance();
        app.sendMessage("Concatenated mapping analysis...");

        ParametersTable table = Vision.getConfig().showDialog(
                                                              "Concatenated Mapping Analysis", "Concatenated Mapping Analysis", app.getMainFrame());
        if (table == null) return;

        final String concatenatedDatasetString = table.getFileParameter("Raw Datasets to Map");
        final String rootOutputDirectory = table.getFileParameter("Output Folder");
        final Boolean noRefit = table.getBooleanParameter("Refit Clusters");
        final String configFile = table.getFileParameter("Config File");

        // verify that there were multiple datasets passed in, otherwise there is nothing to map!
        DataFileStringParser dsp = new DataFileStringParser(concatenatedDatasetString);
        if (dsp.getDatasets().length < 2) {
            Exception e = new Exception("You must specify more than one input dataset to do mapping!");
            Vision.reportException(e);
            return;
        }

        Thread runner = new Thread() {
                public void run() {
                    try {

                        // run mapping calculation
                        performMapping(
                                       concatenatedDatasetString, rootOutputDirectory, 
                                       noRefit, new Config(configFile));

                        // clear progress bar
                        if (Vision.isGUIBased()) {
                            Vision.getInstance().sendMessage("");
                            Vision.getInstance().setProgress(0);
                        }
                    } catch (Exception e) {
                        Vision.reportException(e);
                    }
                }
            };
        runner.start();
    }

    public static void performMapping(String concatenatedDatasetString,
                                      String rootOutputDirectory,       boolean noRefitClusters, Config config) {

        // mapping: project spikes using model file of another data set
        System.out.println("Setting up mapping analysis ... ");
        // always use NWPCA (that's the point of this)
        boolean disableWhitening = false;
        // never compute EIs for the concatenated run (call vision-auto-ei-only manually)
        boolean computeEIs = true;
        // parse the input dataset paths
        DataFileStringParser dsp = new DataFileStringParser(concatenatedDatasetString);
        // save path information
                
        // Vincent Deo - Stanford University - 10/14/2015
        // Removing all names because of overly long folder names
        // String allNames = "";
        String[] paths = dsp.getDatasets();
        String[] names = new String[paths.length];
        String[] outputFolders = new String[paths.length];

        // get the name of each dataset the we're processing (e.g. data000)
        for (int dataset = 0; dataset < names.length; dataset++) {
            File current = new File(paths[dataset]);
            names[dataset] = current.getName();
            /*
              if (dataset == 0) {
              allNames = "-from-" + names[dataset];
              } else {
              allNames = allNames + "_" + names[dataset];
              }
            */
        }
                
        try {
            // 1. RUN NEURON IDENTIFICATION ON CONCATENATED RUN
            // ===================================================

            System.out.println("1/3. Analyzing concatenated dataset...");
            if (!DEBUG) {

                // construct output directory if necessary
                File outDir = new File(rootOutputDirectory);
                if (!outDir.mkdir() && !outDir.exists()) 
                    throw new Exception("Could not create output directory!");
                                
                // analyze the concatenated run
                RunScript.neuronFinding(concatenatedDatasetString, config, rootOutputDirectory, disableWhitening, computeEIs);
            }
                        
            // 2. REPROJECT EACH DATASET IN CONCATENATED PC SPACE
            // ===================================================

            System.out.println("2/3. Projecting each dataset in concatenated space...");
            for(int dataset = 0; dataset < paths.length; dataset++) {

                // construct the subfolders that will hold the output for each dataset
                outputFolders[dataset] = rootOutputDirectory + File.separator + names[dataset]; // + allNames;
                File subfolder = new File(outputFolders[dataset]);
                if (!subfolder.mkdir() && !subfolder.exists()) 
                    throw new Exception("Could not create subfolders in output directory!");

                // reproject each dataset in concatenated space  (only save the projections from this call)
                if (!DEBUG) {  
                    boolean prjOnly = true;
                    RunScript.mappingAnalysis(new File(rootOutputDirectory).getAbsolutePath(),
                                              outputFolders[dataset], new File(paths[dataset]).getAbsolutePath(), config, prjOnly);
                }
            }

            // 3. APPLY AND REFIT CLUSTERS FOR EACH DATASET
            // ===================================================

            // load concatenated model file. we are going to fit it to each dataset individually
            ClusteringModelFile concatenatedModel = new ClusteringModelFile(rootOutputDirectory + File.separator +
                                                                            new File(rootOutputDirectory).getName() + VisionParams.MODEL_FILE_EXTENSION);

            // load concatenated neurons file so we can obtain the TTLs and header for use later
            NeuronFile concatenatedNeurons = new NeuronFile(rootOutputDirectory + File.separator +
                                                            new File(rootOutputDirectory).getName() + VisionParams.NEURON_FILE_EXTENSION);

            System.out.println("3/3. Processing individual datasets...");
            for(int dataset = 0; dataset < paths.length; dataset++) {
                // dataset output directory
                String outputFolder = outputFolders[dataset];

                // output files will be named outputName + a file extension
                String outputName = outputFolder + File.separator + new File(outputFolder).getName();

                // get the output directory that contains a projections file generated in step 2
                String projectionsPath = outputName + VisionParams.PROJ_FILE_EXTENSION;
                File pf = new File(projectionsPath);
                if (! pf.canRead()) {
                    throw new Exception("Could not read a projections file: " + pf.getName());
                }
                // open the prj file and spikes file. create the new model and neurons files
                ProjectionsFile newProjections = new ProjectionsFile(outputName + VisionParams.PROJ_FILE_EXTENSION);
                SpikeFile newSpikes = new SpikeFile(outputName + VisionParams.SPIKES_FILE_EXTENSION);
                NeuronFile newNeurons = new NeuronFile(outputName + VisionParams.NEURON_FILE_EXTENSION,
                                                       concatenatedNeurons.getHeader(), concatenatedNeurons.getHeaderCapacity(),
                                                       newSpikes.getTTLTimes());

                ClusteringModelFile newModel = new ClusteringModelFile(outputName + 
                                                                       VisionParams.MODEL_FILE_EXTENSION, VisionParams.MAX_MODEL_SLOTS, concatenatedModel.getUserHeader());

                // save out a model and neurons file for each dataset and refit if necessary
                System.out.println("Processing " + outputName);
                saveClusters(concatenatedModel, newModel, newSpikes, concatenatedNeurons,
                             newNeurons, newProjections, noRefitClusters);

                // close up the files
                newNeurons.close();
                newModel.close();

                // remove noise and spikes files now that we are done with them
                if (!DEBUG) {
                    File sf = new File(outputName + VisionParams.SPIKES_FILE_EXTENSION);
                    File nf = new File(outputName + VisionParams.NOISE_FILE_EXTENSION);
                    sf.delete();
                    nf.delete();
                }
            }

            concatenatedModel.close();

        } catch (Exception e) {
            Vision.reportException(e);
        }
    }

    /**
     * This function loads up clustering and projections information in order to populate new (but empty) neurons and 
     * model files that are passed in. If "refitClusters" is enabled, the clusters from the old model file will be refitted
     * and saved in the new model file (and neurons). Otherwise, the new model file will be identical to the old model file.
     * The spikes file must contain the spike times from the same dataset as the projections file.
     * 
     * @param oldModel         Model file containing clusters that will be refitted
     * @param newModel         Model to populate with refitted clusters
     * @param newSpikes        Spikes file that corresponds to the projections file passed in
     * @param newNeurons       Empty neurons file that will be populated with the clusters
     * @param newProjections   Projections file that the clusters will be fitted to
     * @param refitClusters    Should clusters be refitted?
     */
    private static void saveClusters(ClusteringModelFile oldModel, ClusteringModelFile newModel, 
                                     SpikeFile newSpikes, NeuronFile oldNeurons, NeuronFile newNeurons, ProjectionsFile newProjections, 
                                     boolean refitClusters) {

        // get the pointer to the main application so we can update the progress bar
        Vision app = Vision.getInstance();
        // save the header from the model file
        VisionHeader modelHeader;
        // get a list of valid electrodes
        int[] extractions = oldModel.getExtractionIds();
        // sort them into ascending numerical order
        Arrays.sort(extractions);
        // figure out how many dimensions we're going to use
        int nDimensions = newProjections.getHeader().nDimensions;
        // start the progress bar
        app.startProgressBar();
        // create the data structures to store projections data
        float[][] data = new float[nDimensions][newProjections.maxSpikesPerElectrode];
        int[] times = new int[newProjections.maxSpikesPerElectrode];
        int[] goodNeuronIDs = oldNeurons.getIDList();
        Arrays.sort(goodNeuronIDs);
        
        try {
            modelHeader = newModel.getUserHeader();

            // save out the clusters on each electrode
            for(int e = 0; e < extractions.length; e++) {
                // update the progress bar
                app.setProgress((int) (100.0 * e / (extractions.length - 1)));

                int electrode = extractions[e];

                // get projections and spike times
                int nSpikes = newProjections.readProjections(electrode, data, times);

                // get the model on the current electrode
                EMModel model = (EMModel) oldModel.getNeuronExtraction(electrode);
                //                              int oldNClusters = model.nClusters;

                // load the clusters and data into the clustering algorithm
                ExpectationMaximization1 em = new ExpectationMaximization1(nDimensions,
                                                                           VisionParams.maxClusters, newProjections.maxSpikesPerElectrode);
                em.reset(data, modelHeader.nEMSpikes, nSpikes, model);

                em.setIdentifier("Electrode " + e + ": ");
                // refit the clusters
                if(refitClusters)
                    em.fit(modelHeader.emLikelihoodDelta, modelHeader.minEMIterations, modelHeader.maxEMIterations);

                // update the model file as necessary
                model.nClusters = em.getClustersCount();
                model.nGaussians = em.getGaussiansCount();
                model.probability = new double[model.nGaussians];
                model.means = new double[model.nGaussians][];
                model.covariances = new double[model.nGaussians][];
                int[] clusterCount = new int[model.nClusters];
                for (int i = 0; i < model.nGaussians; i++) {
                    model.probability[i] = em.getGaussianProbability(i);
                    model.means[i] = em.getMeans(i);
                    model.covariances[i] = em.getCovariances(i);
                }

                // figure out which spikes belong to each cluster on this electrode
                int[][] spikeTimes = new int[model.nClusters][newProjections.maxSpikesPerElectrode];
                if(model.nClusters > 0) {
                    for(int spike = 0; spike < nSpikes; spike++) {
                        int cluster = 0;

                        cluster = em.getCluster(spike);
                        spikeTimes[cluster][clusterCount[cluster]] =  times[spike];
                        clusterCount[cluster]++;

                    }
                }

                // add the clusters to the model file for serialization
                newModel.addExtraction(model);
                // get a list of valid cell ids on this electrode
                int[] ids = model.neuronID;
                // add the clusters to the neurons file for serialization
                for(int id = 0; id < em.getClustersCount(); id++) {
                    // truncate the spikeTimes array to be only as long as the number of spikes
                    int[] currentSpikeTimes = new int[clusterCount[id]];
                    System.arraycopy(spikeTimes[id], 0, currentSpikeTimes, 0, clusterCount[id]);
                    // add this cluster to the neurons file if it's good enough
                                                
                    if(Arrays.binarySearch(goodNeuronIDs, ids[id]) >= 0) {
                        newNeurons.addNeuron(electrode, ids[id], spikeTimes[id], clusterCount[id]);
                    }

                }
            }
        } catch (Exception ex) {
            ex.printStackTrace();
            System.out.println("EM failed. " + ex.getMessage() + ". Miximum reached. Dumitru is sorry.");
            app.endProgressBar();
        }
        // we're done
        app.endProgressBar();
    }

    public static void main(String[] args) throws Exception{
        long t1 = System.currentTimeMillis();

        VisionJSAP jsap = new VisionJSAP( 
                                         ConcatenatedMappingAnalysis.class.getName(), 
                                         new com.martiansoftware.jsap.Parameter[] {
                                             new UnflaggedOption("raw", JSAP.STRING_PARSER, JSAP.REQUIRED, "Raw datasets to analyze."),
                                             new UnflaggedOption("root", JSAP.STRING_PARSER, JSAP.REQUIRED, "Analysis output root folder."),
                                             new Switch("noRefit", 'n', "noRefit", "Do not refit clusters."),
                                             new FlaggedOption( "config", JSAP.STRING_PARSER, "config.xml", JSAP.NOT_REQUIRED, 'c', "config", 
                                                                "Configuration file" )
                                         }
                                          );
        JSAPResult parsedArgs = jsap.parse(args);

        // verify that there were multiple datasets passed in, otherwise there is nothing to map!
        DataFileStringParser dsp = new DataFileStringParser(parsedArgs.getString("raw"));
        if (dsp.getDatasets().length < 2) {
            Exception e = new Exception("You must specify more than one input dataset to do mapping!");
            Vision.reportException(e);
            System.exit(1);
        }

        // warn the user if refitting has been disabled
        if (parsedArgs.getBoolean("noRefit"))
            System.out.println("Cluster refitting has been disabled by the user!");
                
        performMapping(
                       parsedArgs.getString("raw"), parsedArgs.getString("root"), 
                       !parsedArgs.getBoolean("noRefit"), new Config(parsedArgs.getString("config")));
                
        RunScript.printAnalysisTime(t1, System.currentTimeMillis());
        System.exit(1);
    }
}
