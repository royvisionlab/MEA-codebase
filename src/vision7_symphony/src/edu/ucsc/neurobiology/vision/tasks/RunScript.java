package edu.ucsc.neurobiology.vision.tasks;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.io.RawDataFile;
import static edu.ucsc.neurobiology.vision.anf.SpikeFinding.*;
import static edu.ucsc.neurobiology.vision.anf.SpikeFinding.AnalysisToDo.*;
import edu.ucsc.neurobiology.vision.calculations.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * A utility class used in the tasks.package.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 * @owner Matthew Grivich, The Salk Institute
 */
public class RunScript {

    //////////////////////////////////////////////////////////////////////////////////
    //           DO NOT MODIFY THE BELOW METHODS, JUST USE THEM                    //
    //////////////////////////////////////////////////////////////////////////////////

    public static void convertToElectrodeMajor(
                                               String rawDataSorce, String outputFolder, String sigmas) {

        HashMap<String, String> p = new HashMap<String, String>();
        p.put("Save_To_File", outputFolder);
        p.put("Sigmas_File", sigmas);
        p.put("Raw_Data_File", rawDataSorce);
        CalculationManager c = Vision.getInstance().getCalculationManager();
        c.runCalculation("Convert Raw Data To Electrode Major Format", p);
    }


    public static void spikeFinding(
                                    String rawDataFile, String outputPath, String sigma, AnalysisToDo analysisToDo,
                                    boolean unWhitenedCovariances, Config defaults) {

        String cName = "Spike Finding";
        HashMap<String, String> d = defaults.getParameterList(cName);

        HashMap<String, String> p = new HashMap<String, String>();
        p.put("Raw_Data_Source", rawDataFile);
        p.put("Buffer Size (Kb)", d.get("Buffer Size (Kb)"));
        p.put("Buffers Count", d.get("Buffers Count"));
        p.put("Spike Threshold", d.get("Spike Threshold"));
        p.put("Sigma", "" + sigma);
        p.put("TTL Threshold", d.get("TTL Threshold"));
        p.put("Mean Time Constant", d.get("Mean Time Constant"));
        if (d.containsKey("waitForData")) p.put("waitForData", d.get("waitForData"));

        p.put("Diagnostic Plots", "false");
        p.put("saveRawData", "false");

        p.put("Set Electrodes", d.get("Set Electrodes"));
        p.put("Set Electrodes.arrayID", d.get("Set Electrodes.arrayID"));
        p.put("Set Electrodes.arrayPart", d.get("Set Electrodes.arrayPart"));
        p.put("Set Electrodes.arrayNParts", d.get("Set Electrodes.arrayNParts"));
        p.put("Set Electrodes.flipX", d.get("Set Electrodes.flipX"));
        p.put("Set Electrodes.flipY", d.get("Set Electrodes.flipY"));


        if (analysisToDo == AnalysisToDo.DO_NOTHING) {
            p.put("Analysis", "false");
        } else {
            p.put("Analysis", "true");
            p.put("Analysis.Analysis To Do", "" + analysisToDo.ordinal());
            p.put("Analysis.Left Points", d.get("Analysis.Left Points"));
            p.put("Analysis.Right Points", d.get("Analysis.Right Points"));
            p.put("Analysis.Minimization Error", d.get("Analysis.Minimization Error"));
            p.put("Analysis.Spike To Use", d.get("Analysis.Spike To Use"));
            if (!unWhitenedCovariances) {
                p.put("Analysis.Minimum Noise Events", d.get("Analysis.Minimum Noise Events"));
            } else {
                p.put("Analysis.Minimum Noise Events", 0 + "");
            }
            p.put("Analysis.Electrode Usage", "" + d.get("Analysis.Electrode Usage"));
            p.put("Analysis.Output_Path", outputPath);
        }

        Vision.getInstance().getCalculationManager().runCalculation(cName, p);
    }


    public static void neuronFinding(
                                     String rawDataFile, Config config, String outputPath, boolean unwhitenedCovariances, boolean doEI) throws
                                         IOException {

        rawDataNoiseEvaluation(rawDataFile, outputPath);

        String datasetName = new File(outputPath).getName();
        String rawDataNoiseFileName = outputPath + File.separator + datasetName + VisionParams.NOISE_FILE_EXTENSION;
        spikeFinding(rawDataFile, outputPath, rawDataNoiseFileName,     SAVE_SPIKES_AND_COVARINCES, unwhitenedCovariances, config);

        CalculationManager calcManager = Vision.getInstance().getCalculationManager();
        HashMap<String, String> p = new HashMap<String, String>();

        String cName;
        HashMap<String, String> c;
                
        // ------------------------------------------------------------------
        // INTERRUPTED STREAM PATCHING
        // Purpose of the hack below is to recover available data from interrupted
        // streaming analysis.
        // Error handling is done in SpikeBuffer.finishSampleProcessing()
        // But here we need to fix that the number of samples is now invalid
        // (important further for contamination, neuronViewer, etc).
        //
        // spike, cov, ncov file nSamples are still on the initially expected.
        // We correct that after spikefinding and before anything else
        // Shouldn't affect behavior of anything else (Maybe?).
        // Vincent Deo - Stanford University - 10/15/2015
        if (config.getParameterList("Spike Finding").containsKey("waitForData") &&
            Boolean.parseBoolean(config.getParameterList("Spike Finding").get("waitForData"))) {
            RawDataFile rawData = new RawDataFile(rawDataFile);
            int nSamples = rawData.getHeader().getNumberOfSamples();
            rawData.close();

            String[] exts = { VisionParams.SPIKES_FILE_EXTENSION, ".cov", ".ncov" };
            for (String ext : exts) {
                String filePath = outputPath + File.separator + datasetName + ext;
                RandomAccessFile file = new RandomAccessFile(filePath, "rw");
                // Because nSamples is in bytes 21-24, so we need seek offset 20.
                file.seek(20);
                file.writeInt(nSamples);
                file.close();
            }
        }
        // -----------------------------------------------------------------
                
        if (!unwhitenedCovariances) {
            cName = "Noise Whitened Covariances";
            c = config.getParameterList(cName);
            p.clear();
            p.put("Raw Data Path", rawDataFile);
            p.put("Dataset Folder", outputPath);
            p.put("Number of Threads", c.get("Number of Threads"));

            calcManager.runCalculation(cName, p);
        }

        cName = "PCA Neuron Finding: Projections";
        c = config.getParameterList(cName);
        p.clear();
        p.put("Raw_Data_File", rawDataFile);
        p.put("Dataset Folder", outputPath);
        p.put("PCA Dimensions", c.get("PCA Dimensions"));
        if (c.containsKey("nThreads")) {
            p.put("nThreads", c.get("nThreads"));
        } else {
            p.put("nThreads", "3");
        }
        calcManager.runCalculation(cName, p);

        cName = "PCA Neuron Finding: Clustering";
        c = config.getParameterList(cName);
        p.clear();
        p.put("Dataset_Folder", outputPath);
        p.put("From", "1");
        p.put("To", "-1");
        p.put("Generate Report", c.get("Generate Report"));
        p.put("Bins Per Dimension", c.get("Bins Per Dimension"));
        p.put("Clustering Significance", c.get("Clustering Significance"));
        p.put("Minimum Clusters", c.get("Minimum Clusters"));
        p.put("Miximum Clusters", c.get("Miximum Clusters"));
        p.put("Spikes Used For EM", c.get("Spikes Used For EM"));
        p.put("Density Clustering Spike Loss", c.get("Density Clustering Spike Loss"));
        p.put("Min EM Iterations", c.get("Min EM Iterations"));
        p.put("Max EM Iterations", c.get("Max EM Iterations"));
        p.put("EM Likelihood Delta", c.get("EM Likelihood Delta"));
        p.put("Clustering Threads", c.get("Clustering Threads"));
        calcManager.runCalculation(cName, p);

        cName = "Neuron Cleaning";
        c = config.getParameterList(cName);
        p.clear();
        p.put("Neuron_File", outputPath + File.separator + datasetName + ".neurons-raw");
        p.put("Minimun Number of Spikes", c.get("Minimun Number of Spikes"));
        p.put("Maximum Contamination", c.get("Maximum Contamination"));
        p.put("Coincidence Time", c.get("Coincidence Time"));
        p.put("Maximum Correlation", c.get("Maximum Correlation"));
        calcManager.runCalculation("Neuron Cleaning", p);

        if (doEI) {
            cName = "Electrophysiological Imaging Fast";
            c = config.getParameterList(cName);
            p.clear();
            p.put("Dataset Folder", outputPath);
            p.put("Raw Data File", rawDataFile);
            p.put("Mean Time Constant", "0.01");
            p.put("Left Samples",      c.get("Left Samples"));
            p.put("Right Samples",     c.get("Right Samples"));
            p.put("Spikes To Average", c.get("Spikes To Average"));
            p.put("Threads",           c.get("Threads"));
            calcManager.runCalculation(cName, p);
        }
    }



    public static final void mappingAnalysis(String masterFolder, String datasetFolder, String rawDataFile, Config defaults, boolean projectionsOnly) { 
        String datasetName = new File(datasetFolder).getName();
        String sigmaFileName = datasetFolder + File.separator + datasetName + ".noise";
                
        new File(datasetFolder).mkdirs();
                
        rawDataNoiseEvaluation(rawDataFile, datasetFolder);
        spikeFinding(rawDataFile, datasetFolder, sigmaFileName, SAVE_SPIKES, false, defaults);
                
        HashMap<String, String> p = new HashMap<String, String>();
        p.put("Master Dataset Folder", masterFolder);
        p.put("Dataset Folder", datasetFolder);
        p.put("Raw Data File", rawDataFile);
        p.put("Projections Only", Boolean.toString(projectionsOnly));
        Vision.getInstance().getCalculationManager().runCalculation("PCA Neuron Finding: Mapping", p);
    }

        
    public static void rawDataFiltering(String rawDataSorce, String outputFolder) {
        HashMap<String, String> p = new HashMap<String, String>();
        p.put("Output Folder", outputFolder);
        p.put("Raw_Data_File", rawDataSorce);
        p.put("callibrateGains", "" + false);
        CalculationManager c = Vision.getInstance().getCalculationManager();
        c.runCalculation("Raw Data Filtering", p);
    }


    public static void rawDataNoiseEvaluation(String rawDataSorce, String outputFolder) {
        HashMap<String, String> p = new HashMap<String, String>();
        p.put("Raw_Data_File", rawDataSorce);
        p.put("Save_To_File", outputFolder);
        CalculationManager c = Vision.getInstance().getCalculationManager();
        c.runCalculation("Raw Data Noise Evaluation", p);
    }

        
    public static void calcAuxParams(
                                     Config movieXML, String datasetPath) throws IOException {
        HashMap<String, String> map = movieXML.getParameterList("Calculate Auxiliary Parameters");

        //              String datasetName = new File(datasetPath).getName();

        map.put("File_Path", datasetPath);

        Vision.getInstance().getCalculationManager().runCalculation(
                                                                    "Calculate Auxiliary Parameters", map);
    }

    public static void createWhiteNoiseMovie(
                                             Config movieXML, String datasetPath) throws IOException {

        HashMap<String, String> map = movieXML.getParameterList("Make White Noise Movie");
        //              String datasetName = new File(datasetPath).getName();
        map.put("filePath", datasetPath + File.separator);
        Vision.getInstance().getCalculationManager().runCalculation(
                                                                    "Make White Noise Movie", map);
    }


    /**
     * Do not modify this method since it is critical for the work of the
     * fast diagnostics.
     *
     * @param rawDataFileName String
     * @param analysisRootFolder String
     * @param createSubfolders boolean
     * @param movieDefinitionFile String
     * @param defaults Config
     * @param unWhitenedCovariances boolean
     * @param doEI boolean
     * @throws IOException
     */
    public static void whiteNoiseDatasetAnalysis(
                                                 String rawDataFileName, String analysisRootFolder, boolean createSubfolders, Config movieXML,
                                                 Config config, boolean unWhitenedCovariances, boolean doEI) throws IOException {
        // complete the raw data name

        // fix output folder
        File f = new File(rawDataFileName);
                
        if (createSubfolders) {
            analysisRootFolder += File.separator + f.getParentFile().getName() +
                File.separator + f.getName();
            File toMake = new File(analysisRootFolder);
            toMake.mkdirs();
        }
                

        // do the rest
        //              String datasetName = new File(analysisRootFolder).getName();
        neuronFinding(rawDataFileName, config, analysisRootFolder, unWhitenedCovariances, doEI);

        try {
            RunScript.createWhiteNoiseMovie(movieXML, analysisRootFolder);
        } catch (IOException e) {
            throw new IOException("The loading of the white noise movie file failed.");
        }
        RunScript.calcAuxParams(movieXML, analysisRootFolder);
        staCalculation(analysisRootFolder, "", config);
        makeParametersFile(analysisRootFolder, config, true, true, doEI);
    }

        
    /**
     * Do not modify this method since it is critical for the work of the
     * fast diagnostics.
     *
     * @param analysisRootFolder String
     * @param createSubfolders boolean
     * @param movieDefinitionFile String
     * @param defaults Config
     * @throws IOException
     */
    public static void staAnalysis(
                                   String analysisRootFolder, Config movieXML, Config config) throws IOException {

        //Determine the name of various files that may exist in the root folder
        String fName = new File(analysisRootFolder).getName();
        File neurons = new File(analysisRootFolder + 
                                File.separator + fName + VisionParams.NEURON_FILE_EXTENSION);
        File globals = new File(analysisRootFolder + 
                                File.separator + fName + VisionParams.GLOBALS_FILE_EXTENSION);
        File ei = new File(analysisRootFolder + 
                           File.separator + fName + VisionParams.EI_FILE_EXTENSION);
        
        boolean existEi = false;
        
        //Stop this computation if a neurons file does not exist
        if (!neurons.canRead())
            throw new IOException("Cannot read input neurons file!");
        
        //Create an empty globals file if necessary
        if (!globals.canRead())
            throw new IOException("Cannot read input globals file!");
        
        //If eis already exist, tell the parameters file computation
        if (ei.canRead())
            existEi = true;
        
        try {
            RunScript.createWhiteNoiseMovie(movieXML, analysisRootFolder);
        } catch (IOException e) {
            throw new IOException("The loading of the white noise movie file failed.");
        }
        
        RunScript.calcAuxParams(movieXML, analysisRootFolder);
        staCalculation(analysisRootFolder, "", config);
        makeParametersFile(analysisRootFolder, config, true, true, existEi);
    }


    public static void rawMovieDatasetAnalysis(
                                               String rawDataFileName, String analysisRootFolder, boolean createSubfolders,
                                               String rawMovieFile, Config defaults, boolean unwhitenedCovariances, boolean doEI) throws IOException {



        // fix output folder
        File f = new File(rawDataFileName);
        if(createSubfolders) {
            analysisRootFolder += File.separator + f.getParentFile().getName() + File.separator + f.getName();
            File toMake = new File(analysisRootFolder);
            toMake.mkdirs();
        }

        //do the rest
        String datasetName = new File(analysisRootFolder).getName();


        new File(analysisRootFolder).mkdirs();


        neuronFinding(rawDataFileName, defaults, analysisRootFolder, unwhitenedCovariances, doEI);

        createRunTimeParamsFile(analysisRootFolder, 5.8, 0.0, 0.0, 1);

        staCalculation(analysisRootFolder, rawMovieFile, defaults);
        makeParametersFile(analysisRootFolder, defaults, true, true, doEI);
    }




    public static void createRunTimeParamsFile(
                                               String folder, double pixelSize, double xOffset, double yOffset,
                                               int refreshInterval) {

        HashMap<String, String> p = new HashMap<String, String>();
        p.put("File_Path", folder);
        p.put("PixelSize", "" + pixelSize);
        p.put("xOffset", "" + xOffset);
        p.put("yOffset", "" + yOffset);
        p.put("refreshInterval", "" + refreshInterval);
        Vision.getInstance().getCalculationManager().runCalculation(
                                                                    "Create Run Time Parameters File", p);
    }


    public static void staCalculation(String folder, String moviePath, Config defaults) {
        String cName = "STA Calculation";
        HashMap<String, String> d = defaults.getParameterList(cName);

        HashMap<String, String> p = new HashMap<String, String>();
        p.put("File_Path", folder);
        p.put("Raw_Movie_Path", moviePath);
        p.put("STA Depth", d.get("STA Depth"));
        p.put("STA Offset", d.get("STA Offset"));
        p.put("Spikes To Calculate", d.get("Spikes To Calculate"));
        p.put("STAs Calculated At Once", d.get("STAs Calculated At Once"));
        p.put("Double Threaded", d.get("Double Threaded"));
        p.put("Calculate STV", d.get("Calculate STV"));
                
        Vision.getInstance().getCalculationManager().runCalculation(cName, p);
    }

    // Just use the config-xml; no funny stuff!
    // Okay, actually it does require overriding the MainFilePath, but that's the only funny stuff!
    public static void makeParametersFile(String folder, Config defaults) {
        String cName = "Make Parameters File";
        HashMap<String, String> d = defaults.getParameterList(cName);
        d.put("MainFilePath", folder);
        Vision.getInstance().getCalculationManager().runCalculation(cName, d);
    }
        
    public static void makeParametersFile(
                                          String folder, Config defaults, boolean fitSTAs, boolean getTimeCourse, boolean fitEIs) {

        String cName = "Make Parameters File";
        HashMap<String, String> d = defaults.getParameterList(cName);
        LinkedHashMap<String, String> p = new LinkedHashMap<String, String>();

        p.put("MainFilePath", folder);
        p.put("nThreads", d.get("nThreads"));

        if (fitSTAs) {
            p.put("STAFitCalculator", "true");
            p.put("STAFitCalculator.fitContours", "false");
            if (d.containsKey("STAFitCalculator.downsampleFactor")) {
                p.put("STAFitCalculator.downsampleFactor", d.get("STAFitCalculator.downsampleFactor"));
            } else {
                p.put("STAFitCalculator.downsampleFactor", "1");
            }
        }

        if (fitEIs) {
            p.put("EIFitCalculator", "true");
        }

        if (getTimeCourse) {
            p.put("TimeCourseCalculator", d.get("TimeCourseCalculator"));
            p.put("TimeCourseCalculator.significance",
                  d.get("TimeCourseCalculator.significance"));
            p.put("TimeCourseCalculator.nTemporalSubfilters",
                  d.get("TimeCourseCalculator.nTemporalSubfilters"));
            p.put("TimeCourseCalculator.latestZero",
                  d.get("TimeCourseCalculator.latestZero"));
        }

        p.put("AutoCalculator", "true");
        p.put("AutoCalculator.timeRange", "100"); // ms
        p.put("AutoCalculator.binning", "0.5");

        p.put("ContaminationIndexCalculator", "true");

        Vision.getInstance().getCalculationManager().runCalculation(cName, p);
    }


    public static HashMap<String, String> loadMovieDefinitionFile(String movieDefinitionFile) throws
        IOException {

        LineNumberReader r = new LineNumberReader(new FileReader(movieDefinitionFile));
        HashMap<String, String> map = new HashMap<String, String>();
        String s;
        while ( (s = r.readLine()) != null) {
            int n = s.indexOf("=");
            String name = s.substring(0, n).trim();
            String value = s.substring(n + 1, s.length()).trim();
            map.put(name, value);
        }
        r.close();

        if (!map.get("MovieType").equals("White Noise")) {
            throw new IllegalArgumentException(
                                               "The movie definition file (.mdf) provided does not describe a White Noise Movie");
        }

        return map;
    }
        
    public static void printAnalysisTime(long t1, long t2) {
        System.out.print("\nThe whole analysis took: ");
        double t = (t2 - t1) / 1000.0;
        if (t < 60) {
            System.out.println(t + " sec.");
        } else if (t < 3600) {
            System.out.println(StringUtil.format(t / 60.0, 1) + " min.");
        } else {
            System.out.println(StringUtil.format(t / 3600.0, 1) + " hours.");
        }
    }

}
