package edu.ucsc.neurobiology.vision.snf;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;

import edu.ucsc.neurobiology.vision.Vision;
import edu.ucsc.neurobiology.vision.calculations.AbstractCalculation;
import edu.ucsc.neurobiology.vision.io.RawDataFile;
import edu.ucsc.neurobiology.vision.tasks.RunScript;
import edu.ucsc.neurobiology.vision.util.StringUtil;
import edu.ucsc.neurobiology.vision.util.VisionParams;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 * @author Matthew Grivich, The Salk Institute
 */
public class SNFPreprocessing extends AbstractCalculation {
    private String rawDataFileName = null, datasetFolder = null;
    boolean filter = false;

//    public static void main(String[] args) {
//    	Config config = null;
//    	try {
//        config = new Config("config.xml");
//    	} catch(Exception e) {
//    		e.printStackTrace();
//    	}
//        String group = "SNFPreprocessing";
//        ParametersTable t = config.getParameterGroup(group);
//
//
//        if (args.length == t.getParametersCount()) {
//            rawDataFileName = args[0];
//            datasetFolder = args[1];
//            filter = Boolean.valueOf(args[2]).booleanValue();
//        } else {
//            System.out.println(
//                "Incorrect number of command line arguments: " + t.getParametersCount() +
//                " required: \n- rawDataFileName \n- outputFolder \n- filter? (true/false)");
//            System.exit(1);
//
//
//        }
//    }
    
    public void startCalculation() throws Exception {
        String datasetName = new File(datasetFolder).getName();

        // create the folder
        File f = new File(datasetFolder);
        f.mkdirs();

        // generate a folder structure file if needed
        String[] rawDataSources = StringUtil.decomposeString(rawDataFileName, ";");
        if (rawDataSources.length > 1) {
            try {
                PrintWriter w = new PrintWriter(new FileWriter(
                    datasetFolder + File.separator + "folder-structure.txt"));
                for (int i = 0; i < rawDataSources.length; i++) {
                    String name = StringUtil.removeExtension(
                        new File(rawDataSources[i]).getName());
                    RawDataFile rawFile = new RawDataFile(new File(rawDataSources[i]));
                    int nSamples = rawFile.getHeader().getNumberOfSamples();

                    w.println(name + ": " + nSamples);
                }
                w.close();
            } catch (IOException e) {
                Vision.reportFatalException("", e);
            }
        }

        // do the preprocessing
        String filteredRawDataFile;
        if (filter) {
            filteredRawDataFile =
                datasetFolder + File.separator + datasetName +
                VisionParams.BIN_FILE_EXTENSION_512;
            RunScript.rawDataFiltering(rawDataFileName, datasetFolder);
        } else {
            filteredRawDataFile = rawDataFileName;
        }

        RunScript.rawDataNoiseEvaluation(filteredRawDataFile, datasetFolder);

        String noiseFileName =
            datasetFolder + File.separator + datasetName +
            VisionParams.NOISE_FILE_EXTENSION;
        RunScript.convertToElectrodeMajor(filteredRawDataFile, datasetFolder,
                                          noiseFileName);

        // delete the filtered file after the preprocessing
        if (filter) {
            new File(filteredRawDataFile).delete();
        }

        
    }
    
    public void setParameters(HashMap<String, String> p) {
        rawDataFileName = p.get("rawDataFileName");
        datasetFolder = p.get("outputFolder");
        filter = Boolean.parseBoolean(p.get("Filter the Raw Data?"));
    }

}
