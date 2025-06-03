package edu.ucsc.neurobiology.vision.tasks;

import java.io.*;

import java.awt.event.*;
import javax.swing.*;

import com.martiansoftware.jsap.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.io.chunk.*;
import edu.ucsc.neurobiology.vision.parameters.*;
import edu.ucsc.neurobiology.vision.util.*;



/**
 * A class that can be run from the command line to perform the whole neuron finding
 * process on a dataset. No stimulus-specific calculations are done. A working parameters
 * file (.params) is still created.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 * @author Matthew Grivich, The Salk Institue
 */
public class NeuronIdentification extends AbstractAction {

    public NeuronIdentification() {
        super("Neuron Identification");
    }


    public void actionPerformed(ActionEvent event) {
        Vision app = Vision.getInstance();

        final ParametersTable table = app.getConfig().showDialog(
                "Neuron Identification", "Neuron Identification", app.getMainFrame());
        if (table == null) {
            return;
        }

        Thread runner = new Thread() {
            public void run() {
                // read the partams
                String rawDataFileName = table.getFileParameter("Raw Data File");
                final String rootFolder = table.getFileParameter("Output Folder");
                final boolean createSubfolders = table.getBooleanParameter("Create Subfolders");
                final String defaultsFile = table.getFileParameter("Defaults File");
                final boolean unwhitenedCovariances = table.getBooleanParameter("Unwhitened Covariances");
                final boolean doEI = table.getBooleanParameter("Calculate EIs");

                neuronID(rawDataFileName, rootFolder, createSubfolders, defaultsFile, unwhitenedCovariances, doEI);
            }
        };
        runner.start();
    }


    private static void neuronID(
            String rawDataFileName, String rootFolder, boolean createSubfolders, String defaultsFile, boolean unwhitenedCovariances, boolean doEI) {

        try {
            // create the needed folders and prepare the folder names
        //	rawDataFileName = IOUtil.getPhysicalFolderPath(rawDataFileName);

            File f = new File(rawDataFileName);
            if (createSubfolders) {
                rootFolder += File.separator + f.getParentFile().getName() +
                File.separator + f.getName();
            }

            new File(rootFolder).mkdirs();

            // start the calculations
            Config defaults = new Config(defaultsFile);
            RunScript.neuronFinding(rawDataFileName, defaults, rootFolder, unwhitenedCovariances, doEI);

            //Create empty globals file, to prevent crashing later.
            String globalsFileName = rootFolder + File.separator + new File(rootFolder).getName() +	".globals";
            GlobalsFile gf = new GlobalsFile(globalsFileName, ChunkFile.READ);
            System.out.println(globalsFileName);
            gf.close();

            //RunScript.createNoMovieGlobalsFile(defaults, rawDataFileName, rootFolder);
            RunScript.makeParametersFile(rootFolder, defaults, false, false, doEI);
        } catch (Exception e) {
            Vision.reportException(e);
        }
    }



    public static void main(String[] args) throws Exception {
        long t1 = System.currentTimeMillis();

        VisionJSAP jsap = new VisionJSAP( 
                NeuronIdentification.class.getName(), 
                new com.martiansoftware.jsap.Parameter[] {
                    new UnflaggedOption("raw", JSAP.STRING_PARSER, JSAP.REQUIRED, "Raw data file."),
                    new UnflaggedOption("root", JSAP.STRING_PARSER, JSAP.REQUIRED, "Root analysis folder."),
                    new Switch("ei", 'e', "ei", "Calculate electrophysiological images?"), 
                    new Switch("genSubs", 'g', "gsub", "Generate subfolders?"),
                    new Switch("unCov", 'u', "uCov", "Unwhitened Covariances?"),
                    new FlaggedOption( "config", JSAP.STRING_PARSER, "config.xml", JSAP.NOT_REQUIRED, 'c', "config", 
                    "Configuration file" )
                }
        );
        JSAPResult parsedArgs = jsap.parse(args);

        neuronID(
            parsedArgs.getString("raw"),
            parsedArgs.getString("root"),
            parsedArgs.getBoolean("genSubs"),
            parsedArgs.getString("config"),
            parsedArgs.getBoolean("unCov"),
            parsedArgs.getBoolean("ei")
        );

        RunScript.printAnalysisTime(t1, System.currentTimeMillis());
        System.exit(1);
    }

}
