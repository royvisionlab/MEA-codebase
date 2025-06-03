package edu.ucsc.neurobiology.vision.tasks;

import java.awt.event.*;

import javax.swing.*;

import com.martiansoftware.jsap.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.parameters.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * A class that can be run from the command line to perform the whole analysis of a raw
 * noise dataset (pink noise, natural power noise, etc).
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class RawMovieDatasetAnalysis
    extends AbstractAction {


    public RawMovieDatasetAnalysis() {
        super("Raw Movie Dataset Analysis");
    }


    public void actionPerformed(ActionEvent event) {
        Vision app = Vision.getInstance();

        ParametersTable table = app.getConfig().showDialog(
                "RawMovieDatasetAnalysis", "Raw Movie Dataset Analysis", app.getMainFrame());
        if (table == null) {
            return;
        }

        final String rawDataFile = table.getFileParameter("Raw Data File");
        final String outputFolder = table.getFileParameter("Output Folder");
        final boolean createSubfolders = table.getBooleanParameter("Create Subfolders");
        final String rawMovieFile = table.getFileParameter("Raw Movie File");
        final String defaultsFile = table.getFileParameter("Defaults File");
        final boolean unwhitenedCovariances = table.getBooleanParameter("Unwhitened Covariances");
        final boolean doEI = table.getBooleanParameter("Calculate EIs");

        Thread runner = new Thread() {
            public void run() {
                try {
                    RunScript.rawMovieDatasetAnalysis(
                            rawDataFile, outputFolder, createSubfolders, rawMovieFile,
                            new Config(defaultsFile), unwhitenedCovariances, doEI);
                } catch (Exception e) {
                    Vision.reportException(e);
                }
            }
        };
        runner.start();
    }


    public static void main(String[] args) throws Exception {
        long t1 = System.currentTimeMillis();


        
        VisionJSAP jsap = new VisionJSAP( 
                RawMovieDatasetAnalysis.class.getName(), 
                new com.martiansoftware.jsap.Parameter[] {
                        new UnflaggedOption("raw", JSAP.STRING_PARSER, JSAP.REQUIRED, "Raw data file."),
                        new UnflaggedOption("root", JSAP.STRING_PARSER, JSAP.REQUIRED, "Root analysis folder."),
                        new UnflaggedOption("rawMovie", JSAP.STRING_PARSER, JSAP.REQUIRED, "Raw movie file. (.rawMovie)"),
                        new Switch("ei", 'e', "ei", "Calculate electrophysiological images?"), 
                        new Switch("genSubs", 'g', "gsub", "Generate subfolders?"),
                        new Switch("unCov", 'u', "uCov", "Unwhitened Covariances?"),
                        new FlaggedOption( "config", JSAP.STRING_PARSER, "config.xml", JSAP.NOT_REQUIRED, 'c', "config", 
                        "Configuration file" )
                }
        );

        JSAPResult parsedArgs = jsap.parse(args);


        RunScript.rawMovieDatasetAnalysis(
                parsedArgs.getString("raw"),
                parsedArgs.getString("root"),
                parsedArgs.getBoolean("genSubs"),
                parsedArgs.getString("rawMovie"),
                new Config(parsedArgs.getString("config")),
                parsedArgs.getBoolean("unCov"),
                parsedArgs.getBoolean("ei")
                );

        RunScript.printAnalysisTime(t1, System.currentTimeMillis());
        System.exit(1);
    }

}
