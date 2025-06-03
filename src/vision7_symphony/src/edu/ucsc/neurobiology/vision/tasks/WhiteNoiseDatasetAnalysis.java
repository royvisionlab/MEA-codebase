package edu.ucsc.neurobiology.vision.tasks;

import java.awt.event.ActionEvent;

import javax.swing.AbstractAction;

import com.martiansoftware.jsap.FlaggedOption;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Switch;
import com.martiansoftware.jsap.UnflaggedOption;

import edu.ucsc.neurobiology.vision.Config;
import edu.ucsc.neurobiology.vision.Vision;
import edu.ucsc.neurobiology.vision.parameters.ParametersTable;
import edu.ucsc.neurobiology.vision.util.VisionJSAP;


/**
 * A class that can be run from the command line to perform the whole analysis of a white
 * noise dataset.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 * @owner Matthew Grivich, The Salk Institute
 */
public class WhiteNoiseDatasetAnalysis
extends AbstractAction {


    public WhiteNoiseDatasetAnalysis() {
        super("White Noise Dataset Analysis");
    }


    public void actionPerformed(ActionEvent event) {
        Vision app = Vision.getInstance();

        ParametersTable table = app.getConfig().showDialog(
                "WhiteNoiseDatasetAnalysis", "White Noise Dataset Analysis", app.getMainFrame());
        if (table == null) {
            return;
        }

        final String rawDataFile            = table.getFileParameter("Raw Data File");
        final String outputFolder           = table.getFileParameter("Output Folder");
        final boolean createSubfolders      = table.getBooleanParameter("Create Subfolders");
        final String movieXMLFile           = table.getFileParameter("Movie XML File");
        final String configFile             = table.getFileParameter("Config File");
        final boolean unwhitenedCovariances = table.getBooleanParameter("Unwhitened Covariances");
        final boolean doEI                  = table.getBooleanParameter("Calculate EIs");

        Thread runner = new Thread() {
            public void run() {
                try {
                    
                    RunScript.whiteNoiseDatasetAnalysis(
                            rawDataFile, outputFolder, createSubfolders, new Config(movieXMLFile),
                            new Config(configFile), unwhitenedCovariances, doEI);
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
            WhiteNoiseDatasetAnalysis.class.getName(), 
            new com.martiansoftware.jsap.Parameter[] {
                new UnflaggedOption("raw",       JSAP.STRING_PARSER, JSAP.REQUIRED, "Raw data file."),
                new UnflaggedOption("root",      JSAP.STRING_PARSER, JSAP.REQUIRED, "Root analysis folder."),
                new UnflaggedOption("movie_XML", JSAP.STRING_PARSER, JSAP.REQUIRED, "Movie xml file."),
                new Switch("ei", 'e', "ei", "Calculate electrophysiological images?"), 
                new Switch("genSubs", 'g', "gsub", "Generate subfolders?"),	
                new Switch("unCov", 'u', "uCov", "Unwhitened Covariances?"),
                new FlaggedOption( "config", JSAP.STRING_PARSER, "config.xml", JSAP.NOT_REQUIRED, 'c', "config", "Configuration file" )
            }
        );
        JSAPResult parsedArgs = jsap.parse(args);

        RunScript.whiteNoiseDatasetAnalysis(
            parsedArgs.getString("raw"),
            parsedArgs.getString("root"),
            parsedArgs.getBoolean("genSubs"),
            new Config(parsedArgs.getString("movie_XML")),
            new Config(parsedArgs.getString("config")),
            parsedArgs.getBoolean("unCov"),
            parsedArgs.getBoolean("ei")
        );

        RunScript.printAnalysisTime(t1, System.currentTimeMillis());
//		System.exit(1);
    }
}