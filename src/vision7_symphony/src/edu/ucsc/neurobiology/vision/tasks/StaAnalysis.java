package edu.ucsc.neurobiology.vision.tasks;

import java.awt.event.*;
import javax.swing.*;
import com.martiansoftware.jsap.*;
import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.parameters.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * This class performs an STA calculation only. It performs a subset of the tasks
 * that are completed in WhiteNoiseDatasetAnalysis
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 * @owner Tim Machado, The Salk Institute
 */
public class StaAnalysis
extends AbstractAction {


    public StaAnalysis() {
        super("STA Analysis");
    }


    public void actionPerformed(ActionEvent event) {
        Vision app = Vision.getInstance();

        ParametersTable table = Vision.getConfig().showDialog(
                "StaAnalysis", "STA Analysis", app.getMainFrame());
        if (table == null) {
            return;
        }

        final String outputFolder = table.getFileParameter("Analysis Folder");
        final String movieXMLFile = table.getFileParameter("Movie XML File");
        final String configFile = table.getFileParameter("Config File");

        Thread runner = new Thread() {
            public void run() {
                try {
                    RunScript.staAnalysis(outputFolder, new Config(movieXMLFile), new Config(configFile));
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
                StaAnalysis.class.getName(), 
                new com.martiansoftware.jsap.Parameter[] {
                        new UnflaggedOption("root", JSAP.STRING_PARSER, JSAP.REQUIRED, "Root analysis folder."),
                        new UnflaggedOption("movie_XML", JSAP.STRING_PARSER, JSAP.REQUIRED, "Movie xml file."),
                        new FlaggedOption("config",            JSAP.STRING_PARSER, "config.xml", JSAP.NOT_REQUIRED, 'c', "config",            "Configuration file" )
                }
        );
        JSAPResult parsedArgs = jsap.parse(args);

        RunScript.staAnalysis(
                parsedArgs.getString("root"),
                new Config(parsedArgs.getString("movie_XML")),
                new Config(parsedArgs.getString("config"))
        );
        RunScript.printAnalysisTime(t1, System.currentTimeMillis());
    }
}
