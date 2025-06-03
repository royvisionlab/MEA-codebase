package edu.ucsc.neurobiology.vision.tasks;

import java.io.*;

import java.awt.event.*;
import javax.swing.*;

import com.martiansoftware.jsap.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.parameters.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * A class that can be run from the command line to perform the mapping of a given dataset
 * on a master dataset.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class MappingAnalysis extends AbstractAction {

    public MappingAnalysis() {
        super("Mapping Analysis");
    }

    public void actionPerformed(ActionEvent event) {
        Vision app = Vision.getInstance();
        
        ParametersTable table = app.getConfig().showDialog("Mapping Analysis", "Mapping Analysis", app.getMainFrame());
        if (table == null) return;

        final String masterFolder =
            ( (FileParameter) table.getParameter("Master Dataset Folder")).getValue();
        final String datasetFolder =
            ( (FileParameter) table.getParameter("Dataset Folder")).getValue();
        final String rawDataFile =
            ( (FileParameter) table.getParameter("Raw Data File")).getValue();
        final String defaults = ( (FileParameter) table.getParameter(
            "Defaults File")).getValue();

        Thread runner = new Thread() {
            public void run() {
                try {
                    RunScript.mappingAnalysis(masterFolder, datasetFolder, rawDataFile, new Config(defaults), false);
                } catch (Exception ex) { Vision.reportException(ex); }
            }
        };
        runner.start();
    }


    public static void main(String[] args) throws Exception {
        long t1 = System.currentTimeMillis();
        
        VisionJSAP jsap = new VisionJSAP( 
                MappingAnalysis.class.getName(), 
                new com.martiansoftware.jsap.Parameter[] {
                        new UnflaggedOption("master", JSAP.STRING_PARSER, JSAP.REQUIRED, "Master dataset name."),
                        new UnflaggedOption("root", JSAP.STRING_PARSER, JSAP.REQUIRED, "Root analysis folder."),
                        new UnflaggedOption("raw", JSAP.STRING_PARSER, JSAP.REQUIRED, "Raw data file."),
                        new Switch("genSubs", 'g', "gsub", "Generate subfolders?"),	
                        new FlaggedOption( "config", JSAP.STRING_PARSER, "config.xml", JSAP.NOT_REQUIRED, 'c', "config", 
                        "Configuration file" )
                }
        );
        JSAPResult parsedArgs = jsap.parse(args);
        
        String rootFolder = parsedArgs.getString("root");
        File f = new File(parsedArgs.getString("raw"));
        if (parsedArgs.getBoolean("genSubs"))
            rootFolder += File.separator + f.getParentFile().getName() + File.separator + f.getName();
        new File(rootFolder).mkdirs();

        Config config = new Config(parsedArgs.getString("config"));
        boolean projectionsOnly = false;		
        RunScript.mappingAnalysis(parsedArgs.getString("master"), rootFolder, parsedArgs.getString("raw"), config, projectionsOnly);

        RunScript.printAnalysisTime(t1, System.currentTimeMillis());
        System.exit(1);
    }

}