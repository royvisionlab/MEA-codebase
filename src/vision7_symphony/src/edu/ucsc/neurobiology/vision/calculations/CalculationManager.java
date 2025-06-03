package edu.ucsc.neurobiology.vision.calculations;

import java.io.IOException;
import java.io.InputStream;
import java.io.FileInputStream;
import java.util.*;

import javax.xml.parsers.ParserConfigurationException;
import org.xml.sax.SAXException;

import com.martiansoftware.jsap.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.parameters.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * This manages a set of AbstractCalculations.  It allows them to be called programatically.
 *
 * @author Charles A. Loomis, University of California, Santa Cruz
 * @author Dumitru Petrusca, University of California, Santa Cruz
 * @author Peter H. Li, The Salk Institute
 * @owner Matthew Grivich, The Salk Institute
 */
public class CalculationManager {
    private final static String rootPackage = "edu.ucsc.neurobiology.vision";
    private final static String baseConfigPath = "/edu/ucsc/neurobiology/vision/baseconfig.xml";
    private final static String paramsDelimiter = "::";
    
    LinkedHashMap<String, String> calculations = new LinkedHashMap<String, String>();
    private String calculationName;
    private long startTime;
    private SynchronizationObject syncObject = new SynchronizationObject();


    public CalculationManager() {
        // add the calculations
        calculations.put("Raw Data Filtering", "anf.RawDataFiltering");
        calculations.put("Raw Data Noise Evaluation", "anf.RawDataNoiseEvaluation");
        calculations.put("Spike Finding", "anf.SpikeFinding");
        
        calculations.put("1", null);
        calculations.put("Noise Whitened Covariances", "anf.NWCovariance");
        calculations.put("PCA Neuron Finding: Projections", "anf.PCANFProjections");
        calculations.put("PCA Neuron Finding: Clustering", "anf.PCANFClustering");
        calculations.put("PCA Neuron Finding: Mapping", "anf.PCANeuronMapping");
        calculations.put("Neuron Cleaning", "anf.NeuronCleaning");
        calculations.put("Join", "anf.Join");
        
        calculations.put("2", null);
        calculations.put("Make Parameters File", "analysis.MakeParametersFile");
        calculations.put("Electrophysiological Imaging", "analysis.PhysiologicalImaging");
        calculations.put("Electrophysiological Imaging Fast", "analysis.Imaging");
        calculations.put("Electrophysiological Imaging Optimized", "analysis.ImagingNew");
        calculations.put("Calculate Auxiliary Parameters", "analysis.AuxiliaryParametersCalculator");
        calculations.put("STA Calculation", "analysis.STACalculation");
        calculations.put("STA Calculation Parallel", "analysis.STACalculationParallel");
        calculations.put("STA Join", "analysis.STAJoin");
        calculations.put("Whiten Natural Power STAs",
                         "analysis.WhitenNaturalSTAsCalculation");
        
        calculations.put("Make White Noise Movie", "stimulus.CreateWhiteNoiseMovie");
        calculations.put("Make OMS Movie", "stimulus.CreateOMSMovie");
        calculations.put("Make Power Law Movie", "stimulus.CreatePowerLawMovie");
 
        calculations.put("4", null);
        calculations.put("SNF Preprocessing", "snf.SNFPreprocessing");
        calculations.put("Serial Neuron Finding", "snf.SerialNeuronFinding");
        calculations.put("Serial Neuron Mapping", "snf.SerialNeuronMapping");

        calculations.put("5", null);
        calculations.put("Extract 61 Dataset from 512 File",
                         "convert.Strip512to61Data");
        calculations.put("Convert Raw Data from OLD to NEW Format",
                         "convert.Convert61to512Format");
        calculations.put("Convert Raw Data from NEW to OLD Format",
                         "convert.Convert512to61Format");
        calculations.put("Convert Raw Data To Electrode Major Format",
                         "convert.ConvertRawDataToElectrodeMajor");
        calculations.put("Extract Raw Data", "convert.ExtractRawData");
        calculations.put("Compress File or Folder Using bzip2", "convert.BZip2Compress");
        calculations.put("Decompress File or Folder Using bzip2", "convert.BZip2Decompress");
        calculations.put("Generate Globals Files", "convert.AddGlobalsFiles");
        calculations.put("Copy Raw Data Header to Globals", "convert.CopyRawDataHeader512ToGlobals");
        
        calculations.put("6", null);
        calculations.put("Generation Currents", "analysis.GenerationCurrents");
    }

    public List<String> getCalculationNames() {
        Set<String> keys = calculations.keySet();
        
        // Remove the number keys, whose values are null; these are placeholders for display purposes
        List<String> keptKeys = new ArrayList<String>();
        for (String key : keys) {
            if (calculations.get(key) != null) {
                keptKeys.add(key);
            }
        }
        
        return keptKeys;
    }

    public String getCurrentCalculationName() {
        return calculationName;
    }


    public boolean isCurrentlyCalculating() {
        return (syncObject.getState() == syncObject.WORKING);
    }


    /*synchronized */ public void runCalculation(final String cName, HashMap<String, String> parameters) {
        runCalculation(cName, parameters, true);
    }


    /*synchronized */ void runCalculation(final String cName, HashMap<String, String> parameters, boolean block) {
        System.out.println("\nRunning: " + cName + ", at " + new Date());
        calculationName = cName;
        startTime = System.currentTimeMillis();
        final AbstractCalculation calculation;

        try {
            Class cClass = Class.forName(rootPackage + "." + calculations.get(cName));
            calculation = (AbstractCalculation) cClass.newInstance();
            calculation.setParameters(parameters);
        } catch (Exception e) {
            if (Vision.isConsoleBased()) {
                Vision.reportFatalException(
                    cName + " could not run because this exception occured:", e);
            } else { // there is a GUI, do nothing
                Vision.reportException(
                    cName + " could not run because this exception occured:", e);
            }

            syncObject.done();
            return;
        }

        // show the diagnostic panel
        if (Vision.isGUIBased()) {
            Vision.getInstance().getCalculationManagerGUI().setDiagnosticsPanel(
                calculation, calculationName);
        }

        syncObject.setWorking();

        Thread calculationThread = new Thread() {
            public void run() {
                try {
                    calculation.startCalculation();
                } catch (Exception e) {
                    if (Vision.isConsoleBased()) {
                        Vision.reportFatalException("Cannot run Calculation: " + cName, e);
                    } else { // there is a GUI, do nothing
                        Vision.reportException(
                            cName + " could not run because this exception occured:", e);
                        syncObject.done();
                    }
                }
            }
        };
        calculationThread.start();

        // wait until the calculation finishes
        if (block) {
            syncObject.waitUntilDone();
        }
    }


    public void calculationDone() {
        double t = (System.currentTimeMillis() - startTime) / 1000.0;
        String text = calculationName + " done at " + new Date() +
                      ". Took " + StringUtil.format(t, 1) + "sec.\n";
        System.out.println(text);

        if (Vision.isGUIBased()) {
            Vision.getInstance().sendMessage(text);
            Vision.getInstance().setProgress(0);
            Vision.getInstance().getCalculationManagerGUI().removeDiagnosticsPanels();
        }

        syncObject.done();
    }


    public static void runCalculation(String cName, String[] args) {

        Config config = Vision.getInstance().getConfig();
        ParametersTable t = config.getParameterGroup(cName);
        HashMap<String, String> p = new HashMap<String, String>();

        if (args.length == t.getParametersCount()) {
            for (int i = 0; i < args.length; i++) {
                p.put(t.getParameter(i).getName(), args[i]);
            }
        } else {
            System.out.println(
                "Incorrect number of command line arguments: " + t.getParametersCount() +
                " required:");
            for (int i = 0; i < t.getParametersCount(); i++) {
                System.out.println(t.getParameter(i).getName());
            }
            System.out.println("Choose parameters from the GUI...");
            t = config.showDialog(cName, cName, null);
            if (t == null) {
                System.out.println(
                    "You did not provide the required input. The program will now exit.");
                System.exit(1);
            } else {
                for (int i = 0; i < t.getParametersCount(); i++) {
                    p.put(t.getParameter(i).getName(), t.getParameter(i).valueAsString());
                }
            }
        }

        Vision.getInstance().getCalculationManager().runCalculation(cName, p, true);

        System.exit(0);
    }


    private static LinkedHashMap<String, String> getCalcParams(String cName, Config config) {
        LinkedHashMap<String, String> params = null;

        // Surprising that this throws an exception, as Config#getParameterList() is written to return null
        // if cName isn't found...
        try {
            params = config.getParameterList(cName);
        } catch (Exception e) {
            System.err.println("\nCalculation '" + cName + "' not found in config-xml file: " + config.fileName);
            System.err.print(validCalcsString());
            System.exit(1);
        }
        return params;
    }

    private static String validCalcsString() {
        String vcs = "\nAvailable Calculations: \n";
        for (String c :	Vision.getInstance().getCalculationManager().getCalculationNames()) {
            vcs += " " + c + "\n";
        }
        return vcs;
    }
    
//    private static boolean isSubgroup(String groupName, String subgroupName) {
//    	return subgroupName.length() > groupName.length() && subgroupName.startsWith(groupName + ".");
//    }
    
//    private static Set<String> getSubgroups(HashMap<String, String> params, String groupName) {
//    	Set<String> subgroups = new HashSet<String>();
//    	for (String p : params.keySet()) {
//    		if (isSubgroup(groupName, p)) subgroups.add(p);
//    	}
//    	return subgroups;
//    }
    
    public static void main(String[] args) throws Exception {    	
        VisionJSAP jsap = new VisionJSAP( 
            CalculationManager.class.getName(), 
            new com.martiansoftware.jsap.Parameter[] {	
                new FlaggedOption(  "oldconfig",  JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, 'c', "config",          "Configuration file." ),
                new FlaggedOption(  "baseconfig", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, 'b', "base-config",     "Base Configuration file." ),
                new FlaggedOption(  "config",     JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, 'd', "newconfig",       "Defaults Configuration file." ),
                new Switch("listcalcs",  'l', "list-calcs",  "List available calculations."),
                new Switch("listparams", 'p', "list-params", "List valid parameters."),
                new Switch("version",    'V', "version",     "Print version."),
                new Switch("verbose",    'v', "verbose",     "Verbose output."),
                new UnflaggedOption("calc",       JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.REQUIRED,     JSAP.NOT_GREEDY, "Calculation"),
                new UnflaggedOption("params", 	  JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.GREEDY,     "Calculation parameters")
            }
        );
        JSAPResult parsedArgs = jsap.parseWithMangledNegs(args);
        
        if (parsedArgs.getBoolean("listcalcs")) {
            System.out.println(validCalcsString());
            System.exit(0);
        }
        if (parsedArgs.getBoolean("version")) {
            String versionString = "Vision v" + VisionParams.versionString();
            String buildDateString = Vision.buildDateString();
            if (buildDateString != null) versionString += ", build " + buildDateString;
            System.out.println(versionString);
            System.exit(0);
        }
        
        String group = parsedArgs.getString("calc");        
        LinkedHashMap<String, String> params = buildParams(parsedArgs);        
        System.out.println("Running: " + group);
        Vision.getInstance().getCalculationManager().runCalculation(group, params, true);
    }

    private static LinkedHashMap<String, String> buildParams(JSAPResult parsedArgs) throws Exception {
        String group = parsedArgs.getString("calc");
        String[] calcParams = VisionJSAP.demangleNegs(parsedArgs.getStringArray("params"));

        if (parsedArgs.contains("oldconfig")) return buildParamsOld(group, calcParams, parsedArgs);
        else						  		  return buildParamsNew(group, calcParams, parsedArgs);
    }
    
    private static LinkedHashMap<String, String> buildParamsNew(String calc, String[] calcParams, JSAPResult parsedArgs) throws ParserConfigurationException, SAXException, IOException {    	
        // Was the embedded source for list of valid parameters overridden?
        InputStream baseConfigStream;
        if (parsedArgs.contains("baseconfig")) {
            String bcp = parsedArgs.getString("baseconfig"); 
            baseConfigStream = new FileInputStream(bcp);
        } else {
            String bcp = baseConfigPath;
            baseConfigStream = Vision.class.getResourceAsStream(bcp);
        }

        // Get list of valid parameters
        Config baseConfig = new Config(baseConfigStream);
        baseConfigStream.close();
        LinkedHashMap<String, String> baseParams = getCalcParams(calc, baseConfig);
        
        if (parsedArgs.getBoolean("listparams")) {
            System.out.println(validParamsString(calc, baseParams));
            System.exit(0);
        }

        // The params and where they came from
        LinkedHashMap<String, String> params       = new LinkedHashMap<String, String>();
        LinkedHashMap<String, String> paramSources = new LinkedHashMap<String, String>(); 

        // Pull param values from base file; this allows the passed config to have only the values to be overridden
        for (String p : baseParams.keySet()) {
            params.put(p, baseParams.get(p));
            paramSources.put(p, "baseconfig-xml");
        }
        
        // Load params from given config.xml (e.g. movie.xml or primate.xml) to override base
        if (parsedArgs.contains("config")) {
            String configPath = parsedArgs.getString("config");
            if (parsedArgs.getBoolean("verbose")) System.out.println("\nConfig-xml: " + configPath);
            
            Config config = new Config(configPath);
            LinkedHashMap<String, String> configParams = getCalcParams(calc, config);

            for (String p : configParams.keySet()) {
                if (baseParams.containsKey(p)) {
                    params.put(p, configParams.get(p));
                    paramSources.put(p, "config-xml");
                } else {
                    System.err.println("\nUnrecognized calculation parameter '" + p + "' given in config-xml: " + configPath);
                    System.err.print("\nValid parameters:");
                    System.err.println(validParamsString(calc, baseParams));
                    System.exit(1);
                }
            }
        }
        
        // Load params from command-line arguments, process named ones, hold unnamed ones to process in next loop
        LinkedList<String> unnamedParams = new LinkedList<String>();
        for (int i = 0; i < calcParams.length; i++) {
            String[] parsedParam = calcParams[i].split(paramsDelimiter, 2);

            if (parsedParam.length < 2) {
                unnamedParams.add(calcParams[i]);
                continue;
            }
            
            if (!baseParams.containsKey(parsedParam[0])) {
                System.err.println("\nUnrecognized calculation parameter on command line: " + parsedParam[0]);
                System.err.println(validParamsString(calc, baseParams));
                System.exit(1);
            }

            params.put(parsedParam[0], parsedParam[1]);
            paramSources.put(parsedParam[0], "command-line-named");
        }
        
        // Process unnamed command-line args in order they appear in XML file
        Iterator<String> paramNames = baseParams.keySet().iterator();
        for (int i = 0; i < unnamedParams.size(); i++) {
            String pn = paramNames.next();
            params.put(pn, unnamedParams.get(i));
            paramSources.put(pn, "command-line-unnamed");
        }
        
        if (parsedArgs.getBoolean("verbose")) {
            System.out.println("Parameter settings:");
            for (String p : params.keySet()) {
                System.out.println("  " + p + ": " + params.get(p) + " (" + paramSources.get(p) + ")");
            }
        }
        
        return params;
    }
    
    private static String validParamsString(String calc, LinkedHashMap<String, String> baseParams, LinkedHashMap<String, String> exampleParams) {
        String validParamsStr = "\nValid parameters for " + calc + ":";
        for (String p : baseParams.keySet()) {
            validParamsStr += "\n  " + p;
            if (exampleParams.containsKey(p) && exampleParams.get(p).length() > 0) {
                validParamsStr += " (e.g.: " + exampleParams.get(p) + ")";
            }
        }
        return validParamsStr;
    }
    
    private static String validParamsString(String calc, LinkedHashMap<String, String> baseParams) {
        return validParamsString(calc, baseParams, baseParams);
    }
        
    public static LinkedHashMap<String, String> buildParamsOld(String group, String[] calcParams, JSAPResult parsedArgs) throws IOException, SAXException, ParserConfigurationException {
        Config config = new Config(parsedArgs.getString("oldconfig"));        
        LinkedHashMap<String, String> params = getCalcParams(group, config);
        if (calcParams.length == params.size()) {
            int i = 0;
            for (String name : params.keySet()) {
                params.put(name, calcParams[i]);
                System.out.println("Parameter '" + name + "' is set to: " + calcParams[i]);
                i++;
            }
        } else {
            String error = "\nIncorrect number of calculation parameters. Required:";
            error += validParamsString(group, params);
            System.err.println(error);
            System.exit(1);
        }
        
        return params;
    }

}