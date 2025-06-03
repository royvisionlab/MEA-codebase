package edu.ucsc.neurobiology.vision.snf;

import java.util.HashMap;

import edu.ucsc.neurobiology.vision.calculations.AbstractCalculation;

public class SerialNeuronMapping extends AbstractCalculation {
    String masterDatasetFolder, datasetFolder, originalRawFileName;
    int automationLevel;
    
    public void setParameters(HashMap<String, String> parameters) {
        masterDatasetFolder = (String) parameters.get("MasterFolder");
        datasetFolder = (String) parameters.get("WorkFolder");
        originalRawFileName = (String) parameters.get("BINFile");
        automationLevel = (int) Double.parseDouble( (String) parameters.get(
        "AutomationLevel"));


    }

    public void startCalculation() throws Exception {
        SerialMappingPanel snp = new SerialMappingPanel
        (masterDatasetFolder, datasetFolder,originalRawFileName, automationLevel);

    }

}
