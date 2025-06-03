package edu.ucsc.neurobiology.vision.snf;

import java.util.HashMap;

import edu.ucsc.neurobiology.vision.calculations.AbstractCalculation;

public class SerialNeuronFinding extends AbstractCalculation {

    String workFolder, binFile;
    public void setParameters(HashMap<String, String> parameters)
    throws Exception {

        workFolder = (String) parameters.get("WorkFolder");
        binFile = (String) parameters.get("BINFile");

    }


    public void startCalculation() throws Exception {
        SerialNeuronPanel snp = new SerialNeuronPanel(workFolder, binFile);
        
    }

}
