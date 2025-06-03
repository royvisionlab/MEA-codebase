package edu.ucsc.neurobiology.vision.convert;

import java.util.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.calculations.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class Strip512to61Data
    extends AbstractCalculation {

    private String rawFileName, saveToName;
    private int centerElectrode;


    public void startCalculation() throws Exception {
        Vision app = Vision.getInstance();

        int[] electrodes = FileConvert.get61ElectrodeSubset(centerElectrode);
        FileConvert.extractElectrodes(rawFileName, saveToName, electrodes, 0);

        app.sendMessage("Strip512to61Data: Done.");
        app.getCalculationManager().calculationDone();
    }


    public void setParameters(HashMap parameters) {
        rawFileName = (String) parameters.get("Raw_Data_File");
        saveToName = (String) parameters.get("Save_To_File");
        centerElectrode = Integer.parseInt( (String) parameters.get("electrode"));
    }
}
