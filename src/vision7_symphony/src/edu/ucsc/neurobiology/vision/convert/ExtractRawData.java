package edu.ucsc.neurobiology.vision.convert;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.calculations.*;
import edu.ucsc.neurobiology.vision.electrodemap.*;
import edu.ucsc.neurobiology.vision.io.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ExtractRawData extends AbstractCalculation {

    private String rawFileName, saveToName;
    private RawDataHeader header;

    private RawDataFile rawDataFile;
    private ElectrodeMap map;
    int electrodesCount;
    double time; // seconds
    double timeToSkip = 60; // seconds


    public void startCalculation() throws Exception {
        Vision app = Vision.getInstance();

        rawDataFile = new RawDataFile(new File(rawFileName));
        header = rawDataFile.getHeader();
        map = ElectrodeMapFactory.getElectrodeMap(header.getArrayID());
        int nSamples = (int) (time * 20000);

        ArrayList<Integer> usedElectrodes = new ArrayList<Integer>();
        short[] rawData = new short[nSamples];
        for (int n = 0; n < electrodesCount; n++) {
            int electrode;
            while (true) {
                electrode = (int) (Math.random() * map.getNumberOfElectrodes() - 1);
                if (!map.isDisconnected(electrode) &&
                    !usedElectrodes.contains(new Integer(electrode))) {
                    break;
                }
            }
            usedElectrodes.add(new Integer(electrode));

            app.sendMessage("ExtractRawData: Getting Raw Data for " + electrode);
            rawDataFile.getData(electrode, (int) (timeToSkip * 20000), rawData);

            app.sendMessage("ExtractRawData: Saving Raw Data for " + electrode);
            String name = "electrode" + electrode + ".txt";
            PrintWriter writer = new PrintWriter(new FileWriter(name));
            for (int i = 0; i < rawData.length; i++) {
                writer.println(rawData[i]);
            }
            writer.close();
        }

        app.sendMessage("ExtractRawData: Done.");
        app.getCalculationManager().calculationDone();
    }


    public void setParameters(HashMap<String, String> parameters) {
        rawFileName = (String) parameters.get("Raw_Data_File");
        saveToName = (String) parameters.get("Save_To_File");
        time = Double.parseDouble( (String) parameters.get("time"));
        electrodesCount = Integer.parseInt( (String) parameters.get("electrodesCount"));
    }
}
