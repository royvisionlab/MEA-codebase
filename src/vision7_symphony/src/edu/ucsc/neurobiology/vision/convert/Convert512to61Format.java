package edu.ucsc.neurobiology.vision.convert;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.calculations.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class Convert512to61Format
    extends AbstractCalculation {

    private String rawDataFileName, outputFolder;
    private boolean includeHeader;


    public void startCalculation() throws IOException {
        String datasetName =
            StringUtil.removeExtension(new File(rawDataFileName).getName());
        String newName = outputFolder + File.separator +
                         datasetName + "-16bit" + VisionParams.BIN_FILE_EXTENSION_512;

        System.out.println(rawDataFileName);
        System.out.println(newName);

        if (newName.equals(rawDataFileName)) {
            System.out.println(
                "Error: the input and output file names cannot be the same");
            Vision.getInstance().getCalculationManager().calculationDone();
            return;
        }

        File f = new File(newName);
        if (f.exists() && f.isFile()) {
            System.out.println(
                "Error: the output file exists. The file will not be overwritten. The program exits...");
            Vision.getInstance().getCalculationManager().calculationDone();
            return;
        }

        FileConvert.newFormat2oldFormat(rawDataFileName, newName, includeHeader);

        Vision.getInstance().getCalculationManager().calculationDone();
    }


    public void setParameters(HashMap<String, String> parameters) {
        rawDataFileName = ( (String) parameters.get("rawDataFileName")).trim();
        outputFolder = ( (String) parameters.get("outputFolder")).trim();
      includeHeader =  new Boolean((String) parameters.get("Include Header")).booleanValue();

    }

}
