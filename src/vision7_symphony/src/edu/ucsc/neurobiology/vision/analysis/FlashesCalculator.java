package edu.ucsc.neurobiology.vision.analysis;

import java.io.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class FlashesCalculator
    implements ParametersCalculator {

    public String datasetName;
    public int nFlashBlocks;


    private FlashesClassification flashesClassification;


    public void init(String rootPath, String mainFilePath) throws Exception {
        flashesClassification = new FlashesClassification(
            new String[] {datasetName}, rootPath, nFlashBlocks, -1);
    }


    public String getName() {
        return "Flashes";
    }


    public String[][] getParameterTypes() {
        return new String[][] { {"flashResponse", "DoubleArray"}
            , {"flashBinSize", "Double"}
        };
    }


    public Object[] getParameterValues(ParameterCalculatorContext c) throws
        IOException {

        return new Object[] {
            flashesClassification.getFlashResponse(c.id),
            new Double(flashesClassification.getFlashBinning())
        };
    }
}
