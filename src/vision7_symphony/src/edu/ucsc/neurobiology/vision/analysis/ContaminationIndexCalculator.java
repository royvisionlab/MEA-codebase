package edu.ucsc.neurobiology.vision.analysis;

import java.io.*;


/**
 *
 * @author Alexander Sher, University of California, Santa Cruz
 */
public class ContaminationIndexCalculator
    implements ParametersCalculator {

    public void init(String rootPath, String mainFilePath) {

    }


    public String getName() {
        return "Contamination Index";
    }


    public String[][] getParameterTypes() {
        return new String[][] { {"contamination", "Double"}
        };
    }


    public Object[] getParameterValues(ParameterCalculatorContext c) {
        // get the contamination index
        double contamination = Double.NaN;

        try {
            contamination = AutocorrelationCalculator.getContamination(c.id, c.neuronFile);
        } catch (IOException e) {
            e.printStackTrace();
        }

        return new Object[] {
            new Double(contamination)
        };
    }
}
