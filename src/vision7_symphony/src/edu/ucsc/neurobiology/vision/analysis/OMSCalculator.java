package edu.ucsc.neurobiology.vision.analysis;

import java.io.*;


/**
 *
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class OMSCalculator
    implements ParametersCalculator {

    public String datasetName;

    private OMSClassification omsClassification;


    public void init(String rootPath, String mainFilePath) throws Exception {
        File refFile = new File(mainFilePath);
        String refDataSetName = refFile.getName();
        omsClassification = new OMSClassification(datasetName, rootPath, refDataSetName);
    }


    public String getName() {
        return "OMS";
    }


    public String[][] getParameterTypes() {
        return new String[][] { {"omsDistances", "DoubleArray"}
            , {"omsSpikes", "DoubleArray"}
            , {"omsHistData", "DoubleArray"}
            , {"omsErrorHistData", "DoubleArray"}
            , {"omsMaxDistance", "Double"}
        };
    }


    public Object[] getParameterValues(ParameterCalculatorContext c) throws
        IOException {

        double x0 = ( (Double) c.getParameter("x0")).doubleValue();
        double y0 = ( (Double) c.getParameter("y0")).doubleValue();
        double[][] omsParams = omsClassification.calculateOMSParams(c.id, x0, y0);

        return new Object[] {
            omsParams[0],
            omsParams[1],
            omsParams[2],
            omsParams[3],
            new Double(omsParams[4][0])
        };
    }
}
