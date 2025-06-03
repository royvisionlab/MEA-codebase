package edu.ucsc.neurobiology.vision.analysis;

import java.io.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ReversingGratingsCalculator
    implements ParametersCalculator {

    public String datasetName;
    public int binsPerPeriod;

    private ReversingGratings revGratings = null;


    public String getName() {
        return "Reversing Gratings";
    }


    public void init(String rootPath, String mainFilePath) throws Exception {
        if (datasetName != null && datasetName.trim().length() != 0) {
            revGratings = new ReversingGratings(
                rootPath + File.separator + datasetName, binsPerPeriod);
        }
    }


    public String[][] getParameterTypes() {
        return new String[][] { {"reversingFrequencies", "DoubleArray"}
            , {"T1reversingF1", "DoubleArray"}
            , {"T1reversingF2", "DoubleArray"}
            , {"T2reversingF1", "DoubleArray"}
            , {"T2reversingF2", "DoubleArray"}
            , {"T3reversingF1", "DoubleArray"}
            , {"T3reversingF2", "DoubleArray"}
            , {"nonLinIndex", "Double"}
        };
    }


    /*
        revGratings.setCurrentNeuron(id);
        double[][][] gratingsInfo = revGratings.getHarmonic(1);
        row[48] = gratingsInfo[0][0];
        row[49] = gratingsInfo[0][1];
        row[50] = gratingsInfo[0][2];
        if (gratingsInfo.length > 1) {
            row[51] = gratingsInfo[1][1];
            row[52] = gratingsInfo[1][2];
        }
        if (gratingsInfo.length > 2) {
            row[53] = gratingsInfo[2][1];
            row[54] = gratingsInfo[2][2];
        }
     */

    public Object[] getParameterValues(ParameterCalculatorContext c) {
        if (revGratings != null) {
            double[][][] gratingsInfo = null; //[temporalPeriods][spatialFreq:F1:F2][spatialPeriods]
            try {
                revGratings.setCurrentNeuron(c.id);
                gratingsInfo = revGratings.getHarmonic(1);
            } catch (IOException ex) {
                ex.printStackTrace();
                return null;
            }
            
            double maxNonlin = 0.0;
            if(gratingsInfo != null) {
                double val = 0;
                for(int temporal=0; temporal< gratingsInfo.length; temporal++) {
                    for(int spatial=0; spatial < gratingsInfo[0][0].length; spatial++) {
                        val = gratingsInfo[temporal][2][spatial]/gratingsInfo[temporal][1][spatial];
                        if(val > maxNonlin) {
                             maxNonlin = val;
                        }
                    }
                    
                }
            }

            return new Object[] {
                gratingsInfo[0][0],
                gratingsInfo[0][1],
                gratingsInfo[0][2],
                (gratingsInfo.length > 1) ? gratingsInfo[1][1] : null,
                (gratingsInfo.length > 1) ? gratingsInfo[1][2] : null,
                (gratingsInfo.length > 2) ? gratingsInfo[2][1] : null,
                (gratingsInfo.length > 2) ? gratingsInfo[2][2] : null,
                (gratingsInfo != null) ? maxNonlin : null
            };
        }

        return new Object[] {
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null
        };
    }

}
