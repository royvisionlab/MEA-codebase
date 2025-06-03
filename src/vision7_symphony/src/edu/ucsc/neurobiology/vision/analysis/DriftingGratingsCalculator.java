package edu.ucsc.neurobiology.vision.analysis;

import java.io.*;


/**
 *
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class DriftingGratingsCalculator
    implements ParametersCalculator {

    public String datasetName;

    private DriftingGratings driftingGratings;

    public void init(String rootPath, String mainFilePath) throws Exception {
        String[] sinusoidsFileName;
        if (datasetName != null && datasetName.trim().length() != 0) {
            sinusoidsFileName = datasetName.split("#", 0);
        } else {
            sinusoidsFileName = null; /*new String[] {};*/
        }

        driftingGratings = new DriftingGratings(sinusoidsFileName, rootPath);
    }


    public String getName() {
        return "Drifting Gratings";
    }


    public String[][] getParameterTypes() {
        return new String[][] { {"spatialFrequencies", "DoubleArray"}
            , {"temporalFrequencies", "DoubleArray"}
            , {"directions", "DoubleArray"}
            , {"gratingResponse", "DoubleArray"}
            , {"averageGrating", "DoubleArray"}

            , {"xOffDS", "Double"}
            , {"yOffDS", "Double"}
            , {"magOffDS", "Double"}
            , {"angOffDS", "Double"}
            , {"xOnDS", "Double"}
            , {"yOnDS", "Double"}
            , {"magOnDS", "Double"}
            , {"angOnDS", "Double"}
            , {"xOS", "Double"}
            , {"yOS", "Double"}
            , {"magOS", "Double"}
            , {"angOS", "Double"}
            , {"spatialMean", "Double"}
            , {"temporalMean", "Double"}
            , {"spatialRMS", "Double"}
            , {"temporalRMS", "Double"}

        };
    }


    public Object[] getParameterValues(ParameterCalculatorContext c) throws
        IOException {
        double[][] directionRates = driftingGratings.calculate(c.id);

        //Calculate Parameters for Off DS cells
        double[] offDSParams = driftingGratings.calculateDSParams(
            directionRates[0], 32.0, 256.0, 16.0, 256.0);

        //Calculate Parameters for On DS cells
        double[] onDSParams = driftingGratings.calculateDSParams(
            directionRates[0], 64.0, 64.0, 64.0, 2048.0);

        //Calculate Parameters for OS cells
        double[] osParams = driftingGratings.calculateDSParams(
            directionRates[0], 0.0, 2048.0, 0.0, 2048.0);
        
        double[] meansAndRms = driftingGratings.calculateMeansAndRms(directionRates[1]);

        return new Object[] {

            driftingGratings.directionChoices[0],
            driftingGratings.directionChoices[1],
            driftingGratings.directionChoices[2],
            directionRates[0],
            directionRates[1],

            new Double(offDSParams[0]),
            new Double(offDSParams[1]),
            new Double(offDSParams[2]),
            new Double(offDSParams[3]),

            new Double(onDSParams[0]),
            new Double(onDSParams[1]),
            new Double(onDSParams[2]),
            new Double(onDSParams[3]),

            new Double(osParams[4]),
            new Double(osParams[5]),
            new Double(osParams[6]),
            new Double(osParams[7]),
            
            new Double(meansAndRms[0]),
            new Double(meansAndRms[1]),
            new Double(meansAndRms[2]),
            new Double(meansAndRms[3])
        };
    }
}
