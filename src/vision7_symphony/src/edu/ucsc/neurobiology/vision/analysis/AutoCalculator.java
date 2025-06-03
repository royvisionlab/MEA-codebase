package edu.ucsc.neurobiology.vision.analysis;

import java.io.*;

import edu.ucsc.neurobiology.vision.plot.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class AutoCalculator
    implements ParametersCalculator {

    public double timeRange = Double.NaN;
    public double binning = Double.NaN;


    public void init(String rootPath, String mainFilePath) {

    }


    public String getName() {
        return "Auto";
    }


    public String[][] getParameterTypes() {
        return new String[][] { {"Auto", "DoubleArray"}
            , {"acfBinning", "Double"}
            , {"nSpikes", "Double"}
            , {"acfMean", "Double"}
            , {"acfRMS", "Double"}
        };
    }


    public Object[] getParameterValues(ParameterCalculatorContext c) throws
        IOException {

        // get the autocorrelation function
        int spikeCount = c.neuronFile.getSpikeCount(c.id);
        DoubleHistogram auto = AutocorrelationCalculator.calculate(
            c.neuronFile.getSpikeTimes(c.id), timeRange, binning);
        double[] autoArray = auto.toArray();
        double acfMean = 0.0;
        double acfRMS = 0.0;
        double norm = 0.0;
        double temp = 0.0;
        for(int i=0; i<autoArray.length; i++) {
            acfMean += binning*i*autoArray[i];
            acfRMS += binning*i*binning*i*autoArray[i];
            norm += autoArray[i];
            
            

        }
        acfMean/=norm;
        
        for(int i=0; i<autoArray.length; i++) {
            temp = binning*i-acfMean;
            acfRMS += temp*temp*autoArray[i];
            norm += autoArray[i];
            
        }
             
        acfRMS = Math.sqrt(acfRMS/norm);

        return new Object[] {
            auto.toArray(),
            new Double(binning),
            new Double(spikeCount),
            new Double(acfMean),
            new Double(acfRMS)
        };
    }
}
