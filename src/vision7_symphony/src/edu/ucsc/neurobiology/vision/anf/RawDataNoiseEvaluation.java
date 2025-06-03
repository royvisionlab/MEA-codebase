package edu.ucsc.neurobiology.vision.anf;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.calculations.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * Evaluates the noise level on every electrode by using the second 5 seconds of the raw
 * data. The spikes are eliminated by an iterative thresholding procedure.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class RawDataNoiseEvaluation
    extends AbstractCalculation {

    public final static int standardSigmaNoise = 20;

    private String rawDataFileName, outputFolder;
    private RawDataHeader header;
    private int nElectrodes;
    private Vision app;
    private RawDataFile rawDataFile;
    short[][] rawData;
    double time = 5; // seconds
    double timeToSkip = 5; // seconds
    int nSamples;


    public void startCalculation() throws IOException {
        app = Vision.getInstance();

        app.sendMessage("RawDataNoiseEvaluation: Getting Raw Data...");
        DataFileStringParser parser = new DataFileStringParser(rawDataFileName);
        rawDataFile = new RawDataFile(new File(parser.getDatasets()[0]));
        double startTime = parser.getStartTimes()[0];
        
        header = rawDataFile.getHeader();
        this.nElectrodes = header.getNumberOfElectrodes();
        nSamples = (int) (time * 20000);
        rawData = new short[nSamples][nElectrodes];
        rawDataFile.getData( (int) ((timeToSkip+startTime) * 20000), rawData);

        // start
        app.sendMessage("RawDataNoiseEvaluation: Calculating rms noise...");

        double[] rmsSigma = new double[nElectrodes];

        app.startProgressBar();
        for (int electrode = 1; electrode < nElectrodes; electrode++) {
            short[] data = new short[nSamples];
            for (int i = 0; i < nSamples; i++) {
                data[i] = rawData[i][electrode];
            }
            rmsSigma[electrode] = getNoise(data, nSamples);

            app.setProgress(100 * electrode / (nElectrodes - 1));
        }
        app.endProgressBar();

        app.sendMessage("Raw Data Noise Evaluation...");
        String setName = new File(outputFolder).getName();
        String sigmasFileName = outputFolder + File.separator + setName +
                                VisionParams.NOISE_FILE_EXTENSION;

        File f = new File(sigmasFileName);
        if (f.exists()) {
            f.delete();
        }
        PrintWriter pw = new PrintWriter(new FileWriter(sigmasFileName));
        for (int i = 0; i < nElectrodes; i++) {
            pw.println(rmsSigma[i]);
        }
        pw.close();

        Vision.getInstance().getCalculationManager().calculationDone();
    }


    public static double getNoise(short[] data, int nSamples) {
        MeanVarianceCalculator mvc =
            new MeanVarianceCalculator(MeanVarianceCalculator.UNBIASED);
        final int samplesToSkip = 0;
        double sigma = Double.POSITIVE_INFINITY;
        double mean = 0;

        while (true) {
            for (int i = samplesToSkip; i < nSamples; i++) {
                if (Math.abs(data[i] - mean) < 3 * sigma) {
                    mvc.add(data[i]);
                }
            }

            double sigmaNew = mvc.getStandardDeviation();
            mean = mvc.getMean();
//            System.out.println("==> " + sigmaNew);
            if (Math.abs(sigma - sigmaNew) <= 0.01) {
                sigma = sigmaNew;
                break;
            } else {
                sigma = sigmaNew;
            }
        }

        return sigma;
    }


    public void setParameters(HashMap parameters) {
        rawDataFileName = (String) parameters.get("Raw_Data_File");
        outputFolder = (String) parameters.get("Save_To_File");
    }
}
