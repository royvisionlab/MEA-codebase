package edu.ucsc.neurobiology.vision.analysis;

import java.io.*;

import java.awt.*;

import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.io.tags.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class FlashesClassification {
    private int nFlashBlocks;
    private double binning; // = 200; //bin size of histograms (samples)
    private double firstTime, lastTime, period;
    private NeuronFile[] neuronFiles;


    /**
     *
     * @param flashFileName String[]
     * @param rootPath String
     * @param nFlashBlocks int
     * @param binning double in seconds
     * @throws IOException
     */
    public FlashesClassification(String[] flashFileName, String rootPath,
                                 int nFlashBlocks, double binning) throws IOException {

        this.nFlashBlocks = nFlashBlocks;
        neuronFiles = new NeuronFile[flashFileName.length];
        for (int i = 0; i < flashFileName.length; i++) {

            if (flashFileName[i] != null && flashFileName[i].trim().length() != 0) {
                String flashFilePath =
                    rootPath + File.separator + flashFileName[i] + File.separator +
                    flashFileName[i] + VisionParams.NEURON_FILE_EXTENSION;
                neuronFiles[i] = new NeuronFile(flashFilePath);
//                System.out.println("Flashes from: " + flashFilePath);
            }
        }

        int[] ttl = neuronFiles[0].getTTLTimes();
        if (ttl == null) {
            throw new IllegalArgumentException("array ttls is null");
        }

        // calculate the average TTL separation
        // the TTLs are equally spaced in full screen flashes
        period = 0;
        int ttlsPerPeriod = Math.round(ttl.length / nFlashBlocks);
        double count = 0;
        for (int i = ttlsPerPeriod; i < ttl.length; i += ttlsPerPeriod) {
            period += ttl[i] - ttl[i - ttlsPerPeriod];
            count++;
        }
        period /= count;
        period = period * ttl.length / nFlashBlocks / ttlsPerPeriod;
        firstTime = ttl[0];
        lastTime = ttl[0] + nFlashBlocks * period;

        if (binning > 0) {
            this.binning = binning * 20000;
        } else {
            this.binning = Math.floor(period / 200.0);
//            System.err.println(this.binning);
        }
    }


    synchronized public double getFlashBinning() {
        return binning;
    }


    synchronized public double[] getFlashResponse(int id) throws IOException {
        return getFlashResponseHistogram(id).toArray();
    }


    synchronized public DoubleHistogram getFlashResponseHistogram(int id) throws
        IOException {
        // Obtain the response to flashes

        DoubleHistogram flasResponseHist =
            getFlashResponse(neuronFiles[0].getSpikeTimes(id), binning);

        for (int n = 0; n < neuronFiles.length; n++) {
            if (neuronFiles[n] != null) {
                DoubleHistogram flasResponseHist2 =
                    getFlashResponse(neuronFiles[n].
                                     getSpikeTimes(id), binning);
                for (int i = 0; i < flasResponseHist.getBinCount(); i++) {
                    flasResponseHist.fillBin(i, flasResponseHist2.getBin(i));
                }
            }
        }

        for (int i = 0; i < flasResponseHist.getBinCount(); i++) {
            flasResponseHist.scaleBin(i, 1.0 / neuronFiles.length);
        }

//        PlotUtilities.showData("", makeFlashPanel(flasResponseHist.toArray()));

        return flasResponseHist;
    }


    /**
     *
     * @param times int[]
     * @param binning double
     * @return DoubleHistogram
     */
    synchronized public DoubleHistogram getFlashResponse(
        int[] times, double binning) {

        if (times == null) {
            throw new IllegalArgumentException("array times is null");
        }

        DoubleHistogram periodHistogram = new DoubleHistogram(
            "", 0, period / 20000.0, binning / 20000.0);

        for (int i = 0; i < times.length; i++) {
            if (times[i] > firstTime && times[i] < lastTime) {
                periodHistogram.fill( ( (times[i] - firstTime) % period) / 20000.0, 1);
            }
        }

        // perform the average
        periodHistogram.scale(1.0 / nFlashBlocks);

        // convert to spike rate
        periodHistogram.scale(20000.0 / binning);

        return periodHistogram;
    }


    synchronized static public PlotPanel makeFlashPanel(ParametersFile paramsFile, int id) {
        return makeFlashPanel(paramsFile.getArrayCell(id, "flashResponse"),
                              paramsFile.getDoubleCell(id, "flashBinSize"));
    }


    synchronized public static PlotPanel makeFlashPanel(double[] flashResponse,
        double binning) {
        PlotPanel p = new PlotPanel("flashPanel");
        if (flashResponse != null) {
            ScatterPlot scatter = new ScatterPlot();
            for (int i = 0; i < flashResponse.length; i++) {
                scatter.add( (double) i * binning / 20000.0, flashResponse[i]);
            }
            ScatterPlotStyle style = new ScatterPlotStyle();
            style.setConnectingPoints(true);
            style.setSymbolType(SymbolType.NONE);
            p.addData(scatter, style);
        }
        p.setLabels("Time (s)", "Normalized Spike Rate");
        p.autoscale();

        return p;
    }


    synchronized public static PlotPanel makeFlashesPanel(ParametersFile paramsFile,
        int[] neuronsInClass) {
        PlotPanel flashesPanel = new PlotPanel("flashesPanel");
        ScatterPlotStyle style = new ScatterPlotStyle(
            SymbolType.NONE, 1, Color.black, true, Color.black, .25f);
        if (neuronsInClass.length != 0) {
            double binning = paramsFile.getDoubleCell(neuronsInClass[0], "flashBinSize");
            double[] flash;
            for (int i = 0; i < neuronsInClass.length; i++) {
                flash = ( (DoubleArrayTag) paramsFile.getCell(neuronsInClass[i],
                    "flashResponse")).getArray();

                if (flash != null) {
                    double sum = 0;
                    for (int j = 0; j < flash.length; j++) {
                        sum += flash[j];
                    }
                    for (int j = 0; j < flash.length; j++) {
                        flash[j] /= sum;

                    }

                    ScatterPlot scatter = new ScatterPlot();
                    for (int j = 0; j < flash.length; j++) {
                        scatter.add( (double) j * binning / 20000.0, flash[j]);
                    }

                    flashesPanel.addData(scatter, style);

                }
            }
            flashesPanel.setLabels("Time (s)", "Spike Rate (Hz)");
            flashesPanel.autoscale();
            return flashesPanel;
        } else {
            return new PlotPanel();
        }
    }


}
