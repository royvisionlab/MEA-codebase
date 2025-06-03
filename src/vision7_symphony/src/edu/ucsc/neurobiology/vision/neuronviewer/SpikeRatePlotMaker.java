package edu.ucsc.neurobiology.vision.neuronviewer;

import java.io.*;

import java.awt.*;
import javax.swing.tree.*;

import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class SpikeRatePlotMaker
    extends PlotMaker {

    static int maxT = -1;
    public double binning = 10;
    private static double samplingFrequency = -1;
    private HistogramStyle spikeRateStyle = new HistogramStyle("Spike Rate");

    public SpikeRatePlotMaker() {
        super("Spike Rate", CLASS_AND_NEURON_PLOT);
    }


    public void initialize(NeuronViewer viewer) {
        super.initialize(viewer);

        int[] id = neuronFile.getIDList();
        maxT = neuronFile.getNumberOfSamples();
        samplingFrequency = neuronFile.getSamplingFrequency();

    }


    public Component makePlot(IntegerList list, int plotType, TreePath classPath) {
        try {
            return makeSpikeRatePanel(list.toArray());
        } catch (IOException ex) {
            return null;
        }
    }


    private PlotPanel makeSpikeRatePanel(int[] neuronsInClass) throws IOException {

        PlotPanel p = new PlotPanel("spike rate");
        DoubleHistogram spikeHist = getSpikeRateHistogram(
            neuronFile, neuronsInClass, 0, maxT / samplingFrequency,
            binning);

        double[] values = spikeHist.toArray();
        double[] x = new double[values.length];
        for (int i = 0; i < values.length; i++) {
            x[i] = i * binning + binning / 2;
        }

        ScatterPlot scatter = new ScatterPlot(x, values, new double[values.length]);
        p.addData(scatter, new ScatterPlotStyle(
            SymbolType.DISK, 0, Color.black, true, Color.black, 1f));
        p.autoscale();
        double[] r = p.getRange();
        r[2] = 0;
        p.setRange(r);
        p.setLabels("Time (s)", "Spike Rate (Hz)");
        return p;

    }


    /**
     *
     * @param nf NeuronFile
     * @param id int[]
     * @param t1 double in seconds
     * @param timeBinning double in seconds
     * @return DoubleHistogram
     * @throws IOException
     */
    public static DoubleHistogram getSpikeRateHistogram(
        NeuronFile nf, int[] id, double t1, double t2, double timeBinning) throws
        IOException {
        DoubleHistogram spikeRate = new DoubleHistogram("", t1, t2, timeBinning);

        for (int i = 0; i < id.length; i++) {
            int[] tList = nf.getSpikeTimes(id[i]);

            for (int j = 0; j < tList.length; j++) {
                double t = tList[j] / samplingFrequency;
                if (t >= t1) {
                    spikeRate.fill(t, 1);
                }
            }
        }

        spikeRate.scale(1.0 / (timeBinning * id.length));

        return spikeRate;
    }


    public static DoubleHistogram[] getSpikeRateHistogram(SpikeFile spikes,
        double spikeRateBinning) throws IOException {

        double samplingFrequency = 20000.0;
        final int nElectrodes = spikes.getNumberOfElectrodes();
        DoubleHistogram[] spikeRate = new DoubleHistogram[nElectrodes];
//        DoubleHistogram[] aplitudeHist = new DoubleHistogram[nElectrodes];
//        DoubleHistogram[] widthHist = new DoubleHistogram[nElectrodes];
        int lastTTL = spikes.getLastTTL();
        for (int i = 0; i < nElectrodes; i++) {
            spikeRate[i] = new DoubleHistogram("Sr(t)", 0, lastTTL / samplingFrequency,
                                               spikeRateBinning);
//            aplitudeHist[i] = new DoubleHistogram("Amp", 0, 2048, 1);
//            widthHist[i] = new DoubleHistogram("Width", 0, 20, 1);
        }

        // go through the data frame by frame and calculate the spike rate Sr(t)
        SpikeIterator spikeIterator = spikes.iterator();
        while (spikeIterator.hasNext()) {
            Spike spike = (Spike) spikeIterator.next();
            int electrode = spike.electrode;

            spikeRate[electrode].fill(spike.time / samplingFrequency, 1);
//            aplitudeHist[electrode].fill(spike.getAmplitude(), 1);
//            widthHist[electrode].fill(spike.getWidth(), 1);
        }

        return spikeRate;

    }


}
