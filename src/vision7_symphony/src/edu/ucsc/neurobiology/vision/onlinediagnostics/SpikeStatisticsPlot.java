package edu.ucsc.neurobiology.vision.onlinediagnostics;

import java.beans.*;

import java.awt.*;
import javax.swing.*;
import javax.swing.event.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.electrodemap.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.parameters.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * Displays the spike statistics plot. The first plot displays the spike rate and the
 * other the spike amnplitude on 2  user-selectable electrodes as well as the averages
 * over all electrodes.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class SpikeStatisticsPlot
    implements SpikeListener {

    //    private double tMax;

    private final int nElectrodes;
    final double timePerSample;
    private double timeBase;
    private double timeBinning, factor;
    private int binIndex;

    private PlotPanel spikeRatePlot, spikeAmplitudePlot;
    private int redElectrode = 0, blueElectrode = 0;
    private JInternalFrame frame;
    private int lastUpdateTime = -10 * 20000;
    private double refreshRateHz = 1;

    private HistogramStyle redStyle, blueStyle, blackStyle;

    // rates
    private int[] nSpikes;
    private DoubleHistogram[] spikeRates;
    private DoubleHistogram spikeRateAverage;

    // amplitudes
    private double[] amplitudes;
    private DoubleHistogram[] spikeAmplitudes;
    private DoubleHistogram spikeAmplitudeAverage;


    /**
     *
     * @param spikeFinder SpikeFinder
     * @param electrodeMap ElectrodeMap
     * @param timePerSample double
     * @param timeBinning double bin size in seconds
     * @param tMax double
     */
    public SpikeStatisticsPlot(
        final SpikeFinder spikeFinder, ElectrodeMap electrodeMap, double timePerSample,
        double timeBinning, double tMax) {

        this.nElectrodes = electrodeMap.getNumberOfElectrodes();
        this.timePerSample = timePerSample;
//        this.tMax = tMax;
        this.nSpikes = new int[nElectrodes];
        this.amplitudes = new double[nElectrodes];
        this.spikeRates = new DoubleHistogram[nElectrodes];
        this.spikeAmplitudes = new DoubleHistogram[nElectrodes];
        for (int i = 0; i < nElectrodes; i++) {
            spikeRates[i] = new DoubleHistogram("", 0, tMax, timeBinning);
            spikeAmplitudes[i] = new DoubleHistogram("", 0, tMax, timeBinning);
        }
        spikeRateAverage = new DoubleHistogram("", 0, tMax, timeBinning);
        spikeAmplitudeAverage = new DoubleHistogram("", 0, tMax, timeBinning);
        this.factor = 1 / timeBinning;
        this.timeBinning = 20000 * timeBinning; // convert to milliseconds
        redStyle = new HistogramStyle(
            "Hsitogram 1", HistogramStyle.OutlineType.LINEAR, Color.red, 1,
            false, Color.yellow, false, Color.black, 1);
        blueStyle = new HistogramStyle(
            "Hsitogram 2", HistogramStyle.OutlineType.LINEAR, Color.blue, 1,
            false, Color.yellow, false, Color.black, 1);
        blackStyle = new HistogramStyle(
            "Average Hsitogram", HistogramStyle.OutlineType.LINEAR, Color.black, 1,
            false, Color.yellow, false, Color.black, 1);

        // make histogram drawing fast
        redStyle.isFastAndSimple = true;
        blueStyle.isFastAndSimple = true;
        blackStyle.isFastAndSimple = true;

        this.timeBase = 0;
        this.binIndex = 0;
        this.spikeRatePlot = new PlotPanel();
        spikeRatePlot.setRange(0, tMax, 0, 100);
        spikeRatePlot.setLabels("Time (s)", "Spike Rate (spikes/s)");
        spikeRatePlot.setBackground(Color.white);

        this.spikeAmplitudePlot = new PlotPanel();
        spikeAmplitudePlot.setRange(0, tMax, 0, 500);
        spikeAmplitudePlot.setLabels("Time (s)", "Spike Amplitude (ADC)");
        spikeAmplitudePlot.setBackground(Color.white);

        // create the internal frame
        JPanel p = new JPanel(new GridLayout(2, 0));
        p.add(spikeRatePlot);
        p.add(spikeAmplitudePlot);
        frame = Vision.getInstance().createFrame(p, getController(), null,
                                                 "Spike Rate");
        frame.addInternalFrameListener(new InternalFrameAdapter() {
            public void internalFrameClosed(InternalFrameEvent e) {
                spikeFinder.removeSpikeListener(SpikeStatisticsPlot.this);
            }
        });

        spikeFinder.addSpikeListener(this);
        updateElectrodes();
    }


    public void processSpike(Spike spike) {
//        double time = spike.time * timePerSample;

        if (spike.time > timeBase + timeBinning) { // we passed over the time window
            double averageSpikeRate = 0, averageAmplitude = 0;
            for (int i = 1; i < nElectrodes; i++) {
                double n = nSpikes[i] * factor;
                spikeRates[i].fillBin(binIndex, n);
                averageSpikeRate += n;

                double amp = (nSpikes[i] != 0) ? amplitudes[i] / nSpikes[i] : 0;
                spikeAmplitudes[i].fillBin(binIndex, amp);
                averageAmplitude += amp;
            }
            averageSpikeRate /= nElectrodes;
            spikeRateAverage.fillBin(binIndex, averageSpikeRate);

            averageAmplitude /= nElectrodes;
            spikeAmplitudeAverage.fillBin(binIndex, averageAmplitude);

            // repaint the screen if needed
            if (spike.time - lastUpdateTime > 20000.0 / refreshRateHz) {
                lastUpdateTime = spike.time;
                spikeRatePlot.replotAllData();
                spikeAmplitudePlot.replotAllData();
            }

            timeBase += timeBinning;
            binIndex++;
            for (int i = 0; i < nSpikes.length; i++) {
                nSpikes[i] = 0;
                amplitudes[i] = 0;
            }
        }

        nSpikes[spike.electrode]++;
        amplitudes[spike.electrode] += spike.amplitude;
    }


    public void finishSpikeProcessing() {}


    public void updateElectrodes() {
        spikeRatePlot.removeAllData();
        spikeRatePlot.addData(spikeRates[redElectrode], redStyle);
        spikeRatePlot.addData(spikeRates[blueElectrode], blueStyle);
        spikeRatePlot.addData(spikeRateAverage, blackStyle);
        spikeRatePlot.replotAllData();

        spikeAmplitudePlot.removeAllData();
        spikeAmplitudePlot.addData(spikeAmplitudes[redElectrode], redStyle);
        spikeAmplitudePlot.addData(spikeAmplitudes[blueElectrode], blueStyle);
        spikeAmplitudePlot.addData(spikeAmplitudeAverage, blackStyle);
        spikeAmplitudePlot.replotAllData();
    }


    private JComponent getController() {
        ParametersTable t = new ParametersTable();

        t.addParameter(
            new IntegerParameter("Electrode 1", null, null, redElectrode, 0,
                                 nElectrodes - 1),
            new PropertyChangeListener() {
            public void propertyChange(PropertyChangeEvent e) {
                redElectrode = ( (IntegerParameter) e.getNewValue()).getValue();
                updateElectrodes();
            }
        });

        t.addParameter(
            new IntegerParameter("Electrode 2", null, null, blueElectrode, 0,
                                 nElectrodes - 1),
            new PropertyChangeListener() {
            public void propertyChange(PropertyChangeEvent e) {
                blueElectrode = ( (IntegerParameter) e.getNewValue()).getValue();
                updateElectrodes();
            }
        });

        t.addParameter(
            new DoubleParameter("Refresh Rate (Hz)", null, null, refreshRateHz),
            new PropertyChangeListener() {
            public void propertyChange(PropertyChangeEvent e) {
                refreshRateHz = ( (DoubleParameter) e.getNewValue()).getValue();
                if (refreshRateHz > 10) {
                    refreshRateHz = 10;
                }
            }
        });

        return new JScrollPane(t);
    }


    public void dispose() {
        Vision.getInstance().removeFrame(frame);
    }


}
