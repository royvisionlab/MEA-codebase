package edu.ucsc.neurobiology.vision.onlinediagnostics;

import java.util.*;

import java.awt.*;
import javax.swing.*;
import javax.swing.event.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.plot.*;


/**
 * Displays two summary plots; the first represents the mean of the signal on every
 * electrode; the second represents the r.m.s of the signal on every electrode.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class MeansAndSigmasDiagnosticPlot
    implements SampleListener {

    private final int nElectrodes;
    private float[] sum, sum2, means, sigmas;
    private int nSamples;
    private DoubleHistogram sigmasPlot, meansPlot;
    private PlotPanel sigmasPanel, meansPanel;
    private final float norm = (float) (1.0 / 20000.0);
    HistogramStyle meansStyle = new HistogramStyle("Means");
    HistogramStyle sigmasStyle = new HistogramStyle("Sigmas");
    private JInternalFrame frame;


    public MeansAndSigmasDiagnosticPlot(
        final MultipleCompressedSampleInputStream sampleInputStream, int nElectrodes) {

        this.nElectrodes = nElectrodes;
        sum = new float[nElectrodes];
        sum2 = new float[nElectrodes];
        means = new float[nElectrodes];
        sigmas = new float[nElectrodes];

        meansPlot = new DoubleHistogram("", 0, nElectrodes, 1);
        meansPanel = new PlotPanel();
        meansPanel.setRange(0, nElectrodes, -500, 500);
        meansPanel.setLabels("Electrode", "Mean (ADC)");
        meansPanel.addData(meansPlot, meansStyle);

        sigmasPlot = new DoubleHistogram("", 0, nElectrodes, 1);
        sigmasPanel = new PlotPanel();
        sigmasPanel.setRange(0, nElectrodes, 0, 100);
        sigmasPanel.setLabels("Electrode", "Sigma (ADC)");
        sigmasPanel.addData(sigmasPlot, sigmasStyle);

        nSamples = 0;

        JPanel p = new JPanel(new GridLayout(2, 1));
        p.add(meansPanel);
        p.add(sigmasPanel);
        frame = Vision.getInstance().createFrame(p, null, null, "Means & Sigmas");
        frame.addInternalFrameListener(new InternalFrameAdapter() {
            public void internalFrameClosed(InternalFrameEvent e) {
                sampleInputStream.removeSampleListener(MeansAndSigmasDiagnosticPlot.this);
            }
        });

        sampleInputStream.addSampleListener(this);
    }


    float y;
    int i;
    public final void processSample(short[] sample) {
        for (i = 0; i < nElectrodes; i++) {
            y = sample[i];
            sum[i] += y;
            sum2[i] += y * y;
        }

        nSamples++;
        if (nSamples % 20000 == 0) { //the end of a 1s period
            sigmasPlot.clear();
            meansPlot.clear();
            for (int i = 1; i < nElectrodes; i++) {
                means[i] = sum[i] * norm;
                meansPlot.fillBin(i, means[i]);

                sigmas[i] = (float) Math.sqrt( (sum2[i] - norm * sum[i] * sum[i]) /
                                              (20000 - 1.));
                if (!Double.isNaN(sigmas[i])) {
                    sigmasPlot.fillBin(i, sigmas[i]);
                }
            }
            sigmasPanel.replotAllData();
            meansPanel.replotAllData();
            Arrays.fill(sum, 0);
            Arrays.fill(sum2, 0);
        }
    }


    public void dispose() {
        Vision.getInstance().removeFrame(frame);
    }


    public void finishSampleProcessing() {
    }
}
