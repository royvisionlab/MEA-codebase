package edu.ucsc.neurobiology.vision.onlinediagnostics;

import java.awt.*;
import javax.swing.*;
import javax.swing.event.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.parameters.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class TTLPeriodPlot
    implements SpikeListener {

    private double timeBase = 0;
    private double timeWindow = 4000;
    private DoubleHistogram ttlIntervalHist;
    private PlotPanel plotTTL;
    private IntegerList ttlIntervals;
    private int ttlTime = 0;
    private int nTTLs = 0;
    private double nCleanTTLIntervals = 0;
    private int ttlCapacity = 100;
    private int yRange = 2000;
    private int maxTTL = 0;
    HistogramStyle style = new HistogramStyle("TTL Period");
    private JInternalFrame frame;


    public TTLPeriodPlot(final SpikeFinder spikeFinder) {
        ttlIntervals = new IntegerList();

        this.plotTTL = new PlotPanel();
        plotTTL.setLabels("TTL Number", "Period (Samples)");
        plotTTL.setBackground(Color.white);

        ttlIntervalHist = new DoubleHistogram("TTL Intervals", 0, ttlCapacity, 1);
        plotTTL.addData(ttlIntervalHist, style);

        frame = Vision.getInstance().createFrame(plotTTL, getController(), null,
                                                 "TTL Period");
        frame.addInternalFrameListener(new InternalFrameAdapter() {
            public void internalFrameClosed(InternalFrameEvent e) {
                spikeFinder.removeSpikeListener(TTLPeriodPlot.this);
            }
        });

        spikeFinder.addSpikeListener(this);
    }


    public void processSpike(Spike spike) {
        if (spike.electrode == 0) {
            int dt = spike.time - ttlTime;
            ttlIntervals.add(dt);
            //used to find middle of plot, exclude first TTL because the period tends to be strange
            if (nTTLs > 0) {
                maxTTL = Math.max(maxTTL, dt);
            }
            ttlTime = spike.time;

            nTTLs++;

            // we passed over the time window
            if (ttlTime > timeBase + timeWindow) {
                if (nTTLs >= ttlCapacity) {
                    ttlCapacity *= 2;

                    DoubleHistogram ttlIntervalHistNew = new DoubleHistogram(
                        "TTL Intervals", 0, ttlCapacity, 1);
                    for (int i = 0; i < ttlIntervalHist.getBinCount(); i++) {
                        ttlIntervalHistNew.setBin(i, ttlIntervalHist.getBin(i));
                    }
                    ttlIntervalHist = ttlIntervalHistNew;
                    plotTTL.removeAllData();
                    plotTTL.addData(ttlIntervalHist, style);
                }

                ttlIntervalHist.fillBin(nTTLs, dt);

                plotTTL.setRange(0, ttlCapacity, maxTTL - yRange / 2,
                                 maxTTL + yRange / 2);
                plotTTL.replotAllData();

                timeBase += timeWindow;
            }
        }
    }


    public void finishSpikeProcessing() {}


    private JComponent getController() {
        ParametersTable t = new ParametersTable();

//        t.addParameter(
//            new IntegerParameter("Electrode(blue)", 0, 0, nElectrodes - 1),
//            new PropertyChangeListener() {
//                public void propertyChange(PropertyChangeEvent e) {
//                    blueElectrode = ((IntegerParameter)e.getNewValue()).getValue();
//                    updateElectrodes();
//                }
//            }
//        );
        return new JScrollPane(t);
    }


    public void dispose() {
        Vision.getInstance().removeFrame(frame);
    }


}
