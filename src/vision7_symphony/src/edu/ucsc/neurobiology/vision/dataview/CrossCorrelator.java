package edu.ucsc.neurobiology.vision.dataview;

import java.beans.*;
import java.io.*;

import java.awt.*;
import javax.swing.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.analysis.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.parameters.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class CrossCorrelator
    extends JPanel {

    private NeuronFile neuronFile;
    private PlotPanel panel;
    private HistogramStyle style;
    private int[] idList;
    private int[] t1, t2;
    private int id1, id2;

    private double binWidth = 0.25; // in ms
    private double deltaT = 50; // in ms
    private double coincidenceTime = 0.5; // in ms


    public CrossCorrelator(String neuronFileName) throws IOException {
        super(new BorderLayout());

        neuronFile = new NeuronFile(neuronFileName);
        idList = neuronFile.getIDList();

        style = new HistogramStyle();
        panel = new PlotPanel();
        panel.setLabels("Time Difference (ms)", "Coincidence Fraction");

        this.add(panel, BorderLayout.CENTER);
        Vision.getInstance().createFrame(
            this, getController(), null, "Cross Correlator - " + neuronFileName);
    }


    private JComponent getController() {
        ParametersTable t = new ParametersTable();

        final IntegerParameter p1 = new IntegerParameter("Neuron 1", null, null, 0, 0,
            neuronFile.getNumberOfNeurons() - 1);
        final IntegerParameter p2 = new IntegerParameter("Neuron 2", null, null, 0, 0,
            neuronFile.getNumberOfNeurons() - 1);

        PropertyChangeListener listener = new PropertyChangeListener() {
            public void propertyChange(PropertyChangeEvent e) {
                id1 = idList[p1.getValue()];
                id2 = idList[p2.getValue()];
                set(id1, id2);
            }
        };

        t.addParameter(p1, listener);
        t.addParameter(p2, listener);

        final DoubleParameter binWidthParam = new DoubleParameter(
            "Bin Width", null, null, binWidth);
        t.addParameter(binWidthParam, new PropertyChangeListener() {
            public void propertyChange(PropertyChangeEvent e) {
                if (binWidthParam.getValue() > 0) {
                    binWidth = binWidthParam.getValue();
                    set(id1, id2);
                }
            }
        });

        final DoubleParameter timeIntervalParam = new DoubleParameter(
            "Time Interval", null, null, deltaT);
        t.addParameter(timeIntervalParam, new PropertyChangeListener() {
            public void propertyChange(PropertyChangeEvent e) {
                if (timeIntervalParam.getValue() > 0) {
                    deltaT = timeIntervalParam.getValue();
                    set(id1, id2);
                }
            }
        });

        final DoubleParameter coincidenceTimeParam = new DoubleParameter(
            "Coincidence Time", null, null, coincidenceTime);
        t.addParameter(coincidenceTimeParam, new PropertyChangeListener() {
            public void propertyChange(PropertyChangeEvent e) {
                if (coincidenceTimeParam.getValue() > 0) {
                    coincidenceTime = coincidenceTimeParam.getValue();
                    set(id1, id2);
                }
            }
        });

        return new JScrollPane(t);
    }


    public void set(int id1, int id2) {
        this.id1 = id1;
        this.id2 = id2;

        try {
            t1 = neuronFile.getSpikeTimes(id1);
            t2 = neuronFile.getSpikeTimes(id2);
        } catch (IOException ex) {
            ex.printStackTrace();
        }

        ///////////////////////////////////
        DoubleHistogram ccH = CrossCorrelationCalculator.getCrossCorrelationHistogram(
            t1, t2, (int) (binWidth * 20), (int) (deltaT * 20), null);

        // find those bins that have maximum sum.
        double[] cc = ccH.toArray();
        double maxV = Double.NEGATIVE_INFINITY;
        for (int k = 0; k < cc.length - 2 * coincidenceTime; k++) {
            double v = 0;
            for (int m = 0; m < 2 * coincidenceTime; m++) {
                v += cc[k + m];
            }
            if (v > maxV) {
                maxV = v;
            }
        }
        double p = maxV / Math.min(t1.length, t2.length);
        ///////////////////////////////////


        // calculate coupling strength
        int nBins = ccH.getBinCount();
        int n = 0;
        double mean = 0;
        for (int i = 0; i < nBins / 4; i++) {
            mean += ccH.getBin(i);
            n++;
        }
        for (int i = nBins * 3 / 4; i < nBins; i++) {
            mean += ccH.getBin(i);
            n++;
        }
        mean /= n;

        int nCoupling = (int) (ccH.getBinSum() - mean * nBins);
        double coupling = nCoupling / Math.sqrt(t1.length * t2.length);

        panel.removeAllData();

        DoubleHistogram ccH1 = new DoubleHistogram("", -deltaT, +deltaT, binWidth);
        for (int i = 0; i < ccH.getBinCount(); i++) {
            ccH1.setBin(i, ccH.getBin(i));
        }
        panel.addData(ccH1, style);

        ScatterPlot sp = new ScatterPlot();
        double y = ( (double) t1.length * t2.length) /
                   (neuronFile.getNumberOfSamples() / 20.0 / binWidth);
        sp.add(ccH1.getMin(), y);
        sp.add(ccH1.getMax(), y);
        panel.addData(sp, "NONE 0 black, SOLID 2 red");

        panel.autoscale();
        panel.clearLegend();

        panel.addToLegend("ID " + id1 + " - " + t1.length + " spikes");
        panel.addToLegend("ID " + id2 + " - " + t2.length + " spikes");
        panel.addToLegend( (int) maxV + " spikes in common");
        panel.addToLegend("Correlation: " + StringUtil.format(p, 3));
//        panel.addToLegend("--------------");
//        panel.addToLegend("mean " + mean);
//        panel.addToLegend("nCoupling " + nCoupling);
//        panel.addToLegend("Coupling " + StringUtil.format(coupling, 2));
    }


    /**
     *
     * @param t1 int[]
     * @param t2 int[]
     * @param binWidth int in samples
     * @param deltaT int in samples
     * @return double
     */
    public static double getCoupling(int[] t1, int[] t2, int binWidth, int deltaT) {
        DoubleHistogram ccH = CrossCorrelationCalculator.getCrossCorrelationHistogram(
            t1, t2, binWidth, deltaT, null);

        // calculate coupling strength
        int nBins = ccH.getBinCount();
        int n = 0;
        double mean = 0;
        for (int i = 0; i < nBins / 4; i++) {
            mean += ccH.getBin(i);
            n++;
        }
        for (int i = nBins * 3 / 4; i < nBins; i++) {
            mean += ccH.getBin(i);
            n++;
        }
        mean /= n;

        int nCoupling = (int) (ccH.getBinSum() - mean * nBins);
        return nCoupling / Math.sqrt(t1.length * t2.length);
    }


}
