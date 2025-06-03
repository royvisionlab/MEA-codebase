package edu.ucsc.neurobiology.vision.test;

import java.io.*;
import java.text.*;
import java.util.*;

import java.awt.*;
import java.awt.geom.*;
import javax.swing.*;

import edu.ucsc.neurobiology.vision.electrodemap.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.neuronviewer.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ResponseVisualization {
    int nPlots;
    DoubleHistogram[][] spikeHist;
    double[][] x, y;
    double[] max;
    final double rMax = 100, rMax1 = 150;
    final double spikeRateBinning = 0.1;
    String dir = "visualization";


    public ResponseVisualization(String[] classes) throws IOException {
        nPlots = classes.length + 1;
        System.out.println("nPlots: " + nPlots);
        spikeHist = new DoubleHistogram[nPlots][];
        x = new double[nPlots][];
        y = new double[nPlots][];
        max = new double[nPlots];

        String name = "2003-06-10-0\\data014\\data014";

        // set up the electrode array data
        SpikeFile sf = new SpikeFile(name + ".spikes");
        spikeHist[0] = SpikeRatePlotMaker.getSpikeRateHistogram(sf, spikeRateBinning);
//        PlotUtilities.showData("", spikeHist[0][1], new BasicHistogramStyle());
        ElectrodeMap map = ElectrodeMapFactory.getElectrodeMap(sf.getArrayID());
        x[0] = new double[map.getNumberOfElectrodes()];
        y[0] = new double[map.getNumberOfElectrodes()];
        Point2D.Double p = new Point2D.Double(0, 0);
        for (int i = 0; i < map.getNumberOfElectrodes(); i++) {
            map.getPosition(i, p);
            x[0][i] = -p.x;
            y[0][i] = p.y;
        }

        // set up the neuron data
        NeuronFile nf = new NeuronFile(name + VisionParams.NEURON_FILE_EXTENSION);
        ParametersFile pf = new ParametersFile("2003-06-10-0\\data003\\data003.params");
        HashMap<Integer,? extends Object> classID = pf.getClassIDs();
        HashMap<Integer,Double> xMap = null;
        HashMap<Integer,Double> yMap = null;
        try {
            xMap = pf.evaluate("x", "");
            yMap = pf.evaluate("y", "");
        } catch (CannotEvaluateException ex) {
            ex.printStackTrace();
        }
        for (int plotIndex = 1; plotIndex < nPlots; plotIndex++) {
            ArrayList h = new ArrayList();
            DoubleList xList = new DoubleList();
            DoubleList yList = new DoubleList();
            for (Integer id : classID.keySet()) {
                if (classID.get(id).equals(classes[plotIndex - 1])) {
//                    h.add(SpikeRatePlotMaker.getSpikeRateHistogram(nf, id.intValue(), 0,
//                        spikeRateBinning));
                    xList.add(xMap.get(id).doubleValue() * 60);
                    yList.add(yMap.get(id).doubleValue() * 60);
                }
            }
            spikeHist[plotIndex] = (DoubleHistogram[]) h.toArray(new DoubleHistogram[h.
                size()]);
            x[plotIndex] = xList.toArray();
            y[plotIndex] = yList.toArray();
        }

        // find the calibration coeffs
        for (int plotIndex = 0; plotIndex < nPlots; plotIndex++) {
            max[plotIndex] = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < spikeHist[plotIndex].length; i++) {
                double m = spikeHist[plotIndex][i].getMaxValue();
                if (m > max[plotIndex]) {
                    max[plotIndex] = m;
                }
            }
            System.out.println("max[" + plotIndex + " = " + max[plotIndex]);
        }

        System.out.println("Initialization done");

        DecimalFormat format = new DecimalFormat();
        format.setMinimumIntegerDigits(3);
        format.setMaximumIntegerDigits(3);

        JPanel panel = new JPanel();
        panel.setLayout(new GridLayout(0, 1));
        PlotPanel[] plot = new PlotPanel[nPlots];
        for (int i = 0; i < nPlots; i++) {
            plot[i] = new PlotPanel();
            if (i == 0) {
                plot[i].setRange( -1000, 1000, -500, 500);
            } else {
                plot[i].setRange(920 - 120, 2920 - 120, 460, 1460);
            }
//            plot[i].setSize(600, 300);
            plot[i].setAxisVisible(false);

            panel.add(plot[i]);
        }
//        panel.setSize(300, 450);
        JFrame fr = new JFrame();
        fr.add(panel);
        fr.setBounds(0, 0, 500, 750);
        fr.setVisible(true);
        fr.setVisible(false);

        FunctionStyle style = new FunctionStyle("f");
        style.setLineThickness(2);
        FunctionStyle style1 = new FunctionStyle("f");
        style1.setLineThickness(2);

        File d = new File(dir);
        d.mkdir();

//        for (int frame = 29*2; frame < 2*29 + 29*2; frame++) {
//        for (int frame = 2*33; frame < 2*33 + 33*2; frame++) {
        for (int frame = 20; frame < 20 + 20 * 2; frame++) {
//        for (int frame = 14; frame < 14 + 14*2; frame++) {
            for (int plotIndex = 0; plotIndex < nPlots; plotIndex++) {
                plot[plotIndex].removeAllData();

                for (int neuron = 0; neuron < spikeHist[plotIndex].length; neuron++) {
                    double r;
                    if (plotIndex == 0) {
                        r = spikeHist[plotIndex][neuron].getBin(frame) * rMax /
                            max[plotIndex];
                    } else {
                        r = spikeHist[plotIndex][neuron].getBin(frame) * rMax1 /
                            max[plotIndex];
                        r = Math.max(r, 4);
                    }

                    ParametricEllipse ellipse = new ParametricEllipse(
                        x[plotIndex][neuron], y[plotIndex][neuron], r, r, 0);

                    if (plotIndex == 0) {
                        plot[plotIndex].addData(ellipse, style);
                    } else {
                        plot[plotIndex].addData(ellipse, style1);
                    }
                }
            }

            File f = new File(dir + File.separator + "image" +
                              format.format(frame) + ".png");

//                FileOutputStream fos = new FileOutputStream(f);
//            panel.sav.exportToFile(f, panel, null, null, "");
//                png.exportToFile(fos, panel, f, null);
//                fos.close();
        }
    }


    public static void main(String[] args) throws IOException {
        new ResponseVisualization(new String[] {"OFFDS2", "OFFDS3"});
//        new ResponseVisualization(new String[] {"OFFDS1", "OFFDS4", "OFFDS2", "OFFDS3"});
    }
}
