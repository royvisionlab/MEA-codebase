package edu.ucsc.neurobiology.vision.simulation;

import java.io.*;
import static java.lang.Math.*;

import java.awt.*;
import javax.swing.*;

import edu.ucsc.neurobiology.vision.analysis.*;
import edu.ucsc.neurobiology.vision.gui.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.math.fitting.*;
import edu.ucsc.neurobiology.vision.plot.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class LinearNonlinearModel {

    boolean DEBUG = false;
    JFrame window = new JFrame();
    PlotPanel layoutPanel = new PlotPanel();
    FunctionStyle st = new FunctionStyle("f");

    ScatterPlotStyle blackThin = new ScatterPlotStyle(
        SymbolType.NONE, 1, Color.black, true, Color.black, 1);

    ScatterPlotStyle redThin = new ScatterPlotStyle(
        SymbolType.NONE, 1, Color.black, true, Color.red, 1);
    ScatterPlotStyle redThick = new ScatterPlotStyle(
        SymbolType.NONE, 1, Color.black, true, Color.red, 2);

    ScatterPlotStyle blueThin = new ScatterPlotStyle(
        SymbolType.NONE, 1, Color.black, true, Color.blue, 1);
    ScatterPlotStyle blueThick = new ScatterPlotStyle(
        SymbolType.NONE, 1, Color.black, true, Color.blue, 2);

    ScatterPlotStyle style3 = new ScatterPlotStyle(
        SymbolType.NONE, 1, Color.black, true, Color.green, 2);


    public LinearNonlinearModel() throws IOException {
    }


    public static double[] dotProduct(double[] filter, double[] data) {
        double[] result = new double[data.length];

        for (int n = 0; n < data.length; n++) {
            int firstIndex = n - filter.length + 1;
            for (int i = 0; i < filter.length; i++) {
                if (firstIndex + i >= 0) {
                    result[n] += filter[i] * data[firstIndex + i];
                }
            }
        }

        return result;
    }


    public static void dotProduct1(float[][] staFrames, float[][] movie, double[] current) {
        for (int n = 0; n < movie.length; n++) {
            current[n] = 0;

            // calculate the instantaneous current
            int firstFrame = n - staFrames.length + 1;
            for (int f = 0; f < staFrames.length; f++) {
                if (firstFrame + f >= 0) {
                    for (int i = 0; i < staFrames[f].length; i++) {
                        current[n] += staFrames[f][i] * movie[firstFrame + f][i];
                    }
                }
            }
        }
    }


    public static DoubleHistogram convolute(STATimeFunction1 centerT,
                                            STATimeFunction1 surroundT,
                                            DoubleHistogram stimulus) {

        double dt = stimulus.getBinInterval();
        DoubleHistogram g = new DoubleHistogram(
            "", stimulus.getMin(), stimulus.getMax(), dt);
        System.err.println(g.getBinCount());

        double[] tcCenter = new double[ (int) (0.5 / dt)];
        double[] tcSurround = new double[ (int) (0.5 / dt)];
        int n = tcCenter.length;
        for (int i = 0; i < n; i++) {
            double t = -0.5 + i * dt;
            tcCenter[i] = centerT.getValueAt(t);
            tcSurround[i] = surroundT.getValueAt(t);
        }

//        PlotUtil.showArray("TC", tc);

        for (int i = n; i < g.getBinCount(); i++) {
            double v = 0;
            for (int j = 0; j < n; j++) {
                try {
                    v += tcCenter[j] * stimulus.getBin(i - n + j);
                    v += tcSurround[j] * stimulus.getBin(i - n + j);
                } catch (Exception e) {
                    System.err.println(i);
                    System.err.println(j);
                    System.exit(1);
                }
            }

            // rectify
            double R0 = 40;
            double k = 50;

            v = (v < -R0 / k) ? 0 : k * v + R0; // rectify
//            v = k * v; //current

            g.fillBin(i, v);
        }

        return g;
    }


    public static void main(String[] args) throws IOException {
        int id = 31;

        // parasol
//        double t1 = -36 / 1000.0;
//        double a1 = -1;
//        int n1 = 20;
//
//        double t2 = -80 / 1000.0;
//        double a2 = 0.17;
//        int n2 = 5;

        // BTY
        double t1 = -36 / 1000.0;
        double a1 = -1;
        int n1 = 22;

        double t2 = -72 / 1000.0;
        double a2 = +0.16;
        int n2 = 4;

        double ff = 1.5;
        double surroundWeight = 0.60;
        double farSurroundWeight = 0.40;

        double dt = 0.00834;

        double n = abs( -a1 * pow(E / n1, n1) * MathUtil.fact(n1) * t1
                       - a2 * pow(E / n2, n2) * MathUtil.fact(n2) * t2);
        STATimeFunction1 centerT = new STATimeFunction1(
            a1 / n, t1, a2 / n, t2, n1, n2);

        n = abs( -a1 * pow(E / n1, n1) * MathUtil.fact(n1) * t1 * ff
                - a2 * pow(E / n2, n2) * MathUtil.fact(n2) * t2 * ff);
        STATimeFunction1 surroundT = new STATimeFunction1(
            -a1 * surroundWeight / n, t1 * ff, -a2 * surroundWeight / n, t2 * ff, n1, n2);

        PlotPanel p = new PlotPanel();
        p.setRange( -0.5, 0, -4, 4);
        p.addData(centerT, new FunctionStyle("center"));
        p.addData(surroundT, new FunctionStyle("surround", Color.red, 1));
//        PlotUtil.showData("TC", p);

        DoubleHistogram stimulus = new DoubleHistogram("", 0, 8, dt);
        stimulus.fill(0, 2, 1);
        stimulus.fill(2, 4, 0);
        stimulus.fill(4, 6, -1);
        stimulus.fill(6, 8, 0);

        DoubleHistogram g = convolute(centerT, surroundT, stimulus);

        FlashesClassification flashesClassification = new FlashesClassification(
            new String[] {"data010"}, "f:\\data\\2005-04-26-0", 150, dt);
        DoubleHistogram data = flashesClassification.getFlashResponseHistogram(id);

        DoubleHistogram data1 = new DoubleHistogram("", data.getMin(), data.getMax(),
            data.getBinInterval());
        for (int i = 0; i < data.getBinCount() - 1; i++) {
            data1.setBin(i, data.getBin(i + 1) - data.getBin(i));
        }
        PlotUtil.showData("TC", data1);

        JFrame f = new JFrame();
        f.setLayout(new PlaceLayout());
        PlotPanel p1 = new PlotPanel();
        p1.addData(data, new HistogramStyle("data"));
        p1.addData(g, new HistogramStyle("current", HistogramStyle.OutlineType.LINEAR,
                                         Color.red, 1, false,
                                         Color.white, false, Color.white, 1));
        p1.autoscale();
        p1.padY();
        p1.setXAxisVisible(false);
//        p1.setYAxisVisible(false);
        PlotPanel p2 = new PlotPanel();
        p2.addData(stimulus, new HistogramStyle());
        p2.autoscale();
//        p2.setYAxisVisible(false);
        p2.padY();
        f.add(p1, new PlaceC(0, 0, 0, 0, 1, 0.8));
        f.add(p2, new PlaceC(0, 0.8, 0, 0, 1, 1 - 0.8));
        f.setBounds(100, 100, 600, 400);
        f.setVisible(true);
    }
}
