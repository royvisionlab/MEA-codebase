package edu.ucsc.neurobiology.vision.simulation;

import java.io.*;

import static java.awt.Color.*;
import java.awt.*;

import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.io.chunk.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.math.fitting.*;
import edu.ucsc.neurobiology.vision.neuronviewer.*;
import static edu.ucsc.neurobiology.vision.plot.SymbolType.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.stimulus.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class AverageSTAs {
    NeuronFile nf;
    ParametersFile pf;
    public STAFile sf;
    GlobalsFile globalsFile;

    public AverageSTAs(String fileName) throws IOException {
        nf = new NeuronFile(fileName + VisionParams.NEURON_FILE_EXTENSION);
        pf = new ParametersFile(fileName + ".params");
        sf = new STAFile(fileName + VisionParams.STA_FILE_EXTENSION);
        globalsFile = new GlobalsFile(fileName + ".globals", GlobalsFile.READ);
    }


    public void analyzeAverageSTA(String clas) throws Exception {
        int[] id = pf.getNeuronsInClass(clas);

        int w = sf.getWidth();
        int h = sf.getHeight();
        final double N = 6;
        int frameToFit = 24;

        STAFrame[] newFrames = new STAFrame[sf.getSTADepth()];
        for (int i = 0; i < newFrames.length; i++) {
            newFrames[i] = new STAFrame( (int) (w * N), (int) (h * N), 1, 1);
        }

        DoubleHistogram2D n = new DoubleHistogram2D(
            "", 0, w * N, 0, h * N, 1, 1);

        for (int k = 0; k < id.length; k++) {
            STA sta = sf.getSTA(id[k]);
            double x0 = pf.getDoubleCell(id[k], "x0");
            double y0 = pf.getDoubleCell(id[k], "y0");

            ImageFrame[] oldFrames = new ImageFrame[sf.getSTADepth()];
            for (int f = 0; f < newFrames.length; f++) {
                oldFrames[f] = sta.getFrame(f);
            }

            for (int i = 0; i < N * w; i++) {
                for (int j = 0; j < N * h; j++) {
                    int I = (int) ( (i / N - w / 2 + x0));
                    int J = (int) ( (j / N + h / 2 - y0));
                    if (I >= 0 && I < w && J >= 0 && J < h) {
                        // add this pixel to all frames
                        for (int f = 0; f < newFrames.length; f++) {
                            for (int c = 0; c < 3; c++) {
                                newFrames[f].addToPixel(i, j, c,
                                    oldFrames[f].getPixel(I, J, c));
                            }
                        }
                        n.fillBin(i, j, 1);
                    }
                }
            }
        }

        // normalize the frames
        for (int f = 0; f < newFrames.length; f++) {
            for (int i = 0; i < N * w; i++) {
                for (int j = 0; j < N * h; j++) {
                    for (int c = 0; c < 3; c++) {
                        newFrames[f].setPixel(i, j, c,
                                              (float) (newFrames[f].getPixel(i, j, c) /
                            n.getBin(i, j)));
                    }
                }
            }
        }

        // construct the STA
        STA newSTA = new STA(newFrames, 1);

        // fit the spatial RF
        double[] yy = new double[ (int) (N * w)];
        double[] xx = new double[ (int) (N * w)];
        double[] ss = new double[ (int) (N * w)];
        ImageFrame mainFrame = newSTA.getFrame(frameToFit);
        for (int i = 0; i < yy.length; i++) {
            yy[i] = mainFrame.getPixel(i, (int) (N * h / 2), 1);
            xx[i] = i / N;
            ss[i] = 0.001;
        }

        // fit the Gaussian
        double a = mainFrame.getPixel( (int) (N * w / 2), (int) (N * h / 2), 1);
        DOG1DFunction g = new DOG1DFunction( (1 + 1. / 6.) * a, -a / 6, w / 2, 0.68, 1.38);
        g.setParameterState(2, false);

        for (int i = 0; i < 4; i++) {
            // fit only the sigmas
            g.setParameterState(0, false);
            g.setParameterState(1, false);
            g.setParameterState(3, true);
            g.setParameterState(4, true);
            Fitter.fit1D(g, xx, yy, ss, yy.length);

            // fit only the amplitudes
            g.setParameterState(0, true);
            g.setParameterState(1, true);
            g.setParameterState(3, false);
            g.setParameterState(4, false);
            Fitter.fit1D(g, xx, yy, ss, yy.length);
        }

        System.err.println(g);

        PlotPanel pp = new PlotPanel();
        pp.addData(new ScatterPlot(xx, yy, null),
                   new ScatterPlotStyle(FILLED_SQUARE, 4, black));
        pp.addData(g, new FunctionStyle(""));
        pp.autoscale();
        pp.padY();
        PlotUtil.showData("Spatial Receptive Field", pp, new Rectangle(650, 20, 400, 600));

        double[] s1 = {
                      0, // center
                      3.5 * N * g.getSigma1(), // surround
                      3.5 * N * g.getSigma2(), // far surround
                      4 * N * g.getSigma2(), // far surround
        };

        double[] s2 = {
                      N * g.getSigma1(),
                      5 * N * g.getSigma1(),
                      4 * N * g.getSigma2(),
                      5 * N * g.getSigma2(),
        };

        PlotPanel p = new PlotPanel();

        for (int i = 0; i < s1.length; i++) {
            double[] TC = newSTA.getTimeFilters(
                0.5 + N * w / 2, -0.5 + N * h / 2, s1[i], s2[i], 1);

//            MathUtil.divide(TC[1], MathUtil.sumAbs(TC[1]));
            System.err.println("G" + i + ": " + MathUtil.extreme(TC));

            ScatterPlotStyle sp = new ScatterPlotStyle(
                SQUARE, 3, PlotUtil.getColor(i), true, PlotUtil.getColor(i), 2);
            p.addData(new ScatterPlot(TC), sp);
        }

        ParametricEllipse centerEllipse = new ParametricEllipse(
            0.5 + N * w / 2, -0.5 + N * h / 2, N * g.getSigma1(), N * g.getSigma1(), 0);
        ParametricEllipse surroundEllipse = new ParametricEllipse(
            0.5 + N * w / 2, -0.5 + N * h / 2, N * g.getSigma2(), N * g.getSigma2(), 0);
        PlotPanel staPanel = STAPlotMaker.makeSTAPanel(
            newSTA, false, 2, 0, false, true, false, globalsFile, centerEllipse, surroundEllipse);
        
        for (int i = 0; i < s1.length; i++) {
            staPanel.addData(getSP(newSTA,
                                   0.5 + N * w / 2, -0.5 + N * h / 2, s1[i], s2[i]),
                             new ScatterPlotStyle(FILLED_SQUARE, 2, PlotUtil.getColor(i)));
        }
        p.autoscale();
        PlotUtil.showData("Time Courses", p, new Rectangle(20, 350, 600, 600));
//            .setGridVisible(true);
        PlotUtil.showData("Average STA", staPanel, new Rectangle(20, 20, 600, 300));

        /*
                // make the radial integration curve
                double dR = 3;
                DoubleHistogram radialPlot = new DoubleHistogram("", 0, N * w, dR);
                for (int y = 0; y < N * h; y++) {
                    for (int x = 0; x < N * w; x++) {
                        double R = Math.sqrt( (x - N * w / 2) * (x - N * w / 2) +
                                             (y - N * h / 2) * (y - N * h / 2));
                        radialPlot.fill(R, mainFrame.getPixel(x, y, 1));
                    }
                }
                PlotPanel p1 = PlotUtil.showData(radialPlot, new HistogramStyle());
                p1.setYRange( -15, 15);
                p1.setGridVisible(true);
         */
    }


    public static ScatterPlot getSP(STA sta, double x0, double y0, double s1, double s2) {
        ScatterPlot sp = new ScatterPlot();

        int height = sta.getHeight();
        int width = sta.getWidth();
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                double d = Math.sqrt( (x + 0.5 - x0) * (x + 0.5 - x0) +
                                     (y + 0.5 - y0) * (y + 0.5 - y0));
                if (d >= s1 && d < s2) {
                    sp.add(x + 0.5, y + 0.5);
                }
            }
        }

        return sp;
    }


    public static void main(String[] args) throws Exception {
        final AverageSTAs ntest = new AverageSTAs(
            "f:\\Good Data\\2005-04-26-0\\data009\\data009");

        ntest.analyzeAverageSTA("All/OFF/BTY");
    }
}
