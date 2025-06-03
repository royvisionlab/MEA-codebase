package edu.ucsc.neurobiology.vision.test;

import java.io.*;
import static java.lang.Math.*;

import static java.awt.Color.*;
import java.awt.*;
import javax.swing.*;

import edu.ucsc.neurobiology.vision.analysis.*;
import edu.ucsc.neurobiology.vision.gui.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.math.*;
import static edu.ucsc.neurobiology.vision.neuronviewer.ShowClassificationPlotsAction.*;
import edu.ucsc.neurobiology.vision.neuronviewer.*;
import static edu.ucsc.neurobiology.vision.plot.PlotPanel.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.*;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class LargeCellPaper {
//    static String[] cName = {"All/OFF/Midget", "All/OFF/Parasol", "All/OFF/Amacrine"};
    static String[] cName = {"All/OFF/Midget", "All/OFF/Parasol", "All/OFF/BTY"};
    static double w = 8.2 / 2.54 / 2;
    static double h = 8.2 / 2.54 / 4;
    static int bx = 1, by = 1;
    static String[] className = {"Midget", "Parasol", "BTY"};
    static Font letterFont = Font.decode("Arial Bold 10");
    static Font textFont = Font.decode("Arial Italic 8");
    static Font labelFont = Font.decode("Arial Italic 7");
    static int binsPerPeriod = 30;


    public static void makeFigure3(String[] pName, String[] crgName) throws Exception {
        double f1f2Width = 1.7;
        int n = 3;
        int t = 0;

        PlotPanel[] overall = new PlotPanel[3];
        for (int i = 0; i < overall.length; i++) {
            overall[i] = new PlotPanel();
            overall[i].loadStyle("f1f2", new Config("config.xml"));
        }

        // save F1, F2 plots
        for (int i = 0; i < 3; i++) {
            // ReversingGratings m = new ReversingGratings(crgName[i], binsPerPeriod);
            ParametersFile pf = new ParametersFile(pName[i]);
            double[] freqs = pf.getArrayCell(pf.getIDList()[0], "reversingFrequencies");

            for (int c = 0; c < cName.length; c++) {
                PlotPanel p = ReversingGratingsPlotMaker.getHarmonicsPanel(
                    pf.getNeuronsInClass(cName[c]), pf, true, true, true, new Config("config.xml"));

                double[][] f1 = pf.getArrayAverage1("T1reversingF1", cName[c]);
                double[][] f2 = pf.getArrayAverage1("T1reversingF2", cName[c]);

                ScatterPlotStyle style1 = new ScatterPlotStyle("Fodd " + c,
                    SymbolType.DISK, 1, Color.black, true, Color.black, 0.25f);
                ScatterPlotStyle style2 = new ScatterPlotStyle("Feven " + c,
                    SymbolType.DISK, 1, Color.red, true, Color.red, 0.25f);

                if (i == 1) {
                    style1.setDashPattern("1.5, 1");
                    style2.setDashPattern("1.5, 1");
                } else if (i == 2) {
                    style1.setDashPattern("5, 2");
                    style2.setDashPattern("5, 2");
                }

                overall[c].addData(new ScatterPlot(freqs, f1[0], f1[1]), style1);
                overall[c].addData(new ScatterPlot(freqs, f2[0], f2[1]), style2);

                if (i == 0 && c == 0) {
                    t = p.getYAxis().getThickness();
                }

//                                // save
//                                p.setXAxisVisible(true);
//                                p.setYAxisVisible(c == 0);
//
//                 PlotUtil.showData(pName[i] + ": " + cName[c], p, 800, 400);
//                 p.addBackgroundText("" + (char) ('a' + n), RIGHT, TOP, letterFont);
//                                if (i == 0 && c == 0) {
//                 p.addBackgroundText("<html>F<sub>1</sub>", CENTER, 10, textFont);
//                 p.addBackgroundText("<html>F<sub>2</sub>", CENTER, 80, textFont);
//                                }
//
//                                p.setLabelFont(labelFont);
//                                p.getYAxis().setSpacings( -1, 1);
//                 p.saveAsEPS(new File("fig3-" + n++ +".eps"), f1f2Width, 1.2);
            }
            pf.close(false);
        }

        // save overall F1, F2 plots
        int kk = 3;
        for (int c = 0; c < overall.length; c++) {
            overall[c].setXAxisVisible(true);
            overall[c].setYAxisVisible(c == 0);
            overall[c].addBackgroundText("" + (char) ('a' + kk++), RIGHT, TOP, letterFont);
            if (c == 0) {
                overall[c].getYAxis().setThickness(30);
                overall[c].addBackgroundText("<html>F<sub>1</sub>", CENTER, 5, textFont);
                overall[c].addBackgroundText("<html>F<sub>2</sub>", CENTER, 65, textFont);
            }
            overall[c].setLabelFont(labelFont);
            overall[c].getYAxis().setSpacings( -1, 1);
            overall[c].setYRange(1, 115);
            overall[c].setLabels("Spatial Frequency (cyc/mm)", "Response (spikes/s)");

            overall[c].saveAsEPS(new File("fig3-overall-" + c + ".eps"), f1f2Width, 1.2, true);
        }

        // save the example response figures
        ReversingGratings m = new ReversingGratings(crgName[0], binsPerPeriod);

        int[] id = {1899, 621, 922};
        int[] max = {100, 250, 150};
        int[] startPhase0 = {2, 0, 0};
        int[] startPhase6 = {5, 7, 0};

        n = 0;
        GridLayout g = new GridLayout(2, 4);
        g.setHgap(1);
        g.setVgap(1);
        for (int i = 0; i < id.length; i++) {
            m.setCurrentNeuron(id[i]);

            JPanel p = new JPanel(g);
            p.setBackground(white);

            for (int k = 0; k < 4; k++) {
                p.add(m.getPSTH(0, (startPhase0[i] + k * 2) % 8, max[i]));
            }
            for (int k = 0; k < 4; k++) {
                p.add(m.getPSTH(6, (startPhase6[i] + k * 2) % 8, max[i]));
            }

            JFrame f = new JFrame();
            f.add(p);
            f.setVisible(true);
            GraphicsIO.saveAsEPS(p, "fig3-" + n++ +".eps", f1f2Width * 1.01, 0.7, false);
        }
    }


    private static int saveTCPlot(NeuronViewer nv, int n, double y1, double y2,
                                  double ySpacing, char letter) throws IOException {
        for (int c = 1; c < cName.length; c++) {
            nv.select(cName[c]);
            PlotPanel p = (PlotPanel) nv.getPlot("tc");
            if (p == null) {
                continue;
            }

            p.setLabelFont(labelFont);
            p.setRange( -199, 0, y1, y2);
            p.setLabels("Time before spike (ms)", "STA contrast (a.u.)");
            p.getXAxis().shiftRightLabel = true;
            if (c == 1) {
                p.addBackgroundText("" + letter, LEFT, 2, letterFont);
            }

            // create S vs. LM plots
//            p.removeDataWithStyle(p.getStyleWithName("Red Time Course"));
//            ScatterPlotStyle s = (ScatterPlotStyle) p.getStyleWithName(
//                "Green Time Course");
//            s.setConnectionLineColor(Color.orange);

            p.getYAxis().setFixedTickSpacing(0, ySpacing, 2);
            p.getYAxis().multiplier = 50;
            p.getYAxis().setSpacings(0, 1);

            p.axesBorder.setTopPadding(3);
            p.getXAxis().setFixedTickSpacing(0, 50, 0);
            p.getXAxis().tickMarkAllignment = HorizontalAxis.RIGHT_TICK_ALLIGNMENT;
            p.getXAxis().setMarkToLabelSpacing(0);
            p.saveAsEPS(new File("fig2-" + n++ +".eps"), w, 1.2 * h, false);
        }

        return n;
    }


    public static void makeFigure2(String[] pName) throws Exception {
        int n = 0;

        for (int i = 0; i < 3; i++) {
            NeuronViewer nv = new NeuronViewer(pName[i], false, new Config("config.xml"));

            // save example STA's
            if (i == 0) {
                String[] id = {"All/OFF/Parasol/1058", "All/OFF/BTY/675"};
                for (int k = 0; k < id.length; k++) {
                    nv.select(id[k]);
                    PlotPanel p = (PlotPanel) nv.getPlot("sta");

                    FunctionStyle s = (FunctionStyle) p.getStyleWithName("Gaussian Fit");
                    s.setLineColor(white);

                    p.setLegendVisible(false);
                    if (k == 0) {
                        p.addBackgroundText("" + 'a', LEFT, TOP, letterFont);
                    }
                    p.addBackgroundText(className[k + 1], RIGHT, TOP, textFont);
                    p.setScaleBar(500 / 5.8, 1, 0.05, 0.1);
                    p.axesBorder.setBottomPadding(3);
                    p.setDrawBoxLines(false, false, false, false);
                    p.saveAsEPS(new File("fig2-" + n++ +".eps"), w, h);
                    p.setLegendVisible(true);
                }
            }

            // save mosaics
            for (int c = 1; c < cName.length; c++) {
                nv.select(cName[c]);

                PlotPanel p = (PlotPanel) nv.getPlot("mosaic");
                p.setLegendVisible(false);
                if (c == 1) {
                    p.addBackgroundText("" + (char) ('a' + i + 1), LEFT, TOP, letterFont);
                }

                Num d = nv.getParametersFile().getAverage1(
                    "2*5.8*20*((SigmaX*SigmaY)^0.5)", cName[c]);

                if (i == 0 && c == 1) {
                    p.addBackgroundText("<html>RF Diameter = " + d.toString(0) +
                                        "\u03bcm",
                                        RIGHT, 0, Font.decode("Arial Italic 6"));
                } else {
                    p.addBackgroundText(d.toString(0) + "\u03bcm",
                                        RIGHT, 0, Font.decode("Arial Italic 6"));
                }

                p.setScaleBar(500 / 5.8, 1, 0.05, 0.1);
//                p.setDrawBoxLines(false, false, false, false);
                p.saveAsEPS(new File("fig2-" + n++ +".eps"), w, h);
                p.setLegendVisible(true);
            }
        }

        // save time courses
        NeuronViewer nv = new NeuronViewer(pName[0], false, new Config("config.xml"));
        n = saveTCPlot(nv, n, -0.08, 0.028, .02, 'f');

        // save flashes
        ScatterPlot sp = new ScatterPlot();
        int a = 100;
        sp.add(0, a + 0);
        sp.add(0, a + 10);
        sp.add(2, a + 10);
        sp.add(2, a + 0);
        sp.add(4, a + 0);
        sp.add(4, a + -10);
        sp.add(6, a + -10);
        sp.add(6, a + 0);
        sp.add(8, a + 0);
        sp.add(8, a + 10);

        ScatterPlot sp1 = new ScatterPlot();
        sp1.add(1, 0);
        sp1.add(1, 100);

        for (int c = 1; c < cName.length; c++) {
            ParametersFile pFile = nv.getParametersFile();
            double[] average = pFile.getArrayAverage("flashResponse", cName[c]);
            MathUtil.rotateArray(average, -1);
            double[] x = new double[average.length];
            for (int i = 0; i < x.length; i++) {
                x[i] = i * 8.0 / 201.0;
            }

            PlotPanel p = new PlotPanel();
            p.addData(new ScatterPlot(x, average, null),
                      new ScatterPlotStyle("NONE 0 black, SOLID 0.5 black"));
            p.addData(sp, new ScatterPlotStyle("NONE 0 black, SOLID 0.5 black"));
            p.addData(sp1, new ScatterPlotStyle("NONE 0 black, SOLID 0.5 black"));

            p.setLabelFont(labelFont);
            p.setRange(0, 8, 0, 300);
            p.setLabels("Time (s)", "Response (spikes/s)");
            p.getXAxis().shiftRightLabel = true;
            if (c == 1) {
                p.addBackgroundText("" + 'e', LEFT, TOP, letterFont);
            }

            p.getYAxis().setFixedTickSpacing(0, 100, 2);
            p.getYAxis().setMarkToLabelSpacing(0);
            p.getYAxis().setTickToMarkSpacing(1);

            p.getXAxis().setMarkToLabelSpacing(0);
            p.getXAxis().setTickToMarkSpacing(2);

            p.axesBorder.setTopPadding(3);

            p.saveAsEPS(new File("fig2-" + n++ +".eps"), w, 1.3 * h, false);
        }

        n = saveTCPlot(new NeuronViewer("F:\\Good Data\\2005-09-09-1\\data007", false, new Config("config.xml")),
                       n, -0.125, 0.05, .04, 'g');
    }


    public static Num getNonlinearityIndex(int[] id, ParametersFile pFile, int sPeriod) throws
        IOException {

        MeanVarianceCalculator mvc = new MeanVarianceCalculator();
        for (int index = 0; index < id.length; index++) {
            double[] _f1 = pFile.getArrayCell(id[index], "T1reversingF1");
            double[] _f2 = pFile.getArrayCell(id[index], "T1reversingF2");
            if (_f1 != null && _f2 != null) {
                mvc.add(_f2[sPeriod] / _f1[sPeriod]);
            }
        }
        return new Num(mvc.getMean(), mvc.getMeanVariance());
    }


    public static void makeFigure1(String[] pName) throws Exception {
        NeuronViewer nv = new NeuronViewer(pName[0], false, new Config("config.xml"));

        // save clustering plots
        String TF1 = "<html>Time Filter PC<sub>1</sub> (a.u.)";
        String TF2 = "<html>Time Filter PC<sub>2</sub> (a.u.)";
        int thickness = 20; //11;
        int s1 = 0;
        int s2 = 1;

        nv.loadClassification("f:\\Good Data\\2005-04-26-0\\data009\\c1.txt");
        nv.select("All");
        PlotPanel p = nv.makeHistogram(
            nv, "extreme(RedTimeCourse + GreenTimeCourse + BlueTimeCourse)",
            0.01, nv.rightTree.getSelectionPath());
        p.setLabels("Time Filter Amplitude (a.u.)", "Number of cells");
        p.autoscale();
        p.setYRange(0, 39);
        p.setLegendVisible(false);
        p.setLabelFont(labelFont);
        p.addBackgroundText("" + 'a', RIGHT, TOP, letterFont);
        p.getXAxis().setFixedTickSpacing(0, .2, 1);
        p.getXAxis().multiplier = 5;
        p.setXRange( -2.1 / 5, 2.1 / 5);
        p.getXAxis().setSpacings(s1, 2 * s2);
        p.getYAxis().setSpacings(s1, s2);
        p.addBackgroundText("ON", 80, CENTER, textFont);
        p.addBackgroundText("OFF", LEFT, CENTER, textFont);
        p.getXAxis().setThickness(thickness);
        p.saveAsEPS(new File("fig1-0.eps"), w, 0.75 * w, false);

        // pca plot
        nv.loadClassification("f:\\Good Data\\2005-04-26-0\\data009\\c2.txt");
        nv.select("All/OFF");

        String x1 =
            "pca(norm(RedTimeCourse[30, 59]#GreenTimeCourse[30, 59]#BlueTimeCourse[30, 59]), 0)";
        String x2 =
            "pca(norm(RedTimeCourse[30, 59]#GreenTimeCourse[30, 59]#BlueTimeCourse[30, 59]), 1)";
        p = nv.makeScatterPlot(
            x1, x2, TF1, TF2, nv.rightTree.getSelectionPath());
        ScatterPlotStyle s = (ScatterPlotStyle) p.getStyleWithName("Scatter Plot");
        s.setSymbolSize(2);
        p.setLabelFont(labelFont);
        p.addBackgroundText("" + 'b', RIGHT, TOP, letterFont);
        p.getXAxis().setSpacings(s1, 2 * s2);
        p.getYAxis().setSpacings(s1, s2);

        p.getXAxis().setFixedTickSpacing(0, .1, 1);
        p.getXAxis().multiplier = 10;
        p.getYAxis().setFixedTickSpacing(.1, .05, 2);
        p.getYAxis().multiplier = 20;

        p.getXAxis().setThickness(thickness);
        p.addBackgroundText("Midget", LEFT, TOP, textFont);
        p.addBackgroundText("Parasol", RIGHT, 35, textFont);
        p.addBackgroundText("<html>C<sub>1", CENTER, 50, textFont);
        p.addBackgroundText("<html>C<sub>2", LEFT, 37, textFont);
        p.saveAsEPS(new File("fig1-1.eps"), w, 0.75 * w, false);

        nv.loadClassification("f:\\Good Data\\2005-04-26-0\\data009\\c3.txt");
        nv.select("All/OFF/C2");
        p = NeuronViewer.makeHistogram(nv, sizeExp, 20, nv.rightTree.getSelectionPath());
        p.setLabels("RF Diameter (\u03bcm)", "Number of cells");
        p.setRange(0, 599, 0, 4.5);
        p.setLegendVisible(false);
        p.setLabelFont(labelFont);
        p.addBackgroundText("" + 'c', RIGHT, TOP, letterFont);
        p.getXAxis().setSpacings(s1, 2 * s2);
        p.getYAxis().setSpacings(s1, s2);
        p.getXAxis().setThickness(thickness);
        p.addBackgroundText("<html>C<sub>1B", LEFT, 20, textFont);
        p.addBackgroundText("<html>C<sub>1A", 73, 20, textFont);
        p.saveAsEPS(new File("fig1-2.eps"), w, 0.75 * w, false);
    }


    public static void makeSIFigure1(String[] pName) throws IOException {
        NeuronViewer nv = null;
        try {
            nv = new NeuronViewer(pName[0], false, new Config("config.xml"));
        } catch(Exception e) {
            e.printStackTrace();
        }

        // save example EIs
        nv.loadClassification("f:\\Good Data\\2005-04-26-0\\data009\\c1.txt");
        int n = 3;
        String[] neurons = {"All/1280", "All/1670" /*, "All/158"*/};
        int[][] el = { {336, 327, 89, 169}, {433, 361, 94, 149}, {45, 277, 410, 189}
        };
        int[] range = {199, 199, 350};
        double[] angle = {90, -135, -135};
        int[][] f = { {1, 5, 2, 1}, {1, 5, 2, 2}, {1, 20, 30, 30}
        };
        String[][] text = { {"Cell Body", "Dendrite", "Axon 1", "Axon 2"}, {"Cell Body",
                          "Dendrite", "Axon 1", "Axon 2"}, {"Cell Body", "Axon 1",
                          "Axon 2", "Axon 3"},
        };
        nv.getPlotMaker("EI").isSelected = true;
        for (int i = 0; i < neurons.length; i++) {
            nv.select(neurons[i]);

            PhysiologicalImagePanel p1 = (PhysiologicalImagePanel) nv.getPlot("ei");
            p1.setAmplitudeFactor( (i == 2) ? 10 : 1);
            p1.markerElectrode = -1;
//            p1.addBackgroundText("" + (char) ('a' + n), RIGHT, TOP, letterFont);
//            p1.setScaleBar(200, 1, 0.05, 0.9);
            p1.removeAll();
//            GraphicsIO.saveAsEPS(p1, "fig1-" + n++ +".eps", w * 1.1, h * 1.1);
//            p1.removeAllBackgroundText();

            Font ff = Font.decode("Arial Italic 8");
            p1.setScaleBar(200, 2, 0.05, 0.9);
            if (i == 0) {
                p1.addBackgroundText("Cell Body", 115, 15, ff);
                p1.addBackgroundText("Dendrite", 95, 0, ff);
                p1.addBackgroundText("Axon 1", 75, 25, ff);
                p1.addBackgroundText("Axon 2", LEFT, 37, ff);
            }
            if (i == 1) {
                p1.addBackgroundText("Cell Body", 125, 56, ff);
                p1.addBackgroundText("Dendrite", 93, 23, ff);
                p1.addBackgroundText("Axon 1", 70, 85, ff);
                p1.addBackgroundText("Axon 2", 28, BOTTOM, ff);
            }
            if (i == 2) {
                p1.addBackgroundText("Cell Body", 81, 75, ff);
                p1.addBackgroundText("Axon 1", 45, 20, ff);
                p1.addBackgroundText("Axon 2", 155, 15, ff);
                p1.addBackgroundText("Axon 3", 5, 70, ff);
            }
            for (int k = 0; k < el[i].length; k++) {
                //Deprecated implementation.  Must be redone if needed
 //               p1.addArrow(el[i][k], new Arrow( -1, -1, 2, angle[i], 15, 15, 7, 1f));
            }
            GraphicsIO.saveAsEPS(p1, "fig1-SI-" + (i + 1) + ".eps", 2.9, 2.9 / 2, false);

            // save SI plots
            for (int k = 0; k < el[i].length; k++) {
                PlotPanel pan = new PlotPanel();
                pan.setDoubleBuffered(false);
                pan.addData(p1.getPlotFor(el[i][k]), p1.style);

                pan.getYAxis().setFixedTickSpacing(0, 100 / f[i][k], 0);
                pan.getYAxis().setShowlabel(false);
//                pan.getYAxis().setShowTicks(true);
                pan.getYAxis().setThickness(30);
                ( (AxesBorder) pan.getBorder()).setPadding(0, 0, 0, 0);

                pan.addBackgroundText(text[i][k], RIGHT, TOP, textFont);
//                if (f[i][k] > 1) {
//                    pan.addBackgroundText("x" + f[i][k], LEFT, BOTTOM, textFont);
//                }
                if (k != 3) {
                    pan.getXAxis().setShowlabel(false);
                    pan.getXAxis().setShowTicks(false);
                }
                pan.getXAxis().setMarkToLabelSpacing(0);

                pan.setLabels("Time (ms)", "Amplitude (\u03BCV)");
                pan.setLabelFont(Font.decode("Arial Italic 9"));
                if (i == 2) {
                    pan.setRange( -0.9, 2.1, -range[i] / f[i][k], range[i] / f[i][k]);
                } else {
                    pan.setRange( -0.9, 2.1, -range[i] / f[i][k], range[i] / f[i][k]);
                }
                GraphicsIO.saveAsEPS(pan, "fig1-SI-" + (i + 1) + "-" + (k + 1) + ".eps",
                                     2.1, 0.8, false);
            }
        }
    }


    public static void makeFigure4(String[] pName, String[] crgName) throws Exception {
        double[] max = {75, 63, 4.5};
        int[] n = {2, 1, 1};
        int thick = -1;

        for (int c = 0; c < 3; c++) {
            DoubleHistogram nHist = new DoubleHistogram("", 0, 15, 0.5);
            MeanVarianceCalculator mvc = new MeanVarianceCalculator();

            for (int pIndex = 0; pIndex < crgName.length; pIndex++) {
                ReversingGratings m = new ReversingGratings(crgName[pIndex],
                    binsPerPeriod);
                ParametersFile pFile = new ParametersFile(pName[pIndex]);

                // make n-index hist
                int[] id = pFile.getNeuronsInClass(cName[c]);
                for (int k = 0; k < id.length; k++) {
                    double[] f1 = pFile.getArrayCell(id[k], "T1reversingF1");
                    double[] f2 = pFile.getArrayCell(id[k], "T1reversingF2");
                    double maxN = Double.NEGATIVE_INFINITY;

                    for (int i = 0; i < f1.length; i++) {
                        if (f1[i] > 0 || f2[i] > 0) {
                            if (f2[i] / f1[i] > maxN) {
                                maxN = f2[i] / f1[i];
                            }
                        }
                    }

                    nHist.fill(maxN, 1);
                    mvc.add(maxN);
                }
            }

            HistogramStyle h = new HistogramStyle(
                "", HistogramStyle.OutlineType.RECTANGULAR, black, 0.5f, false,
                black, false, black, 1);
            PlotPanel p = PlotUtil.showData("", nHist, h, 500, 300);

            p.setXRange(0, 15);
            p.setLabelFont(labelFont);
            p.setLabels("Nonlinearity Index", "Number of cells");
            p.getXAxis().setShowlabel(c == 2);
//            p.getXAxis().setShowTickMarks(c == 2);
            p.getYAxis().setThickness(35);

            p.addBackgroundText("" + (char) ('a' + c), 3, -3, letterFont);
            p.addBackgroundText(className[c], 200, 7, textFont);

            p.setYRange(0, max[c]);
            //<html>F<sub>2</sub> /F<sub>1</sub>
//            p.addBackgroundText(
//                "mean +/- s.e.m. = " + StringUtil.format(mvc.getMean(), n[c], 1) + " \u00B1 " +
//                StringUtil.format(mvc.getMeanVariance(), n[c], 1), x, 13, textFont);

            if (c == 0) {
                thick = p.getYAxis().getThickness();
            }
            p.getYAxis().setThickness(thick);
            p.saveAsEPS(new File("fig4-" + c + ".eps"), 2 * w, 0.9, false);
        }
    }


    public static void calculateParameters(String[] pName) throws Exception {
        ParametersFile[] pf = new ParametersFile[pName.length];
        for (int prepIndex = 0; prepIndex < pName.length; prepIndex++) {
            pf[prepIndex] = new ParametersFile(pName[prepIndex]);
        }

        // Latency
        Num[] ypLatency = new Num[pName.length];
        Num[] ymLatency = new Num[pName.length];
        for (int prepIndex = 0; prepIndex < pName.length; prepIndex++) {
            Num mLatency = pf[prepIndex].getAverage1("-rl", "All/OFF/Midget");
            Num pLatency = pf[prepIndex].getAverage1("-rl", "All/OFF/Parasol");
            Num yLatency = pf[prepIndex].getAverage1("-rl", "All/OFF/BTY");
            ypLatency[prepIndex] = yLatency.div(pLatency);
            ymLatency[prepIndex] = yLatency.div(mLatency);

            System.err.println("\nPreparation " + (prepIndex + 1));
            System.err.println("   Midget latency: " + mLatency.toString(1));
            System.err.println("   Parasol latency: " + pLatency.toString(1));
            System.err.println("   BTY latency: " + yLatency.toString(1));
            System.err.println("   BTY/Parasol Latency: " +
                               ypLatency[prepIndex].toString(2));
            System.err.println("   BTY/Midget latency: " +
                               ymLatency[prepIndex].toString(2));
        }

//        System.err.println("\nAverage values");
//        System.err.println("BTY/Midget  Latency: " + Num.average(ymLatency).toString(2));
//        System.err.println("BTY/Parasol Latency: " + Num.average(ypLatency).toString(2));

        System.err.println("\nRMS values");
        System.err.println("BTY/Midget  Latency: " + Num.rms(ymLatency).toString(2));
        System.err.println("BTY/Parasol Latency: " + Num.rms(ypLatency).toString(2));
        /*
                // DOT
                System.err.println("\n\nDOT:");
                Num[] dot = new Num[pName.length];
                for (int cIndex = 0; cIndex < cName.length; cIndex++) {
                    for (int prepIndex = 0; prepIndex < pName.length; prepIndex++) {
//                System.err.println(pName[prepIndex] + ": " + cName[cIndex]);
                        dot[prepIndex] = pf[prepIndex].getAverage1("dot", cName[cIndex]);
                    }

//            System.err.println(cName[cIndex] + ": " + Num.average(dot).toString(2));
                    System.err.println(cName[cIndex] + ": " + Num.rms(dot).toString(2));
                }
         */
        // RF Diam. ratio
        System.err.println("\n\nRF Diameter:");
        Num[] diamRatio = new Num[pName.length];
        Num[] diamRatio1 = new Num[pName.length];
        for (int prepIndex = 0; prepIndex < pName.length; prepIndex++) {
//            Num onP = pf[prepIndex].getAverage1(
//                "2*5.8*20*((SigmaX*SigmaY)^0.5)", "All/ON/Parasol");
            Num offP = pf[prepIndex].getAverage1(
                "2*5.8*20*((SigmaX*SigmaY)^0.5)", "All/OFF/Parasol");
            Num P = pf[prepIndex].getAverage1(
                "2*5.8*20*((SigmaX*SigmaY)^0.5)", "All/OFF/Parasol", "All/ON/Parasol");
            Num offBTY = pf[prepIndex].getAverage1(
                "2*5.8*20*((SigmaX*SigmaY)^0.5)", "All/OFF/BTY");
            diamRatio[prepIndex] = offBTY.div(P);
            diamRatio1[prepIndex] = offBTY.div(offP);

//            System.err.println("ON " + onP.toString(0));
//            System.err.println("OFF " + offP.toString(0));
//            System.err.println("ON + OFF " + P.toString(0));
            System.err.println("BTY/OFF P: " + diamRatio1[prepIndex].toString(2));
//            System.err.println("BTY/ON&OFF P: " + diamRatio[prepIndex].toString(2));
        }

        System.err.println("RMS BTY/OFF P: " + Num.rms(diamRatio1).toString(1));
//        System.err.println("RMS BTY/ON&OFF P: " + Num.rms(diamRatio).toString(1));

        /*
                 // Blue - ratio
                 //           String b = "abs(extreme(BlueTimeCourse))/(abs(extreme(GreenTimeCourse)) + abs(extreme(RedTimeCourse)) + abs(extreme(BlueTimeCourse)))";
//           String b = "extreme(BlueTimeCourse)/(extreme(GreenTimeCourse) + extreme(RedTimeCourse) + extreme(BlueTimeCourse))";
//           String b = "extreme(BlueTimeCourse)/(extreme(GreenTimeCourse) + extreme(RedTimeCourse))";
                 String b = "extreme(BlueTimeCourse) / extreme(RedTimeCourse)";

                 Num P = Num.zero;
                 Num BTY = Num.zero;
                 for (int prepIndex = 0; prepIndex < pName.length; prepIndex++) {
                     Num _P = pf[prepIndex].getAverage1(b, "All/OFF/Parasol");
                     Num _BTY = pf[prepIndex].getAverage1(b, "All/OFF/BTY");

//           System.err.println(pName[prepIndex]);
//           System.err.println(_P.toString(3));
//           System.err.println(_BTY.toString(3));

                     P = P.add(_P);
                     BTY = BTY.add(_BTY);
                 }
                 P = P.div(3);
                 BTY = BTY.div(pName.length);

                 System.err.println("Average");
                 System.err.println(P.toString(3));
                 System.err.println(BTY.toString(3));

                 System.err.println("Ratio: " + BTY.div(P).toString(2));
         */
    }


    public static ScatterPlot[] makePlot(int[] id, ParametersFile pFile, double f) throws
        IOException {

        double spread = 0;

        double[] freqs = pFile.getArrayCell(id[0], "reversingFrequencies");
        final int nFreq = freqs.length;
        MeanVarianceCalculator mvc = new MeanVarianceCalculator();

        // calculate overall scaling factor
        ScatterPlot sp = new ScatterPlot();
        ScatterPlot sp1 = new ScatterPlot();
        for (int sPeriod = 0; sPeriod < nFreq; sPeriod++) {
            mvc.reset();
            for (int index = 0; index < id.length; index++) {
                double[] _f1 = pFile.getArrayCell(id[index], "T1reversingF1");
                double[] _f2 = pFile.getArrayCell(id[index], "T1reversingF2");
                if (_f1 != null) {
                    double dx = random() - 0.5;
                    double x = freqs[sPeriod] * (1 + spread * dx) * f;
                    sp.add(x, _f2[sPeriod] / _f1[sPeriod]);
                    mvc.add(_f2[sPeriod] / _f1[sPeriod]);
                }
            }
            sp1.add(freqs[sPeriod], mvc.getMean(), mvc.getMeanVariance());
        }

        return new ScatterPlot[] {sp, sp1};
    }


    public static void main(String[] args) throws Exception {
        final String[] pName = {
                               "f:\\Good Data\\2005-04-26-0\\data009\\data009.params",
                               "f:\\Good Data\\2005-09-09-1\\data002\\data002.params",
                               "f:\\Good Data\\2005-04-14-0\\data002\\data002.params",
        };
        final String[] crgName = {
                                 "f:\\Good Data\\2005-04-26-0\\data005-mapped-data009",
                                 "f:\\Good Data\\2005-09-09-1\\data003",
                                 "f:\\Good Data\\2005-04-14-0\\data004",
        };

        SwingUtilities.invokeAndWait(new Runnable() {
            public void run() {
                try {
                    
                    makeFigure1(pName);
                    makeFigure2(pName);
                    makeFigure3(pName, crgName);
                    makeFigure4(pName, crgName);

                    makeSIFigure1(pName);

//                    calculateParameters(pName);
                } catch (Exception ex) {
                    ex.printStackTrace();
                }
            }
        });
        System.exit(1);
    }
}
