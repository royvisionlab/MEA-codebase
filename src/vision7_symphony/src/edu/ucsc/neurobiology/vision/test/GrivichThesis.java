package edu.ucsc.neurobiology.vision.test;

import java.io.*;
import java.util.*;

import java.awt.*;
import javax.swing.*;


import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.neuronviewer.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.util.*;
import edu.ucsc.neurobiology.vision.plot.eps.EpsGraphics2D;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreePath;
import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import edu.ucsc.neurobiology.vision.math.fitting.Gaussian1DFunction;
import edu.ucsc.neurobiology.vision.analysis.OMSClassification;
import edu.ucsc.neurobiology.vision.gui.PlaceLayout;
import edu.ucsc.neurobiology.vision.gui.PlaceC;
import edu.ucsc.neurobiology.vision.gui.GraphicsIO;
import edu.ucsc.neurobiology.vision.analysis.ReversingGratings;
import javax.swing.border.Border;
import static edu.ucsc.neurobiology.vision.plot.PlotPanel.*;
import edu.ucsc.neurobiology.vision.*;


/**
 * @author Matthew Grivich, University of California, Santa Cruz

 */
public class GrivichThesis {


    public static final int THESIS = 0;
    public static final int PAPER = 1;
    public static final int POWERPOINT = 2;

    public static int mode = POWERPOINT;



    static Font letterFont = Font.decode("Arial Bold 10");
    static Font textFont = Font.decode("Arial Italic 8");
    static Font labelFont = Font.decode("Arial Italic 8");

    static double[] referenceSimplifiedContourArea;
    static double[] referenceContourArea;
    static double[] referenceRL;
    static double[] referenceAxon;

    static NeuronViewer[] nvs, whitenedNVS;
    static ParametersFile[] parameterFiles, whitenedParameterFiles;
    static String[] whitenedParamFilesString = {
                                               "G:\\data\\2005-05-02-0\\data002W\\data002W.params",
                                               "G:\\data\\2005-05-04-0\\data002W\\data002W.params",
                                               "G:\\data\\2005-08-25-0\\data003W\\data003W.params"};

    static String[] paramFilesString = {"G:\\data\\2005-05-02-0\\data002\\data002.params",
                                       "G:\\data\\2005-05-04-0\\data002\\data002.params",
                                       "G:\\data\\2005-08-25-0\\data003\\data003.params"};


    static String[] whitenerParamFilesString = {
                                               "G:\\data\\2005-08-25-0\\data003\\data003.params",
                                               "G:\\data\\2005-08-25-0\\data003W\\data003W.params",
                                               "G:\\data\\2005-08-25-0\\data002\\data002.params"};


    public static String[] classList = {
                                       "All/On/Small/T2/OMS",
                                       "All/On/Small/T3/OMS",
                                       "All/On/Large/T1/Triphasic",
                                       "All/On/Large/T2/DS",
                                       "All/On/Large/T2/AOMS",
                                       "All/On/Large/T2/Composite",
                                       "All/On/Large/T2/A",
                                       "All/On/Large/T3/UniformityA",
                                       "All/On/Large/T3/UniformityB",
                                       "All/On/Large/T3/Blue",
                                       "All/On/Large/T3/OS",
                                       "All/On/Large/T3/Triphasic",
                                       "All/On/Large/T3/A",
                                       "All/On/Large/T3/B",
                                       "All/On/Large/T3/C",
                                       "All/On/Large/T3/D",
                                       "All/On/Large/T3/E",
                                       "All/On/Large/T3/F",
                                       "All/Off/Small/T2/DS",
                                       "All/Off/Small/T2/OMS",
                                       "All/Off/Small/T3/OMS",
                                       "All/Off/Small/T3/A",
                                       "All/Off/Small/T3/B",
                                       "All/Off/Small/T3/C",
                                       "All/Off/Large/T1/AOMS",
                                       "All/Off/Large/T1/A",
                                       "All/Off/Large/T3/A"
    };

    private static String referenceType = "All/Off/Small/T2/OMS";
//





    public static void makeResultsFigures() throws Exception {
        double dpi = 72.0;
        int width = (int) (4.0 * (dpi));
        int height = (int) (2.0 * (dpi));

        int verticalAxisThickness = 32;
        int rightPadding = 8;
        int topPadding = 4;
        Font labelFont = Font.decode("Arial ITALIC 8");

        ScatterPlotStyle scatterPlotStyle = new ScatterPlotStyle(
            "Scatter Plot", SymbolType.DISK, 3, Color.black, false, Color.black, 1);

//Make On DS plot
        PlotPanel p = new PlotPanel("", false, false, false, false);

        ScatterPlot sp = new ScatterPlot();

        DefaultMutableTreeNode node = nvs[0].select("All/On/Large/T2/DS");
        IntegerList list = InteractiveTree.getNeuronsInClass(node,
            new IntegerList(), true);
        int[] neuronsInClass = list.toArray();
        HashMap<Integer,
            Double> xvals = parameterFiles[0].evaluate("xOnDS", neuronsInClass);
        HashMap<Integer,
            Double> yvals = parameterFiles[0].evaluate("yOnDS", neuronsInClass);

        for (int i = 0; i < neuronsInClass.length; i++) {
            double x = ( (Double) xvals.get(new Integer(neuronsInClass[i]))).doubleValue();
            double y = ( (Double) yvals.get(new Integer(neuronsInClass[i]))).doubleValue();
            double theta = 1.0;
            double xNew = -x * Math.cos(theta) - y * Math.sin(theta);
            double yNew = -x * Math.sin(theta) + y * Math.cos(theta);
            sp.add(xNew, yNew);

        }
        p.addData(sp, scatterPlotStyle);
//         p.addArrow(new Arrow(0, 0, 1, .2, 5, .2, 2, 2f));

//          p.addArrow(new Arrow( 100, 100, 0, 270, 50, 15, 20, 2f));
//          p.addArrow(new Arrow( 50, 50, 0, 0, 50, 15, 20, 2f));

        p.setRange( -.8, .8, -.8, .8);
        p.setAntiallias(true);
        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);

        p.setLabelFont(labelFont);
        p.addBackgroundText("CP5", 40, 30, letterFont, Color.black);
        p.addBackgroundText("CP4", 152, 56, letterFont, Color.black);
        p.addBackgroundText("CP6", 100, 121, letterFont, Color.black);
        p.getXAxis().setLabel("Horizontal Direction Sensitivity (+ is Anterior)");
        p.getYAxis().setLabel("Vertical Direction Sensitivity (+ is Superior)");
//        p.saveAsEPS(new File(
//            "g:\\Publications\\thesis\\figures\\results\\OnDS.eps"), 2, 2, true);

        JPanel panel = new JPanel(new PlaceLayout());
        panel.add(p, new PlaceC(0, 0, 0, 0, .5, 1));

        whitenedNVS[0].showNeuronsInSubclasses = true;
        whitenedNVS[0].select("All/On/Large/T2/DS/3");
        p = (PlotPanel) whitenedNVS[0].getPlot("mosaic");
        p.setDoubleBuffered(false);
        p.axesBorder.setRightPadding(2);
        p.axesBorder.setBottomPadding(2);
        p.setLegendVisible(false);
        p.setScaleBar(500, 1, 0.05, 0.1);
        p.addBackgroundText("CP5", LEFT, TOP, letterFont, Color.black);

        panel.add(p, new PlaceC(.5, .20, 0, 0, .25, .30));

        whitenedNVS[0].select("All/On/Large/T2/DS/2");
        p = (PlotPanel) whitenedNVS[0].getPlot("mosaic");
        p.setDoubleBuffered(false);
        p.axesBorder.setRightPadding(2);
        p.axesBorder.setBottomPadding(2);
        p.setLegendVisible(false);
        p.addBackgroundText("CP4", LEFT, TOP, letterFont, Color.black);
        panel.add(p, new PlaceC(.75, .20, 0, 0, .25, .30));

        whitenedNVS[0].select("All/On/Large/T2/DS/1");
        p = (PlotPanel) whitenedNVS[0].getPlot("mosaic");
        p.setDoubleBuffered(false);
        p.axesBorder.setRightPadding(2);
        p.axesBorder.setBottomPadding(2);
        p.setLegendVisible(false);
        p.addBackgroundText("CP6", LEFT, TOP, letterFont, Color.black);
        panel.add(p, new PlaceC(.625, .504, 0, 0, .25, .30));

        JFrame frame = new JFrame();
        frame.add(panel);
        frame.setSize(6 * 72, 3 * 72);
        frame.setVisible(true);

        RepaintManager.currentManager(frame).setDoubleBufferingEnabled(false);
        EpsGraphics2D g = new EpsGraphics2D(
            "", Font.decode("Arial ITALIC 8"), 0, 0, frame.getWidth(), frame.getHeight());

        g.scale(72.0 / dpi, 72.0 / dpi);
        g.setAccurateTextMode(true);
        frame.paint(g);

        RepaintManager.currentManager(frame).setDoubleBufferingEnabled(true);

        FileWriter writer = new FileWriter(
            "g:\\Publications\\thesis\\figures\\results\\OnDS.eps");

        writer.write(g.toString());
        writer.flush();
        writer.close();

        //Off ds plot
        p = new PlotPanel("", false, false, false, false);

        sp = new ScatterPlot();

        node = nvs[0].select("All/Off/Small/T2/DS");
        list = InteractiveTree.getNeuronsInClass(node,
                                                 new IntegerList(), true);
        neuronsInClass = list.toArray();
        xvals = parameterFiles[0].evaluate("xOffDS", neuronsInClass);
        yvals = parameterFiles[0].evaluate("yOffDS", neuronsInClass);

        for (int i = 0; i < neuronsInClass.length; i++) {
            double x = ( (Double) xvals.get(new Integer(neuronsInClass[i]))).doubleValue();
            double y = ( (Double) yvals.get(new Integer(neuronsInClass[i]))).doubleValue();
            double theta = 1.0;
            double xNew = -x * Math.cos(theta) - y * Math.sin(theta);
            double yNew = -x * Math.sin(theta) + y * Math.cos(theta);
            sp.add(xNew, yNew);

        }
        p.addData(sp, scatterPlotStyle);

        p.setRange( -.7, .7, -.7, .7);
        p.setAntiallias(true);
        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);

        p.setLabelFont(labelFont);
        p.addBackgroundText("CP23", 40, 50, letterFont, Color.black);
        p.addBackgroundText("CP21", 152, 68, letterFont, Color.black);
        p.addBackgroundText("CP22", 100, 20, letterFont, Color.black);
        p.addBackgroundText("CP24", 100, 140, letterFont, Color.black);

        p.getXAxis().setLabel("Horizontal Direction Sensitivity (+ is Anterior)");
        p.getYAxis().setLabel("Vertical Direction Sensitivity (+ is Superior)");

        panel = new JPanel(new PlaceLayout());
        panel.add(p, new PlaceC(0, 0, 0, 0, .5, 1));

        whitenedNVS[0].select("All/Off/Small/T2/DS/2");
        p = (PlotPanel) whitenedNVS[0].getPlot("mosaic");
        p.setDoubleBuffered(false);
        p.axesBorder.setRightPadding(2);
        p.axesBorder.setBottomPadding(2);
        p.setLegendVisible(false);
        p.setScaleBar(500, 1, 0.05, 0.1);
        p.addBackgroundText("CP22", LEFT, TOP, letterFont, Color.black);

        panel.add(p, new PlaceC(.5, .20, 0, 0, .25, .30));

        whitenedNVS[0].select("All/Off/Small/T2/DS/4");
        p = (PlotPanel) whitenedNVS[0].getPlot("mosaic");
        p.setDoubleBuffered(false);
        p.axesBorder.setRightPadding(2);
        p.axesBorder.setBottomPadding(2);
        p.setLegendVisible(false);
        p.addBackgroundText("CP21", LEFT, TOP, letterFont, Color.black);
        panel.add(p, new PlaceC(.75, .20, 0, 0, .25, .30));

        whitenedNVS[0].select("All/Off/Small/T2/DS/1");
        p = (PlotPanel) whitenedNVS[0].getPlot("mosaic");
        p.setDoubleBuffered(false);
        p.axesBorder.setRightPadding(2);
        p.axesBorder.setBottomPadding(2);
        p.setLegendVisible(false);
        p.addBackgroundText("CP23", LEFT, TOP, letterFont, Color.black);
        panel.add(p, new PlaceC(.5, .504, 0, 0, .25, .30));

        whitenedNVS[0].select("All/Off/Small/T2/DS/3");
        p = (PlotPanel) whitenedNVS[0].getPlot("mosaic");
        p.setDoubleBuffered(false);
        p.axesBorder.setRightPadding(2);
        p.axesBorder.setBottomPadding(2);
        p.setLegendVisible(false);
        p.addBackgroundText("CP24", LEFT, TOP, letterFont, Color.black);
        panel.add(p, new PlaceC(.75, .504, 0, 0, .25, .30));

        frame = new JFrame();
        frame.add(panel);
        frame.setSize(6 * 72, 3 * 72);
        frame.setVisible(true);

        RepaintManager.currentManager(frame).setDoubleBufferingEnabled(false);
        g = new EpsGraphics2D(
            "", Font.decode("Arial ITALIC 8"), 0, 0, frame.getWidth(), frame.getHeight());

        g.scale(72.0 / dpi, 72.0 / dpi);
        g.setAccurateTextMode(true);
        frame.paint(g);

        RepaintManager.currentManager(frame).setDoubleBufferingEnabled(true);

        writer = new FileWriter(
            "g:\\Publications\\thesis\\figures\\results\\OffDS.eps");

        writer.write(g.toString());
        writer.flush();
        writer.close();

//Make OS plot

        nvs[1].showNeuronsInSubclasses = false;
        nvs[1].select("/All/On/Large/T3/OS");
        p = (PlotPanel) nvs[1].getPlot("overlayedSinusoids");
        p.setLabelFont(labelFont);
        p.setRange(p.getRange()[0], p.getRange()[1], 0, p.getRange()[3]);
        p.setDoubleBuffered(false);
        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);
        PlotUtil.showData("", p);
        p.saveAsEPS(new File(
            "g:\\Publications\\thesis\\figures\\results\\OS.eps"), 2, 2, true);

        //Make BluenessPlot
        p = new PlotPanel();

        node = nvs[1].select("All/On/Large/T3/UniformityA");
        list = InteractiveTree.getNeuronsInClass(node,
                                                 new IntegerList(), true);
        neuronsInClass = list.toArray();

        HashMap<Integer,
            Double> vals = parameterFiles[1].evaluate("blueness", neuronsInClass);
        DoubleHistogram h = new DoubleHistogram("", .8, 1.5, .025);
        HistogramStyle hs = new HistogramStyle(Color.green, 1);

        for (int i = 0; i < neuronsInClass.length; i++) {
            h.fill( ( (Double) vals.get(new Integer(neuronsInClass[i]))).doubleValue(),
                   1.0);

        }

        p.addData(h, hs);

        node = nvs[1].select("All/On/Large/T3/Blue");
        list = InteractiveTree.getNeuronsInClass(node,
                                                 new IntegerList(), true);
        neuronsInClass = list.toArray();

        vals = parameterFiles[1].evaluate("blueness", neuronsInClass);
        h = new DoubleHistogram("", .8, 1.5, .025);
        hs = new HistogramStyle(Color.blue, 1);

        for (int i = 0; i < neuronsInClass.length; i++) {
            h.fill( ( (Double) vals.get(new Integer(neuronsInClass[i]))).doubleValue(),
                   1.0);

        }

        p.addData(h, hs);
        p.setLabels("Sensitivity to Blue", "Number of Neurons");
        p.setRange(.8, 1.5, 0, 11);
        p.setLabelFont(labelFont);
        p.setDoubleBuffered(false);
        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);

        PlotUtil.showData("", p);
        p.saveAsEPS(new File(
            "g:\\Publications\\thesis\\figures\\results\\blueness.eps"), 2, 2, true);
        //REMEMBER TO SET BIN SIZE TO 10 IN CONFIG.XML.
        node = nvs[0].select("All/Off/Large/T1/A");
        p = (PlotPanel) nvs[0].getPlot("spike rate");

        p.setLabelFont(labelFont);
        p.setDoubleBuffered(false);
        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);

        PlotUtil.showData("", p);
        p.saveAsEPS(new File(
            "g:\\Publications\\thesis\\figures\\results\\oscillatingRetina.eps"), 4, 2, true);

    }


    public static void makeMethodsFigures() throws Exception {
        NeuronViewer nv = new NeuronViewer(
            "G:\\data\\2005-05-02-0\\data002\\data002.params", false, new Config("config.xml"));

        JFrame frame = new JFrame();

        frame.setLayout(new GridBagLayout());

        frame.setLayout(new GridBagLayout());
        GridBagConstraints c = new GridBagConstraints();

        double dpi = 72.0;
        int width = (int) (4.0 * (dpi));
        int height = (int) (2.0 * (dpi));

        int verticalAxisThickness = 32;
        int rightPadding = 8;
        int topPadding = 4;
        Font labelFont = Font.decode("Arial ITALIC 8");
        frame.setSize(width + 8, (int) (height * 1.75 + 34));

        /*STA Example
                NeuronViewer nv3 = new NeuronViewer(
                    "G:\\data\\2005-08-25-0\\data002\\data002.params", false);
                ParametersFile pf3 = new ParametersFile(
                    "G:\\data\\2005-08-25-0\\data002\\data002.params");
                double x0 = ( (Double) pf3.evaluate("x0", new Integer(11))) *
                            (5.8 * 2) * 5;
                double y0 = ( (Double) pf3.evaluate("y0", new Integer(11))) *
                            (5.8 * 2) * 5;
//nv3.select("All/On/Y/T2/B/34");
                nv3.select("All/Off/Small/T2/OMS/11");
         PlotPanel panel = makeSTAPlot(nv3, 0, 0, 0, 3 * 72, (int) (1.8 * 72), labelFont,
                                              letterFont, "", x0, y0, 800);
                panel.setAxisVisible(false);
                panel.setScaleBar(500, 2, 0.05, 0.1);
                Frame f = new Frame();
                frame.setSize(3 * 72 + 25, (int) (1.8 * 72 + 40));
                frame.add(panel);
         frame.setVisible(true); //Click save as eps on this.  Next line does not work
//panel.saveAsEPS(new File("G:\\publications\\thesis\\Figures\\stimuli\\sampleSTA.eps"), 4*72, 2*72, false);

         */
        //Flash example panel
        /*       nv.select("/All/On/OscillatorA/460");
               PlotPanel p = (PlotPanel) nv.getPlot("flashPanel");
               p.setRange(0, 4, 0, 45);
               p.axesBorder.getVerticalAxis().setLabel("Spike Rate (Hz)");
               p.setPreferredSize(new Dimension(width, height));
               p.axesBorder.setRightPadding(rightPadding);
               p.axesBorder.setTopPadding(topPadding);
               p.axesBorder.getVerticalAxis().setThickness(verticalAxisThickness);
               p.axesBorder.getVerticalAxis().setLabelPosition(1);
               p.addBackgroundText("" + 'A', RIGHT, TOP, letterFont);
               p.setLabelFont(labelFont);
               c.gridx = 0;
               c.gridy = 0;
               frame.add(p, c);

               p = new PlotPanel();
               ScatterPlot flashPlot = new ScatterPlot();
               flashPlot.add(0, .48);
               flashPlot.add(1, .48);
               flashPlot.add(1.001, 0);
               flashPlot.add(2, 0);
               flashPlot.add(2.0001, -.48);
               flashPlot.add(3, -.48);
               flashPlot.add(3.0001, 0);
               flashPlot.add(4, 0);

               p.addData(flashPlot, new ScatterPlotStyle("",
         SymbolType.NONE, 3, Color.black, true,
                                                         Color.black, 1));

               p.setRange(0, 4, -.5, .5);

               p.setLegendVisible(false);
               p.setLabels("Time (s)", "Intensity(stimulus units)");
               p.setPreferredSize(new Dimension(width, (int) (height * .75)));
               p.axesBorder.setRightPadding(rightPadding);
               p.axesBorder.setTopPadding(topPadding);
               p.axesBorder.getVerticalAxis().setThickness(verticalAxisThickness);
               p.axesBorder.getVerticalAxis().setLabelPosition(1);
               p.addBackgroundText("" + 'B', RIGHT, TOP, letterFont);
               p.setLabelFont(labelFont);
               c.gridx = 0;
               c.gridy = 1;
               frame.add(p, c);

               frame.setVisible(true);

               RepaintManager.currentManager(frame).setDoubleBufferingEnabled(false);
               EpsGraphics2D g = new EpsGraphics2D(
         "", Font.decode("Arial ITALIC 8"), 0, 0, frame.getWidth(), frame.getHeight());

               g.scale(72.0 / dpi, 72.0 / dpi);
               g.setAccurateTextMode(true);
               frame.paint(g);

               FileWriter writer = new FileWriter(
                   "g:\\Publications\\thesis\\figures\\stimuli\\flashesExample.eps");
         //            RepaintManager.currentManager(frame).setDoubleBufferingEnabled(true);
         */
        /*               nv.select("/All/Off/Small/T2/DS/4/86");
                           PlotPanel p = (PlotPanel) nv.getPlot("GratingResponse");
                           p.setLabelFont(labelFont);

                           p.axesBorder.setRightPadding(rightPadding);
                       p.axesBorder.setTopPadding(topPadding);
                       p.setDoubleBuffered(false);
//Frame f = new Frame();
//         f.add(p);
//         f.setSize(6*72, 2*72);
//         f.setVisible(true);


                      p.saveAsEPS(new File(
         "g:\\Publications\\thesis\\figures\\stimuli\\sinExample1.eps"), 6, 2, false);


                         p = (PlotPanel) nv.getPlot("AverageGrating");
                           p.setLabelFont(labelFont);
                         p.axesBorder.setRightPadding(rightPadding);
                 p.axesBorder.setTopPadding(topPadding);
                 p.setDoubleBuffered(false);


                      p.saveAsEPS(new File(
         "g:\\Publications\\thesis\\figures\\stimuli\\sinExample2.eps"), 6, 2, false);

                         p = (PlotPanel) nv.getPlot("AverageGratingPanel2D");
                         p.setLabelFont(labelFont);
                         p.axesBorder.setRightPadding(rightPadding);
                p.axesBorder.setTopPadding(topPadding);
                p.setDoubleBuffered(false);


                      p.saveAsEPS(new File(
         "g:\\Publications\\thesis\\figures\\stimuli\\sinExample3.eps"), 6, 3, false);
         */

        nv.select("/All/On/Small/T2/OMS/174");
        PlotPanel p = (PlotPanel) nv.getPlot("omsScatter");
        p.setDoubleBuffered(false);
        p.axesBorder.setRightPadding(rightPadding);
        p.setLabels("Distance(\u03bcm)", "Spikes/Run");
        p.setLabelFont(labelFont);
        p.setRange(0, 1000, 0, p.getRange()[3] + 10);

        p.saveAsEPS(new File(
            "g:\\Publications\\thesis\\figures\\stimuli\\omsExample1.eps"), 4, 2, false);
        p.setLabels("Distance(\u03bcm)", "Spikes/Run");

        p = (PlotPanel) nv.getPlot("omsHist");
        p.setDoubleBuffered(false);
        p.axesBorder.setRightPadding(rightPadding);

        p.setLabelFont(labelFont);
        p.setRange(0, 1000, 0, p.getRange()[3] + 10);
        p.setLabels("Distance(\u03bcm)", "Spikes/Run");
        p.saveAsEPS(new File(
            "g:\\Publications\\thesis\\figures\\stimuli\\omsExample2.eps"), 4, 2, false);

        /*         nv.select("/All/On/TwoPeaks/XOMSA");
                p = (PlotPanel) nv.getPlot("omsHistograms");
                p.setLabels("Distance(\u03bcm)", "Normalized Spike Rate");
                      p.setLabelFont(labelFont);
                        p.setRange(0,1000, 0, p.getRange()[3] );

                 p.saveAsEPS(new File(
         "g:\\Publications\\thesis\\figures\\stimuli\\omsExample3.eps"), 4, 2, false);*/


        nv.select("All/Off/Large/T1/AOMS/9");
        p = (PlotPanel) nv.getPlot("F1F2");
        p.autoscale();
        p.setRange(p.getRange()[0], p.getRange()[1], 0, p.getRange()[3] + 10);
        p.setLabelFont(labelFont);
        p.axesBorder.getVerticalAxis().setThickness(verticalAxisThickness);
        p.axesBorder.getVerticalAxis().setLabelPosition(1);

        p.saveAsEPS(new File(
            "g:\\Publications\\thesis\\figures\\stimuli\\contrastReversingExample.eps"),
                    4, 2, false);

        //Sample white STA timecourse
        /*           nv = new NeuronViewer("g:\\data\\2005-08-25-0\\data002\\data002.params", false);
                   nv.select("All/Off/Small/T2/OMS/11");
                   PlotPanel p = (PlotPanel) nv.getPlot("tc");
         ScatterPlotStyle style = (ScatterPlotStyle) p.getStyleWithName("Blue Time Course");
                   style.setConnectionLineThickness(1.0f);
                    style = (ScatterPlotStyle) p.getStyleWithName("Red Time Course");
                   style.setConnectionLineThickness(1.0f);
                    style = (ScatterPlotStyle) p.getStyleWithName("Green Time Course");
                   style.setConnectionLineThickness(1.0f);
//              p.removeDataOfType(STATimeFunction1.class);
                   p.setLabels("Time before spike (ms)", "STA (stimulus units)" );
                   p.setLabelFont(labelFont);
                   p.axesBorder.setRightPadding(rightPadding);
                   p.axesBorder.setTopPadding(topPadding);
                   p.setDoubleBuffered(false);
                   PlotUtil.showData("", p);
                   p.saveAsEPS(new File(
         "g:\\Publications\\thesis\\figures\\stimuli\\whiteExampleSTA.eps"), 4, 2, false);
         */
        //Sample white STV timecourse
        /*         p = (PlotPanel) nv.getPlot("tcV");
                 style = (ScatterPlotStyle) p.getStyleWithName("Blue Time Course");
                 style.setConnectionLineThickness(1.0f);
                 style = (ScatterPlotStyle) p.getStyleWithName("Red Time Course");
                style.setConnectionLineThickness(1.0f);
                 style = (ScatterPlotStyle) p.getStyleWithName("Green Time Course");
                style.setConnectionLineThickness(1.0f);

//              p.removeDataOfType(STATimeFunction1.class);
                 p.setLabels("Time before spike (ms)", "STV (stimulus units)");
                 p.setLabelFont(labelFont);
                 p.axesBorder.setRightPadding(rightPadding);
                 p.axesBorder.setTopPadding(topPadding);
                 p.setDoubleBuffered(false);
                 PlotUtil.showData("", p);
                 p.saveAsEPS(new File(
         "g:\\Publications\\thesis\\figures\\stimuli\\whiteExampleSTV.eps"), 4, 2, false);
         */
    }


    public static PlotPanel makeMosaicPlot(NeuronViewer nv, int rightPadding,
                                           int topPadding, int verticalAxisThickness,
                                           int width, int height, Font labelFont,
                                           Font letterFont, String letter) {

        PlotPanel p = (PlotPanel) nv.getPlot("mosaic");
        p.setDoubleBuffered(false);
        p.setLegendVisible(false);
        p.addBackgroundText(letter, LEFT, TOP, letterFont, Color.red);
        p.setAxisVisible(true);
        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);
        p.axesBorder.getVerticalAxis().setThickness(verticalAxisThickness);
        p.axesBorder.getVerticalAxis().setLabelPosition(1);
        p.setPreferredSize(new Dimension(width, height));
        p.setLabelFont(labelFont);
        p.clearLegend();
        return p;

    }


    public static PlotPanel makeSTAPlot(NeuronViewer nv,
                                        int rightPadding, int topPadding,
                                        int verticalAxisThickness, int width, int height,
                                        Font labelFont, Font letterFont, String letter,
                                        double x0, double y0, double viewRadius) {

        PlotPanel p = (PlotPanel) nv.getPlot("sta");
        p.setDoubleBuffered(false);
        p.addBackgroundText(letter, LEFT, TOP, letterFont, Color.white);
        p.setRange(x0 - viewRadius, x0 + viewRadius, y0 - viewRadius / 2,
                   y0 + viewRadius / 2);
        p.setAxisVisible(true);
        p.setLabels("x (\u03BCm)", "y (\u03BCm)");
        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);
        p.axesBorder.getVerticalAxis().setThickness(verticalAxisThickness);
        p.axesBorder.getVerticalAxis().setLabelPosition(1);
        p.setPreferredSize(new Dimension(width, height));
        p.setLabelFont(labelFont);
        return p;

    }


    public static PlotPanel makeTCPlot(NeuronViewer nv, int rightPadding, int topPadding,
                                       int verticalAxisThickness, int width,
                                       int height, Font labelFont, Font letterFont,
                                       String letter) {

        PlotPanel p = (PlotPanel) nv.getPlot("tc");
        p.setDoubleBuffered(false);
        p.setLabels("Time before spike (ms)", "Normalized STA");
        p.addBackgroundText(letter, LEFT, TOP, letterFont, Color.red);
        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);
        p.axesBorder.getVerticalAxis().setThickness(verticalAxisThickness);
        p.axesBorder.getVerticalAxis().setLabelPosition(1);
        p.setPreferredSize(new Dimension(width, height));
        p.setLabelFont(labelFont);
        return p;

    }


    public static PlotPanel makeTCVPlot(NeuronViewer nv, int rightPadding, int topPadding,
                                        int verticalAxisThickness, int width,
                                        int height, Font labelFont, Font letterFont,
                                        String letter) {
        PlotPanel p = (PlotPanel) nv.getPlot("tcV");
        p.setDoubleBuffered(false);
        p.setLabels("Time before spike (ms)", "Normalized STV");
        p.addBackgroundText(letter, LEFT, TOP, letterFont, Color.red);
        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);
        p.axesBorder.getVerticalAxis().setThickness(verticalAxisThickness);
        p.axesBorder.getVerticalAxis().setLabelPosition(1);
        p.setPreferredSize(new Dimension(width, height));
        p.setLabelFont(labelFont);
        return p;
    }


    public static PlotPanel makeAverageGrating2DPlot(NeuronViewer nv, int rightPadding,
        int topPadding,
        int verticalAxisThickness, int width,
        int height, Font labelFont, Font letterFont, String letter) {

        PlotPanel p = (PlotPanel) nv.getPlot("AverageGratingPanel2D");
        p.setDoubleBuffered(false);
        p.addBackgroundText(letter, LEFT, TOP, letterFont, Color.red);
        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);
        p.axesBorder.getVerticalAxis().setThickness(verticalAxisThickness);
        p.axesBorder.getVerticalAxis().setLabelPosition(1);
        p.setPreferredSize(new Dimension(width, height));

        p.setLabelFont(labelFont);
        return p;

    }


    public static PlotPanel makeACFPlot(NeuronViewer nv, int rightPadding, int topPadding,
                                        int verticalAxisThickness, int width,
                                        int height, Font labelFont, Font letterFont,
                                        String letter) {
        PlotPanel p = (PlotPanel) nv.getPlot("acf");
        p.setDoubleBuffered(false);
        p.setLegendVisible(false);
        p.setLabels("Time difference (ms)", "Normalized # of Pairs");
        p.addBackgroundText(letter, LEFT, TOP, letterFont, Color.red);
        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);
        p.axesBorder.getVerticalAxis().setThickness(verticalAxisThickness);
        p.axesBorder.getVerticalAxis().setLabelPosition(1);
        p.setPreferredSize(new Dimension(width, height));
        p.setLabelFont(labelFont);
        return p;

    }


    public static PlotPanel makeF1F2Plot(NeuronViewer nv, int rightPadding,
                                         int topPadding,
                                         int verticalAxisThickness, int width,
                                         int height, Font labelFont, Font letterFont,
                                         String letter) {
        PlotPanel p = (PlotPanel) nv.getPlot("F1F2");
        p.setDoubleBuffered(false);
        p.setLabels("Spatial Frequency (cyc/mm)", "Ave. Spike Rate (Hz)");
        p.addBackgroundText(letter, LEFT, TOP, letterFont, Color.red);
        p.setAxesType(AxisType.LOG10, AxisType.LINEAR);
        p.autoscale();
        p.setRange(p.getRange()[0], p.getRange()[1], p.getRange()[2],
                   p.getRange()[3] + p.getRange()[3] * .2);
        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);
        p.axesBorder.getVerticalAxis().setThickness(verticalAxisThickness);
        p.axesBorder.getVerticalAxis().setLabelPosition(1);
        p.setPreferredSize(new Dimension(width, height));
        p.setLabelFont(labelFont);
        return p;

    }


    public static PlotPanel makeOMSPlot(NeuronViewer nv, int rightPadding, int topPadding,
                                        int verticalAxisThickness, int width,
                                        int height, Font labelFont, Font letterFont,
                                        String letter) {
        PlotPanel p = (PlotPanel) nv.getPlot("omsHistograms");
        p.setDoubleBuffered(false);
        p.setLabels("Distance (\u03bcm)", "Normalized Spike Rate");
        p.addBackgroundText(letter, LEFT, TOP, letterFont, Color.red);

        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);
        p.axesBorder.getVerticalAxis().setThickness(verticalAxisThickness);
        p.axesBorder.getVerticalAxis().setLabelPosition(1);
        p.setPreferredSize(new Dimension(width, height));
        p.setLabelFont(labelFont);
        return p;

    }


    public static PlotPanel makeFlashesPlot(NeuronViewer nv, int rightPadding,
                                            int topPadding,
                                            int verticalAxisThickness, int width,
                                            int height, Font labelFont, Font letterFont,
                                            String letter) {
        PlotPanel p = (PlotPanel) nv.getPlot("flashesPanel");
        p.setDoubleBuffered(false);
        p.setLabels("Time(s)", "Normalized Spike Rate");
        p.addBackgroundText(letter, LEFT, TOP, letterFont, Color.red);
        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);
        p.axesBorder.getVerticalAxis().setThickness(verticalAxisThickness);
        p.axesBorder.getVerticalAxis().setLabelPosition(1);
        p.setPreferredSize(new Dimension(width, height));
        p.setLabelFont(labelFont);
        return p;
    }


    public static void makeClassFigure(int nv, String className, String neuronToShow) throws
        Exception {
        if (mode == THESIS) makeClassFigureForThesis(nv, className, neuronToShow);
        else if (mode == POWERPOINT) makeClassFigureForPowerPoint(nv, className,
            neuronToShow);
    }


    public static void makeClassFigureForPowerPoint(int nv, String className,
        String neuronToShow) throws
        Exception {

        JFrame frame = new JFrame();

        frame.setLayout(new GridBagLayout());
        GridBagConstraints c = new GridBagConstraints();



        double dpi = 72.0;
        int width = (int) (3.2 * (dpi));
        int height = (int) (2.3 /*1.7*/ * (dpi));
        int staHeight = (int) (1.6 /*1.6*/ * (dpi));
        int verticalAxisThickness = 52;
        int rightPadding = 15;
        int topPadding = 4;
        Font labelFont = Font.decode("Arial ITALIC 10");
        Font letterFont = Font.decode("Arial Bold 12");
        frame.setSize(width * 3 + 8, height * 3 + 34);

        whitenedNVS[nv].showNeuronsInSubclasses = true;
        whitenedNVS[nv].select(className);
        c.gridx = 0;
        c.gridy = 0;
        frame.add(makeMosaicPlot(whitenedNVS[nv], rightPadding, topPadding,
                                 verticalAxisThickness,
                                 width, staHeight, labelFont, letterFont, "A"), c);

        whitenedNVS[nv].showNeuronsInSubclasses = false;

        whitenedNVS[nv].select(className + "/" + neuronToShow);

        c.gridx = 1;
        c.gridy = 0;
        double x0 = ( (Double) parameterFiles[nv].evaluate("x0", new Integer(neuronToShow))) *
                    5.8 * 5;
        double y0 = ( (Double) parameterFiles[nv].evaluate("y0", new Integer(neuronToShow))) *
                    5.8 * 5;

        frame.add(makeSTAPlot(whitenedNVS[nv], rightPadding, topPadding,
                              verticalAxisThickness,
                              width, staHeight, labelFont, letterFont, "B", x0, y0, 800),
                  c);



        whitenedNVS[nv].select(className);

        c.gridx = 2;
        c.gridy = 0;
        frame.add(makeTCPlot(whitenedNVS[nv], rightPadding, topPadding,
                             verticalAxisThickness,
                             width, height, labelFont, letterFont, "C"), c);

        nvs[nv].select(className);

      c.gridx = 0;
      c.gridy = 1;
      frame.add(makeTCPlot(nvs[nv], rightPadding, topPadding,
                           verticalAxisThickness,
                           width, height, labelFont, letterFont, "D"), c);


        c.gridx = 1;
        c.gridy = 1;
        frame.add(makeTCVPlot(nvs[nv], rightPadding, topPadding,
                              verticalAxisThickness,
                              width, height, labelFont, letterFont, "E"), c);

        c.gridx = 2;
        c.gridy = 1;
        frame.add(makeAverageGrating2DPlot(nvs[nv], rightPadding,
                                           topPadding, verticalAxisThickness, width,
                                           height, labelFont, letterFont, "F"), c);

        c.gridx = 0;
        c.gridy = 2;
        frame.add(makeACFPlot(nvs[nv], rightPadding, topPadding,
                              verticalAxisThickness,
                              width, height, labelFont, letterFont, "G"), c);

        c.gridx = 1;
        c.gridy = 2;
        frame.add(makeF1F2Plot(nvs[nv], rightPadding, topPadding,
                               verticalAxisThickness,
                               width, height, labelFont, letterFont, "H"), c);

        c.gridx = 2;
        c.gridy = 2;
        frame.add(makeOMSPlot(nvs[nv], rightPadding, topPadding,
                              verticalAxisThickness,
                              width, height, labelFont, letterFont, "I"), c);



        frame.setVisible(true);

        RepaintManager.currentManager(frame).setDoubleBufferingEnabled(false);
        BufferedImage bi = new BufferedImage(
            frame.getWidth(), frame.getHeight(), BufferedImage.TYPE_INT_RGB);
        Graphics2D g2 = bi.createGraphics();
        g2.setColor(Color.white);
        g2.fillRect(0, 0, frame.getWidth(), frame.getHeight());
        frame.paint(g2);
        ImageIO.write(bi, "PNG",
                      new File("g:\\Publications\\defense\\classes\\" +
                               className.replace("/", "") +
                               ".png"));

        RepaintManager.currentManager(frame).setDoubleBufferingEnabled(true);

//        frame.setVisible(false);
//        frame.dispose();

    }


    public static void makeClassFigureForThesis(int nv, String className,
                                                String neuronToShow) throws
        Exception {
        NeuronViewer currentNeuronViewer;
        JFrame frame = new JFrame();

        frame.setLayout(new GridBagLayout());
        GridBagConstraints c = new GridBagConstraints();

        double dpi = 72.0;
        int width = (int) (2.8 * (dpi));
        int height = (int) (1.45 /*1.7*/ * (dpi));
        int staHeight = (int) (1.6 * (dpi));
        int verticalAxisThickness = 42;
        int rightPadding = 9;
        int topPadding = 10;
        Font labelFont = Font.decode("Arial ITALIC 8");
        Font letterFont = Font.decode("Arial Bold 10");
        frame.setSize(width * 2 + 8, staHeight + height * 4 + 34);

        currentNeuronViewer = whitenedNVS[nv];
        currentNeuronViewer.showNeuronsInSubclasses = true;
        currentNeuronViewer.select(className);
        c.gridx = 0;
        c.gridy = 0;
        frame.add(makeMosaicPlot(currentNeuronViewer, rightPadding, topPadding,
                                 verticalAxisThickness,
                                 width, staHeight, labelFont, letterFont, "A"), c);
        currentNeuronViewer.showNeuronsInSubclasses = false;

        currentNeuronViewer.select(className + "/" + neuronToShow);

        c.gridx = 1;
        c.gridy = 0;

        double x0 = ( (Double) whitenedParameterFiles[nv].evaluate("x0",
            new Integer(neuronToShow))) * 5.8 * 5;
        double y0 = ( (Double) whitenedParameterFiles[nv].evaluate("y0",
            new Integer(neuronToShow))) * 5.8 * 5;
//       System.out.println(x0);
//       System.out.println(y0);
        frame.add(makeSTAPlot(currentNeuronViewer, rightPadding, topPadding,
                              verticalAxisThickness,
                              width, staHeight, labelFont, letterFont, "B", x0, y0, 800),
                  c);

        currentNeuronViewer.select(className);

        c.gridx = 0;
        c.gridy = 1;
        frame.add(makeTCPlot(currentNeuronViewer, rightPadding, topPadding,
                             verticalAxisThickness,
                             width, height, labelFont, letterFont, "C"), c);

        currentNeuronViewer = nvs[nv];
        currentNeuronViewer.showNeuronsInSubclasses = false;
        currentNeuronViewer.select(className);
        c.gridx = 1;
        c.gridy = 1;
        frame.add(makeTCPlot(currentNeuronViewer, rightPadding, topPadding,
                             verticalAxisThickness,
                             width, height, labelFont, letterFont, "D"), c);

        c.gridx = 0;
        c.gridy = 2;
        frame.add(makeTCVPlot(currentNeuronViewer, rightPadding, topPadding,
                              verticalAxisThickness,
                              width, height, labelFont, letterFont, "E"), c);

        c.gridx = 1;
        c.gridy = 2;
        frame.add(makeAverageGrating2DPlot(currentNeuronViewer, rightPadding,
                                           topPadding, verticalAxisThickness, width,
                                           height, labelFont, letterFont, "F"), c);

        c.gridx = 0;
        c.gridy = 3;
        frame.add(makeACFPlot(currentNeuronViewer, rightPadding, topPadding,
                              verticalAxisThickness,
                              width, height, labelFont, letterFont, "G"), c);

        c.gridx = 1;
        c.gridy = 3;
        frame.add(makeF1F2Plot(currentNeuronViewer, rightPadding, topPadding,
                               verticalAxisThickness,
                               width, height, labelFont, letterFont, "H"), c);

        c.gridx = 0;
        c.gridy = 4;
        frame.add(makeOMSPlot(currentNeuronViewer, rightPadding, topPadding,
                              verticalAxisThickness,
                              width, height, labelFont, letterFont, "I"), c);

        c.gridx = 1;
        c.gridy = 4;
        frame.add(makeFlashesPlot(currentNeuronViewer, rightPadding, topPadding,
                                  verticalAxisThickness,
                                  width, height, labelFont, letterFont, "J"), c);
        //Create white filler for the last spot.
//        JPanel jp = new JPanel();
//        jp.setPreferredSize(new Dimension(width, height));
//        jp.setBackground(Color.white);
//        c.gridx = 0;
//        c.gridy = 6;
//        frame.add(jp, c);

        frame.setVisible(true);

        RepaintManager.currentManager(frame).setDoubleBufferingEnabled(false);
        EpsGraphics2D g = new EpsGraphics2D(
            "", Font.decode("Arial ITALIC 8"), 0, 0, frame.getWidth(), frame.getHeight());

        g.scale(72.0 / dpi, 72.0 / dpi);
        g.setAccurateTextMode(true);
        frame.paint(g);

        RepaintManager.currentManager(frame).setDoubleBufferingEnabled(true);

        FileWriter writer = new FileWriter(
            "g:\\Publications\\thesis\\figures\\classes\\" + className.replace("/", "") +
            ".eps");

        writer.write(g.toString());
        writer.flush();
        writer.close();

        frame.setVisible(false);
        frame.dispose();

    }

public static void makeWhitenerFigure(String className, String neuronToShow) throws Exception {
    if(mode==THESIS) {
        makeWhitenerFigureForThesis(className, neuronToShow);
    } else {
        makeWhitenerFigureForPowerPoint(className, neuronToShow);
    }
}

public static void makeWhitenerFigureForPowerPoint(String className, String neuronToShow) throws
    Exception {
    NeuronViewer[] wnvs = new NeuronViewer[whitenerParamFilesString.length];
    ParametersFile[] wpfs = new ParametersFile[whitenerParamFilesString.length];

    for (int i = 0; i < wpfs.length; i++) {
        wpfs[i] = new ParametersFile(
            whitenerParamFilesString[i]);
        wnvs[i] = new NeuronViewer(whitenerParamFilesString[i], false, new Config("config.xml"));

    }

    JFrame frame = new JFrame();

    frame.setLayout(new GridBagLayout());
    GridBagConstraints c = new GridBagConstraints();

//    double dpi = 72.0;
//    int width = (int) (2.4 * (dpi));
//    int height = (int) (1.1 /*1.7*/ * (dpi));
//    int staHeight = (int) (1.2 * (dpi));
//    int verticalAxisThickness = 39;
//    int rightPadding = 9;
//    int topPadding = 6;
//    Font labelFont = Font.decode("Arial ITALIC 8");
//    frame.setSize(width * 2 + 8, staHeight * 3 + height * 3 + 34);

    double dpi = 72.0;
    int width = (int) (3.2 * (dpi));
    int height = (int) (1.6 /*1.7*/ * (dpi));
    int staHeight = (int) (1.6 /*1.6*/ * (dpi));
    int verticalAxisThickness = 50;
    int rightPadding = 15;
    int topPadding = 4;
    Font labelFont = Font.decode("Arial ITALIC 10");
    Font letterFont = Font.decode("Arial Bold 12");
    frame.setSize(width * 3 + 8, height * 3 + 34);


    //Make mosaic plot
    wnvs[0].showNeuronsInSubclasses = false;
    wnvs[0].select(className);

    c.gridx = 0;
    c.gridy = 0;
    frame.add(makeMosaicPlot(wnvs[0], rightPadding, topPadding,
                             verticalAxisThickness,
                             width, staHeight, labelFont, letterFont, "A")
              , c);

    wnvs[0].showNeuronsInSubclasses = false;

    //Make STA plot

    wnvs[0].select(className + "/" + neuronToShow);

    c.gridx = 1;
    c.gridy = 0;
    double x0 = ( (Double) wpfs[0].evaluate("x0", new Integer(neuronToShow))) *
                5.8 * 5;
    double y0 = ( (Double) wpfs[0].evaluate("y0", new Integer(neuronToShow))) *
                5.8 * 5;
    frame.add(makeSTAPlot(wnvs[0], rightPadding, topPadding,
                          verticalAxisThickness,
                          width, staHeight, labelFont, letterFont, "B", x0, y0, 900),
              c);

    wnvs[0].select(className);
    c.gridx = 2;
    c.gridy = 0;
    frame.add(makeTCPlot(wnvs[0], rightPadding, topPadding,
                        verticalAxisThickness,
                        width, height, labelFont, letterFont, "C"), c);





    //Make mosaic plot, second params
    wnvs[1].showNeuronsInSubclasses = false;
    wnvs[1].select(className);

    c.gridx = 0;
    c.gridy = 1;
    frame.add(makeMosaicPlot(wnvs[1], rightPadding, topPadding,
                             verticalAxisThickness,
                             width, staHeight, labelFont, letterFont, "D")
              , c);

    wnvs[1].showNeuronsInSubclasses = false;

    //Make STA plot


    wnvs[1].select(className + "/" + neuronToShow);


    c.gridx = 1;
    c.gridy = 1;
    frame.add(makeSTAPlot(wnvs[1], rightPadding, topPadding,
                          verticalAxisThickness,
                          width, staHeight, labelFont, letterFont, "E", x0, y0, 900),
              c);

    wnvs[1].select(className);
    c.gridx = 2;
  c.gridy = 1;
    frame.add(makeTCPlot(wnvs[1], rightPadding, topPadding,
                      verticalAxisThickness,
                      width, height, labelFont, letterFont, "F"), c);



    //Make mosaic plot, third params file
    wnvs[2].showNeuronsInSubclasses = false;
    wnvs[2].select(className);

    c.gridx = 0;
    c.gridy = 2;
    frame.add(makeMosaicPlot(wnvs[2], rightPadding, topPadding,
                             verticalAxisThickness,
                             width, staHeight, labelFont, letterFont, "G")
              , c);

    wnvs[2].showNeuronsInSubclasses = false;

    //Make STA plot

     wnvs[2].select(className + "/" + neuronToShow);
    c.gridx = 1;
    c.gridy = 2;
    frame.add(makeSTAPlot(wnvs[2], rightPadding, topPadding,
                          verticalAxisThickness,
                          width, staHeight, labelFont, letterFont, "H", x0, y0, 900),
              c);

    wnvs[2].select(className);
    c.gridx = 2;
    c.gridy = 2;
    frame.add(makeTCPlot(wnvs[2], rightPadding, topPadding,
                      verticalAxisThickness,
                      width, height, labelFont, letterFont, "I"), c);

         frame.setVisible(true);

    RepaintManager.currentManager(frame).setDoubleBufferingEnabled(false);

          BufferedImage bi = new BufferedImage(
              frame.getWidth(), frame.getHeight(), BufferedImage.TYPE_INT_RGB);
          Graphics2D g2 = bi.createGraphics();
          g2.setColor(Color.white);
          g2.fillRect(0, 0, frame.getWidth(), frame.getHeight());
          frame.paint(g2);
          ImageIO.write(bi, "PNG",
                        new File("g:\\Publications\\defense\\figures\\" +
                                 className.replace("/", "") +
                                 ".png"));

          RepaintManager.currentManager(frame).setDoubleBufferingEnabled(true);

//          frame.setVisible(false);
//          frame.dispose();


}


    public static void makeWhitenerFigureForThesis(String className, String neuronToShow) throws
        Exception {
        NeuronViewer[] wnvs = new NeuronViewer[whitenerParamFilesString.length];
        ParametersFile[] wpfs = new ParametersFile[whitenerParamFilesString.length];

        for (int i = 0; i < wpfs.length; i++) {
            wpfs[i] = new ParametersFile(
                whitenerParamFilesString[i]);
            wnvs[i] = new NeuronViewer(whitenerParamFilesString[i], false, new Config("config.xml"));

        }

        JFrame frame = new JFrame();

        frame.setLayout(new GridBagLayout());
        GridBagConstraints c = new GridBagConstraints();

        double dpi = 72.0;
        int width = (int) (2.4 * (dpi));
        int height = (int) (1.1 /*1.7*/ * (dpi));
        int staHeight = (int) (1.2 * (dpi));
        int verticalAxisThickness = 39;
        int rightPadding = 9;
        int topPadding = 6;
        Font labelFont = Font.decode("Arial ITALIC 8");
        frame.setSize(width * 2 + 8, staHeight * 3 + height * 3 + 34);

        //Make mosaic plot
        wnvs[0].showNeuronsInSubclasses = false;
        wnvs[0].select(className);

        c.gridx = 0;
        c.gridy = 0;
        frame.add(makeMosaicPlot(wnvs[0], rightPadding, topPadding,
                                 verticalAxisThickness,
                                 width, staHeight, labelFont, letterFont, "A")
                  , c);

        wnvs[0].showNeuronsInSubclasses = false;

        //Make STA plot

        wnvs[0].select(className + "/" + neuronToShow);

        c.gridx = 1;
        c.gridy = 0;
        double x0 = ( (Double) wpfs[0].evaluate("x0", new Integer(neuronToShow))) *
                    5.8 * 5;
        double y0 = ( (Double) wpfs[0].evaluate("y0", new Integer(neuronToShow))) *
                    5.8 * 5;
        frame.add(makeSTAPlot(wnvs[0], rightPadding, topPadding,
                              verticalAxisThickness,
                              width, staHeight, labelFont, letterFont, "B", x0, y0, 900),
                  c);

        wnvs[0].select(className);
        c.gridx = 0;
        c.gridy = 1;
        frame.add(makeTCPlot(wnvs[0], rightPadding, topPadding,
                            verticalAxisThickness,
                            width, height, labelFont, letterFont, "C"), c);





        //Make mosaic plot, second params
        wnvs[1].showNeuronsInSubclasses = false;
        wnvs[1].select(className);

        c.gridx = 0;
        c.gridy = 2;
        frame.add(makeMosaicPlot(wnvs[1], rightPadding, topPadding,
                                 verticalAxisThickness,
                                 width, staHeight, labelFont, letterFont, "D")
                  , c);

        wnvs[1].showNeuronsInSubclasses = false;

        //Make STA plot


        wnvs[1].select(className + "/" + neuronToShow);

        c.gridx = 1;
        c.gridy = 0;

        c.gridx = 1;
        c.gridy = 2;
        frame.add(makeSTAPlot(wnvs[1], rightPadding, topPadding,
                              verticalAxisThickness,
                              width, staHeight, labelFont, letterFont, "E", x0, y0, 900),
                  c);

        wnvs[1].select(className);
        c.gridx = 0;
      c.gridy = 3;
        frame.add(makeTCPlot(wnvs[1], rightPadding, topPadding,
                          verticalAxisThickness,
                          width, height, labelFont, letterFont, "F"), c);



        //Make mosaic plot, third params file
        wnvs[2].showNeuronsInSubclasses = false;
        wnvs[2].select(className);

        c.gridx = 0;
        c.gridy = 4;
        frame.add(makeMosaicPlot(wnvs[2], rightPadding, topPadding,
                                 verticalAxisThickness,
                                 width, staHeight, labelFont, letterFont, "G")
                  , c);

        wnvs[2].showNeuronsInSubclasses = false;

        //Make STA plot

         wnvs[2].select(className + "/" + neuronToShow);
        c.gridx = 1;
        c.gridy = 4;
        frame.add(makeSTAPlot(wnvs[2], rightPadding, topPadding,
                              verticalAxisThickness,
                              width, staHeight, labelFont, letterFont, "H", x0, y0, 900),
                  c);

        wnvs[2].select(className);
        c.gridx = 0;
      c.gridy = 5;
        frame.add(makeTCPlot(wnvs[2], rightPadding, topPadding,
                          verticalAxisThickness,
                          width, height, labelFont, letterFont, "I"), c);

             frame.setVisible(true);

        RepaintManager.currentManager(frame).setDoubleBufferingEnabled(false);
        EpsGraphics2D g = new EpsGraphics2D(
            "", Font.decode("Arial ITALIC 8"), 0, 0, frame.getWidth(), frame.getHeight());

        g.scale(72.0 / dpi, 72.0 / dpi);
        g.setAccurateTextMode(true);

        frame.paint(g);

        RepaintManager.currentManager(frame).setDoubleBufferingEnabled(true);

        FileWriter writer = new FileWriter(
            "g:\\Publications\\thesis\\figures\\results\\" + className.replace("/", "") +
            ".eps");

        writer.write(g.toString());
        writer.flush();
        writer.close();
//        frame.setVisible(true);
//        frame.dispose();

    }


    public static void makeSeparationQualityFigure(String[] whitenerParamFilesString,
        int preparationNumber) throws
        Exception {

        NeuronViewer[] wnvs = new NeuronViewer[whitenerParamFilesString.length];
        for (int i = 0; i < whitenerParamFilesString.length; i++) {
            wnvs[i] = new NeuronViewer(whitenerParamFilesString[i], false, new Config("config.xml"));
        }

        GridBagConstraints c = new GridBagConstraints();

        double dpi = 72.0;
        int width = (int) (3.0 * (dpi));
        int height = (int) (3.0 * (dpi));

        int verticalAxisThickness = 39;
        int rightPadding = 9;
        int topPadding = 4;
        Font labelFont = Font.decode("Arial ITALIC 8");

        PlotPanel p = new PlotPanel();

        wnvs[0].select("All");
        CalculateClassSeparationAction action = new CalculateClassSeparationAction();
        action.initialize(wnvs[0]);
        TreePath path = wnvs[0].rightTree.getSelectionPath();

        double[] separations = action.calculateClassSeparation(path);
        int bins = 1000;
        double max = 10;
        double binSize = max / bins;
        double[] x = new double[bins];
        double[] values = new double[bins];
        for (int i = 0; i < x.length; i++) {
            x[i] = i * binSize;
            for (int j = 0; j < separations.length; j++) {
                if (separations[j] > x[i]) values[i]++;
            }
        }

        MathUtil.divide(values, MathUtil.max(values));
        ScatterPlot scatter = new ScatterPlot(x, values, new double[x.length]);
        ScatterPlotStyle style = new ScatterPlotStyle("Scatter 1", SymbolType.NONE, 0,
            Color.black, true, Color.black, 1f);
        p.addData(scatter, style);
        scatter = new ScatterPlot(new double[] {4, 5}, new double[] {1.08, 1.08}, null);
        p.addData(scatter, style);
        p.addBackgroundText("60 min., natural-power", 120, 5,
                            Font.decode("Arial Italic 8"));

        wnvs[1].select("All");
        action = new CalculateClassSeparationAction();
        action.initialize(wnvs[1]);
        path = wnvs[1].rightTree.getSelectionPath();
        separations = action.calculateClassSeparation(path);
        x = new double[bins];
        values = new double[bins];
        for (int i = 0; i < x.length; i++) {
            x[i] = i * binSize;
            for (int j = 0; j < separations.length; j++) {
                if (separations[j] > x[i]) values[i]++;
            }
        }
        MathUtil.divide(values, MathUtil.max(values));

        scatter = new ScatterPlot(x, values, new double[values.length]);
        style = new ScatterPlotStyle("Scatter 2", SymbolType.NONE, 0, Color.black, true,
                                     Color.black, 1f);
        style.setDashPattern("5,5");
        p.addData(scatter, style);
        scatter = new ScatterPlot(new double[] {4, 5}, new double[] {1.047, 1.047}, null);
        p.addData(scatter, style);
        p.addBackgroundText("30 min., natural-power", 120, 15,
                            Font.decode("Arial Italic 8"));

        wnvs[2].select("All");
        action = new CalculateClassSeparationAction();
        action.initialize(wnvs[2]);
        path = wnvs[2].rightTree.getSelectionPath();
        separations = action.calculateClassSeparation(path);
        x = new double[bins];
        values = new double[bins];
        for (int i = 0; i < x.length; i++) {
            x[i] = i * binSize;
            for (int j = 0; j < separations.length; j++) {
                if (separations[j] > x[i]) values[i]++;
            }
        }
        MathUtil.divide(values, MathUtil.max(values));

        scatter = new ScatterPlot(x, values, new double[values.length]);
        style = new ScatterPlotStyle("Scatter 3", SymbolType.NONE, 0, Color.black, true,
                                     Color.black, 1f);
        style.setDashPattern("2,2");
        p.addData(scatter, style);
        scatter = new ScatterPlot(new double[] {4, 5}, new double[] {1.015, 1.015}, null);
        p.addData(scatter, style);
        p.addBackgroundText("30 min., white", 120, 25, Font.decode("Arial Italic 8"));

        p.setLabels("Separation Parameter", "Fraction of Paired Types Separated");
        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);
//          p.axesBorder.getVerticalAxis().setThickness(verticalAxisThickness);
//          p.axesBorder.getVerticalAxis().setLabelPosition(1);
        p.setPreferredSize(new Dimension(width, height));
        p.setLabelFont(labelFont);
        p.autoscale();

        p.setRange(0, 10, .5 /*p.getRange()[2]*/, 1.1);

        PlotUtil.showData("", p);

        p.saveAsEPS(new File(
            "g:\\Publications\\thesis\\figures\\results\\separation" +
            preparationNumber + ".eps"), 3, 3, false);

    }


    public static void makeClassFigures() throws Exception {
//        makeClassFigure(0, "All/On/Small/T2/OMS", "67");
        makeClassFigure(2, "All/On/Small/T3/OMS", "264");
        makeClassFigure(2, "All/On/Large/T1/Triphasic", "5");
        makeClassFigure(0, "All/On/Large/T2/DS/3", "24");
        makeClassFigure(0, "All/On/Large/T2/AOMS", "12");
        makeClassFigure(0, "All/On/Large/T2/Composite", "253");
        makeClassFigure(2, "All/On/Large/T2/A", "188");
        makeClassFigure(0, "All/On/Large/T3/UniformityA", "175");
        makeClassFigure(2, "All/On/Large/T3/UniformityB", "408");
        makeClassFigure(1, "All/On/Large/T3/Blue", "508");
        makeClassFigure(1, "All/On/Large/T3/OS", "240");
        makeClassFigure(2, "All/On/Large/T3/Triphasic", "491");
        makeClassFigure(0, "All/On/Large/T3/A", "142");
        makeClassFigure(2, "All/On/Large/T3/B", "180");
        makeClassFigure(0, "All/On/Large/T3/C", "370");
        makeClassFigure(0, "All/On/Large/T3/D", "156");
        makeClassFigure(2, "All/On/Large/T3/E", "61");
        makeClassFigure(0, "All/On/Large/T3/F", "258");
        makeClassFigure(0, "All/Off/Small/T2/DS/2", "171");
        makeClassFigure(0, "All/Off/Small/T2/OMS", "28");
        makeClassFigure(2, "All/Off/Small/T3/OMS", "229");
        makeClassFigure(2, "All/Off/Small/T3/A", "99");
        makeClassFigure(0, "All/Off/Small/T3/B", "488");
        makeClassFigure(1, "All/Off/Small/T3/C", "489");
        makeClassFigure(0, "All/Off/Large/T1/AOMS", "16");
        makeClassFigure(0, "All/Off/Large/T1/A", "10");
        makeClassFigure(2, "All/Off/Large/T3/A", "463");

    }


    public static void calculateReferenceParameters() throws Exception {
        for (int i = 0; i < parameterFiles.length; i++) {

            DefaultMutableTreeNode node = whitenedNVS[i].select(referenceType);
            IntegerList list = InteractiveTree.getNeuronsInClass(node,
                new IntegerList(), false);
            int[] neuronsInClass = list.toArray();

            referenceSimplifiedContourArea[i] = calculateClassParameter(neuronsInClass, i,
                "2*((simpleContourArea/3.14159)^.5)"
                                                /* "5.8*5*2*((SigmaX*SigmaY)^.5)"*/).x;
            System.out.println("Simplified: " + i + " " +
                               calculateClassParameter(neuronsInClass, i,
                "2*((simpleContourArea/3.14159)^.5)").toString(true));
//            System.out.println("Gauss: " +calculateClassParameter(neuronsInClass, i,
//
//                                "5.8*5*2*((SigmaX*SigmaY)^.5)").x);
//
            referenceContourArea[i] = calculateClassParameter(neuronsInClass, i,
                "2*((contourArea/3.14159)^.5)"
                                      /* "5.8*5*2*((SigmaX*SigmaY)^.5)"*/).x;
            System.out.println("Contour: " + i + " " +
                               calculateClassParameter(neuronsInClass, i,
                "2*((contourArea/3.14159)^.5)").toString(true));
            referenceRL[i] = calculateClassParameter(neuronsInClass, i, "rl").x;
            System.out.println("rl: " + i + " " +
                               calculateClassParameter(neuronsInClass, i, "rl").toString(true));

        }

        referenceAxon[0] = .835;
        referenceAxon[1] = .914;
        referenceAxon[2] = .987;
    }


    public static Num calculateClassParameter(String className, int parameterFile,
                                              String expression) {
        if (!expression.matches("Axon Speed")) {
            throw new IllegalStateException();
        }

        if (className.matches("All/On/Small/T2/OMS")) {
            switch (parameterFile) {
                case 0:
                    return new Num(.885, .018);
                case 1:
                    return new Num(.819, .019);
                case 2:
                    return new Num(.886, .025);
                default:
                    throw new IllegalStateException("Preparation " + parameterFile +
                        " not found for class " + className + ".");
            }
        } else if (className.matches("All/On/Small/T3/OMS")) {
            switch (parameterFile) {
                case 0:
                    return new Num(.622, .002);
                case 2:
                    return new Num(.684, .019);
                default:
                    throw new IllegalStateException("Preparation " + parameterFile +
                        " not found for class " + className + ".");
            }
        } else if (className.matches("All/On/Large/T3/E")) {
            switch (parameterFile) {
                case 2:
                    return new Num(.616, .029);
                default:
                    throw new IllegalStateException("Preparation " + parameterFile +
                        " not found for class " + className + ".");
            }
        } else if (className.matches("All/On/Large/T2/DS")) {
            switch (parameterFile) {
                case 0:
                    return new Num(.783, .033);
                case 1:
                    return new Num(.891, .064);
                case 2:
                    return new Num(1.008, .031);
                default:
                    throw new IllegalStateException("Preparation " + parameterFile +
                        " not found for class " + className + ".");
            }
        } else if (className.matches("All/On/Large/T1/Triphasic")) {
            switch (parameterFile) {
                case 2:
                    return new Num(1.142, .066);
                default:
                    throw new IllegalStateException("Preparation " + parameterFile +
                        " not found for class " + className + ".");
            }
        } else if (className.matches("All/On/Large/T2/A")) {
            switch (parameterFile) {
                case 0:
                    return new Num(.775, .003);
                case 2:
                    return new Num(.75, .044);
                default:
                    throw new IllegalStateException("Preparation " + parameterFile +
                        " not found for class " + className + ".");
            }
        } else if (className.matches("All/On/Large/T2/AOMS")) {
            switch (parameterFile) {
                case 0:
                    return new Num(.902, .026);
                case 1:
                    return new Num(.9, .026);
                case 2:
                    return new Num(.866, .042);
                default:
                    throw new IllegalStateException("Preparation " + parameterFile +
                        " not found for class " + className + ".");
            }
        } else if (className.matches("All/On/Large/T2/Composite")) {
            switch (parameterFile) {
                case 0:
                    return new Num(.826, .034);
                case 1:
                    return new Num(.908, .027);
                default:
                    throw new IllegalStateException("Preparation " + parameterFile +
                        " not found for class " + className + ".");
            }
        } else if (className.matches("All/On/Large/T3/UniformityA")) {
            switch (parameterFile) {
                case 0:
                    return new Num(.775, .021);
                case 1:
                    return new Num(.743, .016);
                default:
                    throw new IllegalStateException("Preparation " + parameterFile +
                        " not found for class " + className + ".");
            }
        } else if (className.matches("All/On/Large/T3/Blue")) {
            switch (parameterFile) {
                case 0:
                    return new Num(.63, .031);
                case 1:
                    return new Num(.639, .018);
                default:
                    throw new IllegalStateException("Preparation " + parameterFile +
                        " not found for class " + className + ".");
            }
        } else if (className.matches("All/On/Large/T3/OS")) {
            switch (parameterFile) {
                case 1:
                    return new Num(.634, .032);
                default:
                    throw new IllegalStateException("Preparation " + parameterFile +
                        " not found for class " + className + ".");
            }
        } else if (className.matches("All/On/Large/T3/Triphasic")) {
            switch (parameterFile) {

                case 1:
                    return new Num(0, Double.POSITIVE_INFINITY);
                case 2:
                    return new Num(.631, .056);
                default:
                    throw new IllegalStateException("Preparation " + parameterFile +
                        " not found for class " + className + ".");
            }
        } else if (className.matches("All/On/Large/T3/UniformityB")) {
            switch (parameterFile) {
                case 2:
                    return new Num(.448, .066);
                default:
                    throw new IllegalStateException("Preparation " + parameterFile +
                        " not found for class " + className + ".");
            }
        } else if (className.matches("All/On/Large/T3/A")) {
            switch (parameterFile) {
                case 0:
                    return new Num(.707, .02);
                case 1:
                    return new Num(.665, .079);
                case 2:
                    return new Num(.695, .021);
                default:
                    throw new IllegalStateException("Preparation " + parameterFile +
                        " not found for class " + className + ".");
            }
        } else if (className.matches("All/On/Large/T3/C")) {
            switch (parameterFile) {
                case 0:
                    return new Num(.521, .019);
                case 2:
                    return new Num(.676, .021);
                default:
                    throw new IllegalStateException("Preparation " + parameterFile +
                        " not found for class " + className + ".");
            }
        } else if (className.matches("All/On/Large/T3/B")) {
            switch (parameterFile) {
                case 2:
                    return new Num(.795, .032);
                default:
                    throw new IllegalStateException("Preparation " + parameterFile +
                        " not found for class " + className + ".");
            }
        } else if (className.matches("All/On/Large/T3/D")) {
            switch (parameterFile) {
                case 0:
                    return new Num(.496, .032);
                case 1:
                    return new Num(.517, .037);
                case 2:
                    return new Num(.641, .023);
                default:
                    throw new IllegalStateException("Preparation " + parameterFile +
                        " not found for class " + className + ".");
            }
        } else if (className.matches("All/On/Large/T3/F")) {
            switch (parameterFile) {
                case 0:
                    return new Num(.653, .054);
                case 1:
                    return new Num(.638, .036);
                default:
                    throw new IllegalStateException("Preparation " + parameterFile +
                        " not found for class " + className + ".");
            }
        } else if (className.matches("All/Off/Small/T2/DS")) {
            switch (parameterFile) {
                case 0:
                    return new Num(.825, .016);
                case 1:
                    return new Num(.883, .015);
                case 2:
                    return new Num(.864, .017);
                default:
                    throw new IllegalStateException("Preparation " + parameterFile +
                        " not found for class " + className + ".");
            }
        } else if (className.matches("All/Off/Small/T2/OMS")) {
            switch (parameterFile) {
                case 0:
                    return new Num(.835, .012);
                case 1:
                    return new Num(.914, .025);
                case 2:
                    return new Num(.987, .014);
                default:
                    throw new IllegalStateException("Preparation " + parameterFile +
                        " not found for class " + className + ".");
            }
        } else if (className.matches("All/Off/Small/T3/A")) {
            switch (parameterFile) {
                case 0:
                    return new Num(.613, .025);
                case 1:
                    return new Num(0, Double.POSITIVE_INFINITY);
                case 2:
                    return new Num(.745, .032);
                default:
                    throw new IllegalStateException("Preparation " + parameterFile +
                        " not found for class " + className + ".");
            }
        } else if (className.matches("All/Off/Large/T3/A")) {
            switch (parameterFile) {
                case 2:
                    return new Num(.587, .053);
                default:
                    throw new IllegalStateException("Preparation " + parameterFile +
                        " not found for class " + className + ".");
            }
        } else if (className.matches("All/Off/Small/T3/B")) {
            switch (parameterFile) {
                case 0:
                    return new Num(.661, .052);
                case 2:
                    return new Num(.676, .011);
                default:
                    throw new IllegalStateException("Preparation " + parameterFile +
                        " not found for class " + className + ".");
            }
        } else if (className.matches("All/Off/Small/T3/OMS")) {
            switch (parameterFile) {
                case 0:
                    return new Num(.513, .032);
                case 2:
                    return new Num(.585, .02);
                default:
                    throw new IllegalStateException("Preparation " + parameterFile +
                        " not found for class " + className + ".");
            }
        } else if (className.matches("All/Off/Large/T1/AOMS")) {
            switch (parameterFile) {
                case 0:
                    return new Num(1.123, .036);
                case 1:
                    return new Num(1.015, .036);
                case 2:
                    return new Num(1.182, .049);
                default:
                    throw new IllegalStateException("Preparation " + parameterFile +
                        " not found for class " + className + ".");
            }
        } else if (className.matches("All/Off/Large/T1/A")) {
            switch (parameterFile) {
                case 0:
                    return new Num(.924, .03);
                case 1:
                    return new Num(1.015, .033);
                default:
                    throw new IllegalStateException("Preparation " + parameterFile +
                        " not found for class " + className + ".");
            }
        } else if (className.matches("All/Off/Small/T3/C")) {
            switch (parameterFile) {
                case 0:
                    return new Num(.575, .054);
                case 1:
                    return new Num(.542, .028);
                default:
                    throw new IllegalStateException("Preparation " + parameterFile +
                        " not found for class " + className + ".");
            }
        }

        else {
            throw new IllegalStateException("Class " + className +
                                            " not found in axon speed database.");
        }

    }


    public static Num calculateClassParameter(int[] neuronsInClass, int parameterFile,
                                              String expression) throws Exception {

        if (expression.equals("F2 over F1")) {

            return calculateNonlinearity(neuronsInClass,
                                         whitenedParameterFiles[parameterFile]);
        } else if (expression.equals("OMS")) {
            return OMSClassification.calculateOMSParameter(neuronsInClass,
                parameterFiles[parameterFile], 0, 150, 500, 800);

        } else {

            MeanVarianceCalculator mvc = new MeanVarianceCalculator(
                MeanVarianceCalculator.UNBIASED);

            HashMap values = whitenedParameterFiles[parameterFile].evaluate(expression,
                neuronsInClass);
            for (int j = 0; j < neuronsInClass.length; j++) {
                double value = ( (Double) values.get(neuronsInClass[j])).doubleValue();
                if (! (new Double(value).isNaN())) {
                    mvc.add(value);
                }
            }

            Num value = new Num(mvc.getMean(), mvc.getMeanVariance());
            //           System.out.println("N: " + neuronsInClass.length);
            //           System.out.println("Value: " + value);
            return value;
        }
    }


    public static void calculateParameters(String[] classNames) throws Exception {

        int rightPadding = 8;
        int topPadding = 4;

        ScatterPlot scats[] = new ScatterPlot[6];
        Gaussian1DFunction funcs[][] = new Gaussian1DFunction[6][27];
        for (int i = 0; i < scats.length; i++) {
            scats[i] = new ScatterPlot();

        }

        JPanel jpanel = new JPanel();
        jpanel.setLayout(new GridLayout(3, 2));
        DoubleHistogram hist = new DoubleHistogram("", 0, 20, 2);

        int currentRow = 0, currentColumn = 1;
        Table table2 = new Table(100, 11);
        table2.setTitle(0, "");

        table2.setTitle(currentColumn++, "N");
        table2.setTitle(currentColumn++, "Contour Diameter");
        table2.setTitle(currentColumn++, "Simplified Contour Diameter");
        table2.setTitle(currentColumn++, "Axon Speed");
        table2.setTitle(currentColumn++, "RL");
        table2.setTitle(currentColumn++, "SRM");
        table2.setTitle(currentColumn++, "F2/F1");
        table2.setTitle(currentColumn++, "OMS");
        table2.setTitle(currentColumn++, "blueness");
        table2.setTitle(currentColumn++, "magOS");

        for (int i = 0; i < classNames.length; i++) {
            ArrayList<Num> contourList = new ArrayList<Num> ();
            ArrayList<Num> simplifiedContourList = new ArrayList<Num> ();
            ArrayList<Num> axonList = new ArrayList<Num> ();
            ArrayList<Num> srmList = new ArrayList<Num> ();
            ArrayList<Num> f2f1List = new ArrayList<Num> ();
            ArrayList<Num> omsList = new ArrayList<Num> ();
            ArrayList<Num> rlList = new ArrayList<Num> ();
            ArrayList<Num> blueList = new ArrayList<Num> ();
            ArrayList<Num> osList = new ArrayList<Num> ();

            int nNeurons = 0;
            for (int j = 0; j < parameterFiles.length; j++) {
                table2.setCell(currentRow, 0, classNames[i] + "   " + (j + 1));
                currentColumn = 1;
                DefaultMutableTreeNode node = whitenedNVS[j].select(classNames[i]);
                if (node != null) {
                    IntegerList list = new IntegerList();
                    //If on a DS class, get the children folders
                    if (classNames[i].contains("DS")) {
                        for (int k = 0; k < node.getChildCount(); k++) {
                            DefaultMutableTreeNode dsNode = (DefaultMutableTreeNode) node.
                                getChildAt(k);
                            list.add(InteractiveTree.getNeuronsInClass(dsNode,
                                new IntegerList(), false));
                        }
                    } else {
                        list = InteractiveTree.getNeuronsInClass(node,
                            new IntegerList(), false);
                    }
                    int[] neuronsInClass = list.toArray();
                    nNeurons += neuronsInClass.length;
                    if (neuronsInClass.length != 0) {

                        table2.setCell(currentRow, currentColumn++, neuronsInClass.length);

                        Num num = calculateClassParameter(neuronsInClass, j,
                            "2*((contourArea/3.14159)^.5)"
                                  /* "5.8*5*2*((SigmaX*SigmaY)^.5)"*/).div(
                                      referenceContourArea[j]);
                        contourList.add(num);
                        table2.setCell(currentRow, currentColumn++, num.toString(true));

                        num = calculateClassParameter(neuronsInClass, j,
                            "2*((simpleContourArea/3.14159)^.5)"
                              /* "5.8*5*2*((SigmaX*SigmaY)^.5)"*/).div(
                                  referenceSimplifiedContourArea[j]);
                        simplifiedContourList.add(num);
                        table2.setCell(currentRow, currentColumn++, num.toString(true));

                        num = calculateClassParameter(classNames[i], j,
                            "Axon Speed").div(referenceAxon[j]);
                        //If quantity was measured.
                        if (num.err() < Double.POSITIVE_INFINITY) {
                            axonList.add(num);
                            table2.setCell(currentRow, currentColumn++, num.toString(true));
                        } else {
                            table2.setCell(currentRow, currentColumn++, "N/A");
                        }

                        num = calculateClassParameter(neuronsInClass, j,
                            "rl").div(referenceRL[j]);
                        rlList.add(num);
                        table2.setCell(currentRow, currentColumn++, num.toString(true));
                        num = calculateClassParameter(neuronsInClass, j,
                            "srm");
                        srmList.add(num);
                        table2.setCell(currentRow, currentColumn++, num.toString(true));

                        num = calculateClassParameter(neuronsInClass, j, "F2 over F1");
                        f2f1List.add(num);
                        table2.setCell(currentRow, currentColumn++, num.toString(true));

                        num = calculateClassParameter(neuronsInClass, j, "OMS");

                        omsList.add(num);
                        table2.setCell(currentRow, currentColumn++, num.toString(true));

                        num = calculateClassParameter(neuronsInClass, j,
                            "blueness");
                        blueList.add(num);
                        table2.setCell(currentRow, currentColumn++, num.toString(true));

                        num = calculateClassParameter(neuronsInClass, j,
                            "magOS");
                        osList.add(num);
                        table2.setCell(currentRow++, currentColumn++, num.toString(true));

                    } else {
                        System.err.println("No neurons for class " + classNames[i] +
                                           " found in preparation " + j);
                    }
                } else System.err.println("Class " + classNames[i] +
                                          " not found in preparation " + (j + 1));
            }

            hist.fill(Num.averageRescaleErrors(f2f1List).x, 1.0);

            int classNumber = i + 1;
            if (classNumber > 4) classNumber += 2;
            if (classNumber > 21) classNumber += 3;

            Num num = Num.averageRescaleErrors(simplifiedContourList);
            if (classNames[i].equals(referenceType)) {
                num = new Num(num.x, 0.000);
            }
            scats[0].add(num.x, classNumber, 0, num.err());
            funcs[0][i] = new Gaussian1DFunction(1, num.x, num.err());

            num = Num.averageRescaleErrors(axonList);
            if (classNames[i].equals(referenceType)) {
                num = new Num(num.x, 0.000);
            }
            scats[1].add(num.x, classNumber, 0, num.err());
            funcs[1][i] = new Gaussian1DFunction(1, num.x, num.err());

            num = Num.averageRescaleErrors(rlList);
            if (classNames[i].equals(referenceType)) {
                num = new Num(num.x, 0.000);
            }
            scats[2].add(num.x, classNumber, 0, num.err());
            funcs[2][i] = new Gaussian1DFunction(1, num.x, num.err());

            num = Num.averageRescaleErrors(srmList);
            scats[3].add(num.x, classNumber, 0, num.err());
            funcs[3][i] = new Gaussian1DFunction(1, num.x, num.err());

            num = Num.averageRescaleErrors(f2f1List);
            scats[4].add(num.x, classNumber, 0, num.err());
            funcs[4][i] = new Gaussian1DFunction(1, num.x, num.err());

            num = Num.averageRescaleErrors(omsList);
            scats[5].add(num.x, classNumber, 0, num.err());
            funcs[5][i] = new Gaussian1DFunction(1, num.x, num.err());

            currentColumn = 1;
            table2.setCell(currentRow, 0, classNames[i] + " Tot");
            table2.setCell(currentRow, currentColumn++, nNeurons);
            table2.setCell(currentRow, currentColumn++,
                           Num.averageRescaleErrors(contourList).toString(true));
            table2.setCell(currentRow, currentColumn++,
                           Num.averageRescaleErrors(simplifiedContourList).toString(true));
            table2.setCell(currentRow, currentColumn++,
                           Num.averageRescaleErrors(axonList).toString(true));
            table2.setCell(currentRow, currentColumn++,
                           Num.averageRescaleErrors(rlList).toString(true));
            table2.setCell(currentRow, currentColumn++,
                           Num.averageRescaleErrors(srmList).toString(true));
            table2.setCell(currentRow, currentColumn++,
                           Num.averageRescaleErrors(f2f1List).toString(true));
            table2.setCell(currentRow, currentColumn++,
                           Num.averageRescaleErrors(omsList).toString(true));
            table2.setCell(currentRow, currentColumn++,
                           Num.averageRescaleErrors(blueList).toString(true));
            table2.setCell(currentRow++, currentColumn++,
                           Num.averageRescaleErrors(osList).toString(true));
        }

        PlotUtil.showData("", hist);
        PlotPanel p = new PlotPanel();
        p.addData(scats[0],
                  new ScatterPlotStyle("Scatter", SymbolType.FILLED_SQUARE, 3,
                                       Color.black));
        FunctionSum f = new FunctionSum(funcs[0]);
        p.addData(f, new FunctionStyle("f", Color.red, 1));

        ScatterPlot scatter = new ScatterPlot();
        scatter.add(1.28, 0);
        scatter.add(1.28, 33);
        ScatterPlotStyle divider = new ScatterPlotStyle("Divider", SymbolType.NONE, 0,
            Color.black, true, Color.black, 1);
        divider.setDashPattern("5, 5");
        p.addData(scatter, divider);

//        p.addBackgroundText("" + 'A', LEFT, TOP, letterFont, Color.black);
        p.setLabels("Relative Simplified-Contour Effective-Diameter", "Type Number");
        p.setLabelFont(labelFont);
        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);
        p.setDoubleBuffered(false);

        p.setRange(.5, 4.5, 0, 33);
        jpanel.add(p);
        p = new PlotPanel();
        p.addData(scats[1],
                  new ScatterPlotStyle("Scatter", SymbolType.FILLED_SQUARE, 3,
                                       Color.black));
        f = new FunctionSum(funcs[1]);

        p.addData(f, new FunctionStyle("f", Color.red, 1));

//       p.addBackgroundText("" + 'B', LEFT, TOP, letterFont, Color.black);

        scatter = new ScatterPlot();
        scatter.add(.87, 0);
        scatter.add(.87, 33);
        p.addData(scatter, divider);

        scatter = new ScatterPlot();
        scatter.add(1.06, 0);
        scatter.add(1.06, 33);
        p.addData(scatter, divider);

        p.setLabels("Relative Axon Speed", "Type Number");
        p.setLabelFont(labelFont);
        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);
        p.setDoubleBuffered(false);
        p.setRange(.3, 1.3, 0, 33);
        jpanel.add(p);

//        PlotPanel axonsOnly = new PlotPanel();
//        axonsOnly.addData(scats[1],
//                          new ScatterPlotStyle("", SymbolType.NONE, 0, Color.black));
//        axonsOnly.setLabels("Axon Speed", "Type Number");
//        axonsOnly.setLabelFont(labelFont);
//        axonsOnly.axesBorder.setRightPadding(rightPadding);
//        axonsOnly.axesBorder.setTopPadding(topPadding);
//        axonsOnly.setDoubleBuffered(false);
//        axonsOnly.setRange(.3, 1.3, 0, 33);
//        axonsOnly.saveAsPNG(new File(
//            "g:\\Publications\\thesis\\figures\\results\\axonsScatters.png"), 3, 3);



        p = new PlotPanel();
        p.addData(scats[2],
                  new ScatterPlotStyle("Scatter", SymbolType.FILLED_SQUARE, 3,
                                       Color.black));
        f = new FunctionSum(funcs[2]);
        p.addData(f, new FunctionStyle("f", Color.red, 1));
//       p.addBackgroundText("" + 'C', LEFT, TOP, letterFont, Color.black);
        p.setLabels("Relative Response Latency", "Type Number");
        p.setLabelFont(labelFont);
        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);
        p.setDoubleBuffered(false);
        p.setRange(.6, 2.0, 0, 33);
        jpanel.add(p);

        p = new PlotPanel();
        p.addData(scats[3],
                  new ScatterPlotStyle("Scatter", SymbolType.FILLED_SQUARE, 3,
                                       Color.black));
        f = new FunctionSum(funcs[3]);
        p.addData(f, new FunctionStyle("f", Color.red, 1));
//       p.addBackgroundText("" + 'D', LEFT, TOP, letterFont, Color.black);
        p.setLabels("Sustained Response Measure", "Type Number");
        p.setLabelFont(labelFont);
        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);
        p.setDoubleBuffered(false);
        p.setRange( -.4, .3, 0, 33);
        jpanel.add(p);

        p = new PlotPanel();
        p.addData(scats[4],
                  new ScatterPlotStyle("Scatter", SymbolType.FILLED_SQUARE, 3,
                                       Color.black));
        f = new FunctionSum(funcs[4]);
        p.addData(f, new FunctionStyle("f", Color.red, 1));
//       p.addBackgroundText("" + 'D', LEFT, TOP, letterFont, Color.black);

        scatter = new ScatterPlot();
        scatter.add(1.25, 0);
        scatter.add(1.25, 20.5);
        scatter.add(3, 20.5);
        scatter.add(3, 33);
        p.addData(scatter, divider);

        p.setLabels("Nonlinearity Index", "Type Number");
        p.setLabelFont(labelFont);
        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);
        p.setDoubleBuffered(false);
        p.setRange(0, 21, 0, 33);
        jpanel.add(p);

        p = new PlotPanel();
        p.addData(scats[5],
                  new ScatterPlotStyle("Scatter", SymbolType.FILLED_SQUARE, 3,
                                       Color.black));
        f = new FunctionSum(funcs[5]);
        p.addData(f, new FunctionStyle("f", Color.red, 1));
//       p.addBackgroundText("" + 'D', LEFT, TOP, letterFont, Color.black);

        scatter = new ScatterPlot();
        scatter.add(1.13, 0);
        scatter.add(1.13, 33);
        p.addData(scatter, divider);

        scatter = new ScatterPlot();
        scatter.add(.88, 0);
        scatter.add(.88, 33);
        p.addData(scatter, divider);

        p.setLabels("Object Motion Sensitivity Index", "Type Number");
        p.setLabelFont(labelFont);
        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);
        p.setDoubleBuffered(false);
        p.setRange(.3, 1.6, 0, 33);
        jpanel.add(p);

        JFrame frame = new JFrame();
        frame.setSize(72 * 6, (int) (72 * 7.5));
        frame.add(jpanel);
        frame.setVisible(true);

        RepaintManager.currentManager(frame).setDoubleBufferingEnabled(false);

        EpsGraphics2D g = new EpsGraphics2D(
            "", Font.decode("Arial ITALIC 8"), 0, 0, frame.getWidth(), frame.getHeight());
        g.scale(1, 1);
        g.setAccurateTextMode(true);
        frame.paint(g);

        FileWriter writer = new FileWriter(
            "g:\\Publications\\thesis\\figures\\results\\paramsScatters.eps");
        writer.write(g.toString());
        writer.flush();
        writer.close();
//      frame.setVisible(false);



//       PlotUtil.showData("", p);
        table2.draw(System.out);

        try {
            table2.excelDraw(new PrintStream(new File(
                "g:\\Publications\\thesis\\params.table")));
        } catch (IOException e) {
            e.printStackTrace();
        }

    }


    public static Num calculateNonlinearity(
        int[] id, ParametersFile pFile) throws
        IOException {
        Num maxNonLinearity = Num.zero;

        double[] _freqs = pFile.getArrayCell(id[0], "reversingFrequencies");
        final int nFreq = _freqs.length;
        for (int tPeriod = 0; tPeriod < 3; tPeriod++) {

            MeanVarianceCalculator[] mvcF1 = new MeanVarianceCalculator[nFreq];
            MeanVarianceCalculator[] mvcF2 = new MeanVarianceCalculator[nFreq];
            for (int periodIndex = 0; periodIndex < nFreq; periodIndex++) {
                mvcF1[periodIndex] = new MeanVarianceCalculator();
                mvcF2[periodIndex] = new MeanVarianceCalculator();
            }

            for (int i = 0; i < id.length; i++) {

                double[] _f1 = pFile.getArrayCell(id[
                                                 i], "T" + (tPeriod + 1) + "reversingF1");
                double[] _f2 = pFile.getArrayCell(id[
                                                 i], "T" + (tPeriod + 1) + "reversingF2");

                if (_f1 != null && _f2 != null) {
                    for (int j = 0; j < _f1.length; j++) {
                        mvcF1[j].add(_f1[j]);
                        mvcF2[j].add(_f2[j]);
                    }
                }
            }

            double[] f1 = new double[nFreq];
            double[] f2 = new double[nFreq];
            double[] f1Err = new double[nFreq];
            double[] f2Err = new double[nFreq];
            for (int sPeriod = 0; sPeriod < nFreq; sPeriod++) {
                f1[sPeriod] = mvcF1[sPeriod].getMean();
                f1Err[sPeriod] = mvcF1[sPeriod].getMeanVariance();
                f2[sPeriod] = mvcF2[sPeriod].getMean();
                f2Err[sPeriod] = mvcF1[sPeriod].getMeanVariance();

                //Calculate maximum Non-linearity index.
                Num NumF1 = new Num(f1[sPeriod], f1Err[sPeriod]);
                Num NumF2 = new Num(f2[sPeriod], f2Err[sPeriod]);
                Num ratio = NumF2.div(NumF1);
                if (ratio.x > maxNonLinearity.x) {
                    maxNonLinearity = ratio;
                }
            }

        }

        return maxNonLinearity;

    }


    public static void makeDetailedF1F2() throws Exception {
        double f1f2Width = 5;
        int n = 3;

        String[] pName = {
                         "g:\\data\\2005-05-02-0\\data002\\data002.params"
        };
        String[] name = {
                        "g:\\data\\2005-05-02-0\\data005"
        };
        String[] cName = {"All/Off/Large/T1/AOMS"};

        int binsPerPeriod = 30;

        // save the example response figures
        ReversingGratings m = new ReversingGratings(name[0], binsPerPeriod);

        int[] id = {9};
        int[] max = {600};
        int[] startPhase1 = {0};
        int[] startPhase3 = {0};
        int[] startPhase5 = {0};

        n = 0;
        GridLayout g = new GridLayout(3, 8);
        g.setHgap(1);
        g.setVgap(1);
        JPanel p = new JPanel(g);
        p.setBackground(Color.white);
        for (int i = 0; i < id.length; i++) {
            m.setCurrentNeuron(id[i]);
            p.removeAll();
            p.setBorder(new B(8, 10, 15, 30, "" + (char) ('B' + n)));

            for (int k = 0; k < 8; k++) {
                p.add(m.getPSTH(0, (startPhase1[i] + k) % 8, max[i]));
            }

            for (int k = 0; k < 8; k++) {
                p.add(m.getPSTH(3, (startPhase3[i] + k) % 8, max[i]));
            }

            for (int k = 0; k < 8; k++) {
                p.add(m.getPSTH(5, (startPhase5[i] + k) % 8, max[i]));
            }

            JFrame f = new JFrame();
            f.add(p);
            f.setVisible(true);
            GraphicsIO.saveAsEPS(p,
                "g:\\Publications\\thesis\\figures\\stimuli\\contrastReversingExample1.eps",
                                 f1f2Width, 2, false);
        }
    }


    static class B
        implements Border {

        int top;
        int left;
        int bottom;
        int right;
        String a;

        public B(int top, int left, int bottom, int right, String a) {
            this.top = top;
            this.left = left;
            this.bottom = bottom;
            this.right = right;
            this.a = a;
        }


        public Insets getBorderInsets(Component c) {
            return new Insets(top, left, bottom, right);
        }


        public boolean isBorderOpaque() {
            return false;
        }


        public void paintBorder(Component c, Graphics gIn, int x, int y, int width,
                                int height) {
            Graphics2D g = (Graphics2D) gIn;
            g.setColor(Color.black);
            g.setFont(letterFont);
            FontMetrics fm = g.getFontMetrics();
            g.setFont(Font.decode("Arial Italic 8"));

            int xx = 20 + left;
            for (int i = 0; i < 8; i++) {

                g.drawString("" + i * 45 + "\u00b0", xx,
                             fm.getAscent() - fm.getDescent() - 1);
                xx += 44;
            }

            xx -= 12;
            g.drawString(".17", xx + 5, 30);
            g.drawString("cyc/mm", xx, 30 + 8);

            g.drawString("1.3", xx + 5, 75);
            g.drawString("cyc/mm", xx, 75 + 8);

            g.drawString("5.4", xx + 5, 120);
            g.drawString("cyc/mm", xx, 120 + 8);

            int yPos = 52;
            g.translate(8, yPos);
            g.rotate( -Math.PI / 2);
            g.drawString("Spike Rate", 0, 0);
            g.rotate(Math.PI / 2);
            g.translate( -8, -yPos);
            yPos += 48;

            g.translate(8, yPos);
            g.rotate( -Math.PI / 2);
            g.drawString("Spike Rate", 0, 0);
            g.rotate(Math.PI / 2);
            g.translate( -8, -yPos);
            yPos += 48;

            g.translate(8, yPos);
            g.rotate( -Math.PI / 2);
            g.drawString("Spike Rate", 0, 0);
            g.rotate(Math.PI / 2);
            g.translate( -8, -yPos);

            g.drawLine(326, 153, 368, 153);
            g.drawString("125 ms", 333, 162);

        }
    }
//Test comment

    public static void main(String[] args) throws Exception {
        UIDefaults defaults = UIManager.getDefaults();
        
  //      NeuronViewer nv = new NeuronViewer("/Volumes/Canguro/data/2000-12-14-1/data051/data051.params", false, new Config("config.xml"));

        SwingUtilities.invokeAndWait(new Runnable() {
            public void run() {
                try {
//                    parameterFiles = new ParametersFile[paramFilesString.length];
//                    nvs = new NeuronViewer[paramFilesString.length];
//                    for (int i = 0; i < parameterFiles.length; i++) {
//                        parameterFiles[i] = new ParametersFile(paramFilesString[i]);
//                        nvs[i] = new NeuronViewer(paramFilesString[i], false);
//                    }
//
//                    whitenedParameterFiles = new ParametersFile[whitenedParamFilesString.
//                                             length];
//                    whitenedNVS = new NeuronViewer[whitenedParamFilesString.length];
//                    for (int i = 0; i < parameterFiles.length; i++) {
//                        whitenedParameterFiles[i] = new ParametersFile(
//                            whitenedParamFilesString[i]);
//                        whitenedNVS[i] = new NeuronViewer(whitenedParamFilesString[i], false);
//                    }
//
//                    referenceSimplifiedContourArea = new double[parameterFiles.length];
//                    referenceContourArea = new double[parameterFiles.length];
//                    referenceRL = new double[parameterFiles.length];
//                    referenceAxon = new double[parameterFiles.length];
//                    calculateReferenceParameters();
//                    calculateParameters(classList);

//                    makeClassFigures();
//                    makeMethodsFigures();
//                    makeResultsFigures();


//                    makeWhitenerFigure("/All/Off/Small/T3/A", "326");
//                    makeWhitenerFigure("/All/On/Large/T3/C", "66");

//                    makeDetailedF1F2();

                    /*                    makeSeparationQualityFigure(new String[] {
                              "G:\\data\\2005-05-02-0\\data002\\data002.params",
                     "G:\\data\\2005-05-02-0\\data002firsthalf\\data002firsthalf.params",
                              "G:\\data\\2005-05-02-0\\data001\\data001.params"}
                                                                    , 1);
                     */
                    /*                      makeSeparationQualityFigure(new String[] {
                                "G:\\data\\2005-05-04-0\\data002\\data002.params",
                     "G:\\data\\2005-05-04-0\\data002firsthalf\\data002firsthalf.params",
                                "G:\\data\\2005-05-04-0\\data001\\data001.params"}
                                                                      , 2);

                                          makeSeparationQualityFigure(new String[] {
                                "G:\\data\\2005-08-25-0\\data003\\data003.params",
                     "G:\\data\\2005-08-25-0\\data003firsthalf\\data003firsthalf.params",
                                "G:\\data\\2005-08-25-0\\data002\\data002.params"},
                                                                      3);
                     */
                } catch (Exception ex) {
                    ex.printStackTrace();
                }
            }
        });
//        System.exit(1);
    }
}
