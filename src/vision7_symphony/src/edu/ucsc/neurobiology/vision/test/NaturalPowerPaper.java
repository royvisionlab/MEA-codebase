package edu.ucsc.neurobiology.vision.test;

import java.io.*;
import java.util.*;

import java.awt.*;

import javax.swing.*;


import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.io.chunk.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.neuronviewer.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.util.*;
import edu.ucsc.neurobiology.vision.plot.eps.EpsGraphics2D;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreePath;
import edu.ucsc.neurobiology.vision.math.fitting.Gaussian1DFunction;
import edu.ucsc.neurobiology.vision.analysis.OMSClassification;
import edu.ucsc.neurobiology.vision.gui.PlaceLayout;
import edu.ucsc.neurobiology.vision.gui.PlaceC;
import javax.swing.plaf.ColorUIResource;
import edu.ucsc.neurobiology.vision.gui.GraphicsIO;
import edu.ucsc.neurobiology.vision.analysis.ReversingGratings;
import javax.swing.border.Border;
import static edu.ucsc.neurobiology.vision.plot.PlotPanel.*;
import edu.ucsc.neurobiology.vision.*;


/**
 * @author Matthew Grivich, University of California, Santa Cruz

 */
public class NaturalPowerPaper {


    public static final int THESIS = 0;
    public static final int PAPER = 1;
    public static final int POWERPOINT = 2;

    public static int mode = PAPER;


    static Font letterFont = Font.decode("Arial Bold 10");
    static Font textFont = Font.decode("Arial Italic 8");
    static Font labelFont = Font.decode("Arial Italic 8");

    static double[] referenceSimplifiedContourArea;
    static double[] referenceContourArea;
    static double[] referenceRL;
    static double[] referenceAxon;
    static double[] referenceAcfMean;
    static double[] referenceSpatialMean;
    static double[] referenceTemporalMean;
    static double[] referenceRLE;

    static NeuronViewer[] nvs, whitenedNVS;
    static ParametersFile[] parameterFiles, whitenedParameterFiles;
    static GlobalsFile[] globalsFiles;
    static String[] whitenedParamFilesString = {
        "e:\\data\\2005-05-02-0\\data002W\\data002W.params",
        "e:\\data\\2005-05-04-0\\data002W\\data002W.params",
    "e:\\data\\2005-08-25-0\\data003W\\data003W.params"};

    static String[] paramFilesString = {"e:\\data\\2005-05-02-0\\data002\\data002.params",
        "e:\\data\\2005-05-04-0\\data002\\data002.params",
    "e:\\data\\2005-08-25-0\\data003\\data003.params"};

    static String[] globalsFilesString = {"e:\\data\\2005-05-02-0\\data002\\data002.globals",
        "e:\\data\\2005-05-04-0\\data002\\data002.globals",
    "e:\\data\\2005-08-25-0\\data003\\data003.globals"};

    static String[] whitenerParamFilesString = {
        "e:\\data\\2005-08-25-0\\data003\\data003.params",
        "e:\\data\\2005-08-25-0\\data003W\\data003W.params",
    "e:\\data\\2005-08-25-0\\data002\\data002.params"};


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

    public static int[] prepList = {0, 2, 2, 0, 0, 0, 2, 0, 2, 1, 1, 2, 0, 2, 0, 0, 2, 0,
        0, 0, 2, 2, 0, 1, 0, 0, 2};
    




    public static String[] dsClassList = {

        "All/On/Large/T2/DS/1",
        "All/On/Large/T2/DS/3",
        "All/On/Large/T2/DS/2",
        "All/On/Large/T3/OS",
        "All/Off/Small/T2/DS/3",
        "All/Off/Small/T2/DS/1",
        "All/Off/Small/T2/DS/2",
        "All/Off/Small/T2/DS/4",
    };

    public static String[] allClassList = {
        "All/On/Small/T2/OMS",
        "All/On/Small/T3/OMS",
        "All/On/Large/T1/Triphasic",
        "All/On/Large/T2/DS/1",
        "All/On/Large/T2/DS/2",
        "All/On/Large/T2/DS/3",
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
        "All/Off/Small/T2/DS/1",
        "All/Off/Small/T2/DS/2",
        "All/Off/Small/T2/DS/3",
        "All/Off/Small/T2/DS/4",
        "All/Off/Small/T2/OMS",
        "All/Off/Small/T3/OMS",
        "All/Off/Small/T3/A",
        "All/Off/Small/T3/B",
        "All/Off/Small/T3/C",
        "All/Off/Large/T1/AOMS",
        "All/Off/Large/T1/A",
        "All/Off/Large/T3/A"

    };
    
    public static int[] allPrepList = {0, 2, 2,0,0,0, 0, 0, 2, 0, 2, 1, 1, 2, 0, 2, 0, 0, 2, 0,
        0,0,0,0, 0, 2, 2, 0, 1, 0, 0, 2};
    


    public static int[] dsPrepList = {0,0,0,1,0,0,0,0};

    private static String referenceType = "All/Off/Small/T2/OMS";






    public static void makeResultsFigures() throws Exception {
        double dpi = 72.0;
        int width = (int) (3.5 * (dpi));
        int height = (int) (4.5 * (dpi));

        int verticalAxisThickness = 0;//32;
        int leftPadding = 20;
        int rightPadding = 8;
        int topPadding = 4;
        int bottomPadding = 0;
        Font labelFont = Font.decode("Arial ITALIC 8");
        ScatterPlotStyle scatterPlotStyle = new ScatterPlotStyle(
                "Scatter Plot", SymbolType.DISK, 3, Color.black, false, Color.black, 1);

//		Make On DS plot
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
            double theta = 0;
            double xNew = x * Math.cos(theta) + y * Math.sin(theta);
            double yNew = - x * Math.sin(theta) + y * Math.cos(theta);
            sp.add(xNew, yNew);

        }
        p.addData(sp, scatterPlotStyle);


        p.setAntiallias(true);
        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);

        p.setLabelFont(labelFont);
        p.addBackgroundText("CP5", 60, 30, letterFont, Color.black);
        p.addBackgroundText("CP4", 130, 45, letterFont, Color.black);
        p.addBackgroundText("CP6", 80, 115, letterFont, Color.black);
        p.getXAxis().setLabel("Horizontal Sensitivity (+ is Anterior)");
        p.getXAxis().setManualTicks(new double[]{-.8, -.4, 0, .4, .8});
        p.getYAxis().setLabel("Vertical Sensitivity (+ is Superior)");
        p.getYAxis().setManualTicks(new double[]{-.8, -.4, 0, .4, .8});
        p.setRange(-.8, .8, -.8, .8);
        p.addBackgroundText("A", LEFT, TOP, letterFont, Color.black);

        JPanel panel = new JPanel(new PlaceLayout());
        JPanel dsCombinedPanel = new JPanel(new PlaceLayout());
        panel.add(p, new PlaceC(0, 0, 0, 0, 1, 2.0/6.0));

        whitenedNVS[0].showNeuronsInSubclasses = true;
        whitenedNVS[0].select("All/On/Large/T2/DS/1");


        p = (PlotPanel) makeMosaicPlot(whitenedNVS[0], leftPadding, rightPadding, topPadding, bottomPadding,
                verticalAxisThickness,
                width, height, labelFont, letterFont, "CP4", Color.black);
        p.setScaleBar(.500, 1, 0.05, 0.1);
        panel.add(p, new PlaceC(0, 2.0/6.0, 0, 0, 1, 1.0/6.0));

        whitenedNVS[0].select("All/On/Large/T2/DS/3");
        p = (PlotPanel) makeMosaicPlot(whitenedNVS[0], leftPadding, rightPadding, topPadding, bottomPadding,
                verticalAxisThickness,
                width, height, labelFont, letterFont, "CP5", Color.black);
        panel.add(p, new PlaceC(0, 3.0/6.0, 0, 0, 1, 1.0 / 6.0));

        whitenedNVS[0].select("All/On/Large/T2/DS/2");
        p = (PlotPanel) makeMosaicPlot(whitenedNVS[0], leftPadding, rightPadding, topPadding, bottomPadding,
                verticalAxisThickness,
                width, height, labelFont, letterFont, "CP6", Color.black);
        panel.add(p, new PlaceC(0, 4.0/6.0, 0, 0, 1, 1.0 / 6.0));
        dsCombinedPanel.add(panel, new PlaceC(0, 0, 0, 0, .32, 1));

        //Off ds plot
    /*	p = new PlotPanel("", false, false, false, false);

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
            double theta = 0;
            double xNew = x * Math.cos(theta) + y * Math.sin(theta);
            double yNew = -x * Math.sin(theta) + y * Math.cos(theta);
            sp.add(xNew, yNew);

        }
        p.addData(sp, scatterPlotStyle);


        p.setAntiallias(true);
        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);

        p.setLabelFont(labelFont);
        p.addBackgroundText("CP23", 40, 70, letterFont, Color.black);
        p.addBackgroundText("CP21", 123, 88, letterFont, Color.black);
        p.addBackgroundText("CP22", 113, 40, letterFont, Color.black);
        p.addBackgroundText("CP24", 88, 118, letterFont, Color.black);

        p.getXAxis().setLabel("Horizontal Sensitivity (+ is Anterior)");
        p.getXAxis().setManualTicks(new double[]{-.8, -.4, 0, .4, .8});
        p.getYAxis().setLabel("Vertical Sensitivity (+ is Superior)");
        p.getYAxis().setManualTicks(new double[]{-.8, -.4, 0, .4, .8});
        p.setRange( -.8, .8, -.8, .8);
        p.addBackgroundText("B", LEFT, TOP, letterFont, Color.black);

        panel = new JPanel(new PlaceLayout());
        panel.add(p, new PlaceC(0, 0, 0, 0, 1, 2.0/6.0));

        whitenedNVS[0].select("All/Off/Small/T2/DS/3");
        p = (PlotPanel) makeMosaicPlot(whitenedNVS[0], leftPadding, rightPadding, topPadding, bottomPadding,
                verticalAxisThickness,
                width, height, labelFont, letterFont, "CP22", Color.black);
        panel.add(p, new PlaceC(0, 2.0/6.0, 0, 0, 1, 1.0 / 6.0));

        whitenedNVS[0].select("All/Off/Small/T2/DS/1");
        p = (PlotPanel) makeMosaicPlot(whitenedNVS[0], leftPadding, rightPadding, topPadding, bottomPadding,
                verticalAxisThickness,
                width, height, labelFont, letterFont, "CP21", Color.black);
        panel.add(p, new PlaceC(0, 3.0/6.0, 0, 0, 1, 1.0 / 6.0));

        whitenedNVS[0].select("All/Off/Small/T2/DS/2");
        p = (PlotPanel) makeMosaicPlot(whitenedNVS[0], leftPadding, rightPadding, topPadding, bottomPadding,
                verticalAxisThickness,
                width, height, labelFont, letterFont, "CP23", Color.black);
        panel.add(p, new PlaceC(0, 4.0/6.0, 0, 0, 1, 1.0 / 6.0));

        whitenedNVS[0].select("All/Off/Small/T2/DS/4");
        p = (PlotPanel) makeMosaicPlot(whitenedNVS[0], leftPadding, rightPadding, topPadding, bottomPadding,
                verticalAxisThickness,
                width, height, labelFont, letterFont, "CP24", Color.black);
        panel.add(p, new PlaceC(0, 5.0/6.0, 0, 0, 1, 1.0 / 6.0));
        dsCombinedPanel.add(panel, new PlaceC(.34, 0, 0, 0, .32, 1));
*/


//		Make OS plot
/*
        nvs[1].showNeuronsInSubclasses = false;
        nvs[1].select("/All/On/Large/T3/OS");
        p = (PlotPanel) nvs[1].getPlot("overlayedSinusoids");

        p.setLabelFont(labelFont);
        p.setRange(p.getRange()[0], p.getRange()[1], 0, p.getRange()[3]);
        p.setDoubleBuffered(false);
        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);
        p.addBackgroundText("C", LEFT, TOP, letterFont, Color.black);

        panel = new JPanel(new PlaceLayout());
        panel.add(p, new PlaceC(0, 0, 0, 0, 1, 2.0/6.0));

        whitenedNVS[1].select("/All/On/Large/T3/OS");
        p = (PlotPanel) makeMosaicPlot(whitenedNVS[1], leftPadding, rightPadding, topPadding, bottomPadding,
                verticalAxisThickness,
                width, height, labelFont, letterFont, "CP13", Color.black);
        panel.add(p, new PlaceC(0, 2.0/6.0, 0, 0, 1, 1.0 / 6.0));

        dsCombinedPanel.add(panel, new PlaceC(.67, 0, 0, 0, .33, 1));

        JFrame frame = new JFrame();

        frame.add(dsCombinedPanel);
        frame.setSize( (int) (7.24 * 72), (int) (7.0 * 72));
        frame.setVisible(true);
        RepaintManager.currentManager(frame).setDoubleBufferingEnabled(false);
        EpsGraphics2D g = new EpsGraphics2D(
                "", Font.decode("Arial ITALIC 8"), 0, 0, frame.getWidth(), frame.getHeight());
        g.setColorDepth(g.RGB);

        g.scale(72.0 / dpi, 72.0 / dpi);
        g.setAccurateTextMode(true);
        frame.paint(g);



        FileWriter writer = new FileWriter(
        "e:\\Publications\\GuineaPaper\\figures\\combinedDS.eps");

        writer.write(g.toString());
        writer.flush();
        writer.close();
        RepaintManager.currentManager(frame).setDoubleBufferingEnabled(true);
        frame = new JFrame();

*/



        /*
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
        p.addBackgroundText("A", LEFT, TOP, letterFont, Color.red);

        //        PlotUtil.showData("", p);
        //        p.saveAsEPS(new File(
        //            "e:\\Publications\\GuineaPaper\\figures\\blueness.eps"), 2, 2, true);

        panel = new JPanel(new PlaceLayout());

        panel.add(p, new PlaceC(0, 0, 0, 0, .33, .5));
        nvs[0].showNeuronsInSubclasses = false;
        nvs[1].showNeuronsInSubclasses = false;
        nvs[2].showNeuronsInSubclasses = false;

        //X-cell
        nvs[2].select("All/On/Small/T3/OMS");
        panel.add(makeF1F2Plot(nvs[2], rightPadding, topPadding,
                               verticalAxisThickness,
                               width, height, labelFont, letterFont, "B"),
                  new PlaceC(.33, 0, 0, 0, .33, .5));

        //Y-cell
        nvs[0].select("All/Off/Large/T1/AOMS");
        panel.add(makeF1F2Plot(nvs[0], rightPadding, topPadding,
                               verticalAxisThickness,
                               width, height, labelFont, letterFont, "C"),
                  new PlaceC(.67, 0, 0, 0, .33, .5));

        //non-OMS
        nvs[2].select("All/On/Large/T1/Triphasic");
       panel.add(makeOMSPlot(nvs[2], rightPadding, topPadding,
                              verticalAxisThickness,
                              width, height, labelFont, letterFont, "D"),
                 new PlaceC(0, .5, 0, 0, .33, .5));

       //OMS
       nvs[0].select("All/On/Small/T2/OMS");
      panel.add(makeOMSPlot(nvs[0], rightPadding, topPadding,
                             verticalAxisThickness,
                             width, height, labelFont, letterFont, "E"),
                new PlaceC(.33, .5, 0, 0, .33, .5));

      //anti-OMS
       nvs[0].select("All/On/Large/T2/AOMS");
      panel.add(makeOMSPlot(nvs[0], rightPadding, topPadding,
                             verticalAxisThickness,
                             width, height, labelFont, letterFont, "F"),
                new PlaceC(.67, .5, 0, 0, .33, .5));




        frame = new JFrame();
        frame.add(panel);
        frame.setSize( (int) (7.24 * 72), (int) (3.5 * 72));
        frame.setVisible(true);

        RepaintManager.currentManager(frame).setDoubleBufferingEnabled(false);
        g = new EpsGraphics2D(
            "", Font.decode("Arial ITALIC 8"), 0, 0, frame.getWidth(), frame.getHeight());
            g.setColorDepth(g.RGB);
        g.scale(72.0 / dpi, 72.0 / dpi);

        g.setAccurateTextMode(true);
        frame.paint(g);



        writer = new FileWriter(
            "e:\\Publications\\GuineaPaper\\figures\\CombinedResults.eps");

        writer.write(g.toString());
        writer.flush();
        writer.close();
        RepaintManager.currentManager(frame).setDoubleBufferingEnabled(true);
         */
        
        panel = new JPanel(new PlaceLayout());

verticalAxisThickness=42;

//Monophasic STA
 nvs[0].select("All/On/Large/T3/F");
 panel.add(makeTCPlot(nvs[0], rightPadding, topPadding,
                        verticalAxisThickness,
                        width, height, labelFont, letterFont, "A"),
           new PlaceC(0, 0, 0, 0, .33, .33));
 //Biphasic STA
 nvs[0].select("All/On/Large/T3/C");
 panel.add(makeTCPlot(nvs[0], rightPadding, topPadding,
                        verticalAxisThickness,
                        width, height, labelFont, letterFont, "B"),
           new PlaceC(.33, 0, 0, 0, .33, .33));

 //Triphasic STA
 nvs[2].select("All/On/Large/T1/Triphasic");
 panel.add(makeTCPlot(nvs[2], rightPadding, topPadding,
                        verticalAxisThickness,
                        width, height, labelFont, letterFont, "C"),
           new PlaceC(.67, 0, 0, 0, .33, .33));

 //Monophasic STV
 nvs[0].select("All/On/Large/T3/D");
panel.add(makeTCVPlot(nvs[0], rightPadding, topPadding,
                       verticalAxisThickness,
                       width, height, labelFont, letterFont, "D"),
          new PlaceC(0, .33, 0, 0, .33, .33));

//Biphasic STV
nvs[0].select("All/On/Large/T3/C");
panel.add(makeTCVPlot(nvs[0], rightPadding, topPadding,
                      verticalAxisThickness,
                      width, height, labelFont, letterFont, "E"),
         new PlaceC(.33, .33, 0, 0, .33, .33));

//Double Peaked STV
nvs[0].select("All/On/Large/T2/AOMS");
panel.add(makeTCVPlot(nvs[0], rightPadding, topPadding,
                      verticalAxisThickness,
                      width, height, labelFont, letterFont, "F"),
         new PlaceC(.67, .33, 0, 0, .33, .33));

//Uniformity STV
nvs[0].select("All/On/Large/T3/UniformityA");
panel.add(makeTCVPlot(nvs[0], rightPadding, topPadding,
                      verticalAxisThickness,
                      width, height, labelFont, letterFont, "G"),
         new PlaceC(0, .67, 0, 0, .33, .33));




 JFrame frame = new JFrame();
 frame.add(panel);
 frame.setSize( (int) (7.24 * 72), (int) (5.25 * 72));
 frame.setVisible(true);

 RepaintManager.currentManager(frame).setDoubleBufferingEnabled(false);
 EpsGraphics2D g = new EpsGraphics2D(
     "", Font.decode("Arial ITALIC 8"), 0, 0, frame.getWidth(), frame.getHeight());
     g.setColorDepth(g.RGB);
 g.scale(72.0 / dpi, 72.0 / dpi);

 g.setAccurateTextMode(true);
 frame.paint(g);



 FileWriter writer = new FileWriter(
     "e:\\Publications\\GuineaPaper\\figures\\TimeCourses.eps");

 writer.write(g.toString());
 writer.flush();
 writer.close();
 RepaintManager.currentManager(frame).setDoubleBufferingEnabled(true);
         

    }


    public static void makeMethodsFigures() throws Exception {
        NeuronViewer nv = new NeuronViewer(
                "e:\\data\\2005-05-02-0\\data002\\data002.params", false, new Config("config.xml"));

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
                    "e:\\data\\2005-08-25-0\\data002\\data002.params", false);
                ParametersFile pf3 = new ParametersFile(
                    "e:\\data\\2005-08-25-0\\data002\\data002.params");
                double x0 = ( (Double) pf3.evaluate("x0", new Integer(11))) *
                            (5.8 * 2) * 5;
                double y0 = ( (Double) pf3.evaluate("y0", new Integer(11))) *
                            (5.8 * 2) * 5;
//nv3.select("All/On/Y/T2/B/34");
                nv3.select("All/Off/Small/T2/OMS/11");
         PlotPanel panel = makeSTAPlot(nv3, 0, 0, 0, 3 * 72, (int) (1.8 * 72), labelFont,
                                              letterFont, "", x0, y0, .800);
                panel.setAxisVisible(false);
                panel.setScaleBar(500, 2, 0.05, 0.1);
                Frame f = new Frame();
                frame.setSize(3 * 72 + 25, (int) (1.8 * 72 + 40));
                frame.add(panel);
         frame.setVisible(true); //Click save as eps on this.  Next line does not work
//panel.saveAsEPS(new File("e:\\publications\\thesis\\Figures\\stimuli\\sampleSTA.eps"), 4*72, 2*72, false);

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
                   "e:\\Publications\\thesis\\figures\\stimuli\\flashesExample.eps");
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
         "e:\\Publications\\thesis\\figures\\stimuli\\sinExample1.eps"), 6, 2, false);


                         p = (PlotPanel) nv.getPlot("AverageGrating");
                           p.setLabelFont(labelFont);
                         p.axesBorder.setRightPadding(rightPadding);
                 p.axesBorder.setTopPadding(topPadding);
                 p.setDoubleBuffered(false);


                      p.saveAsEPS(new File(
         "e:\\Publications\\thesis\\figures\\stimuli\\sinExample2.eps"), 6, 2, false);

                         p = (PlotPanel) nv.getPlot("AverageGratingPanel2D");
                         p.setLabelFont(labelFont);
                         p.axesBorder.setRightPadding(rightPadding);
                p.axesBorder.setTopPadding(topPadding);
                p.setDoubleBuffered(false);


                      p.saveAsEPS(new File(
         "e:\\Publications\\thesis\\figures\\stimuli\\sinExample3.eps"), 6, 3, false);
         */

        nv.select("/All/On/Small/T2/OMS/174");
        PlotPanel p = (PlotPanel) nv.getPlot("omsScatter");
        p.setDoubleBuffered(false);
        p.axesBorder.setRightPadding(rightPadding);
        p.setLabels("Distance(\u03bcm)", "Spikes/Run");
        p.setLabelFont(labelFont);
        p.setRange(0, 1000, 0, p.getRange()[3] + 10);

        p.saveAsEPS(new File(
        "e:\\Publications\\thesis\\figures\\stimuli\\omsExample1.eps"), 4, 2, false);
        p.setLabels("Distance(\u03bcm)", "Spikes/Run");

        p = (PlotPanel) nv.getPlot("omsHist");
        p.setDoubleBuffered(false);
        p.axesBorder.setRightPadding(rightPadding);

        p.setLabelFont(labelFont);
        p.setRange(0, 1000, 0, p.getRange()[3] + 10);
        p.setLabels("Distance(\u03bcm)", "Spikes/Run");
        p.saveAsEPS(new File(
        "e:\\Publications\\thesis\\figures\\stimuli\\omsExample2.eps"), 4, 2, false);

        /*         nv.select("/All/On/TwoPeaks/XOMSA");
                p = (PlotPanel) nv.getPlot("omsHistograms");
                p.setLabels("Distance(\u03bcm)", "Normalized Spike Rate");
                      p.setLabelFont(labelFont);
                        p.setRange(0,1000, 0, p.getRange()[3] );

                 p.saveAsEPS(new File(
         "e:\\Publications\\thesis\\figures\\stimuli\\omsExample3.eps"), 4, 2, false);*/


        nv.select("All/Off/Large/T1/AOMS/9");
        p = (PlotPanel) nv.getPlot("F1F2");
        p.autoscale();
        p.setRange(p.getRange()[0], p.getRange()[1], 0, p.getRange()[3] + 10);
        p.setLabelFont(labelFont);
        p.axesBorder.getVerticalAxis().setThickness(verticalAxisThickness);
        p.axesBorder.getVerticalAxis().setLabelPosition(1);

        p.saveAsEPS(new File(
        "e:\\Publications\\thesis\\figures\\stimuli\\contrastReversingExample.eps"),
        4, 2, false);

        //Sample white STA timecourse
        /*           nv = new NeuronViewer("e:\\data\\2005-08-25-0\\data002\\data002.params", false);
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
         "e:\\Publications\\thesis\\figures\\stimuli\\whiteExampleSTA.eps"), 4, 2, false);
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
         "e:\\Publications\\thesis\\figures\\stimuli\\whiteExampleSTV.eps"), 4, 2, false);
         */
    }


    public static PlotPanel makeMosaicPlot(NeuronViewer nv, int leftPadding, int rightPadding,
            int topPadding, int bottomPadding,
            int verticalAxisThickness,
            int width, int height, Font labelFont,
            Font letterFont, String letter, Color letterColor) {

        PlotPanel p = (PlotPanel) nv.getPlot("mosaic");
        p.setDoubleBuffered(false);
        p.setLegendVisible(false);
        p.addBackgroundText(letter, LEFT, TOP, letterFont, letterColor);
        if(verticalAxisThickness == 0) {
            p.setAxisVisible(false);
        } else {
            p.setAxisVisible(true);
        }
        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);
        p.axesBorder.setBottomPadding(bottomPadding);
        p.axesBorder.setLeftPadding(leftPadding);
        p.axesBorder.getVerticalAxis().setThickness(verticalAxisThickness);
        p.axesBorder.getVerticalAxis().setLabelPosition(1);
        double[] range = p.getRange();
        p.axesBorder.getHorizontalAxis().setManualTicks(new double[] {range[0], (range[0] + range[1])/2, range[1]});
        p.axesBorder.getVerticalAxis().setManualTicks(new double[] {range[2], range[3]});
        p.setPreferredSize(new Dimension(width, height));
        p.setLabelFont(labelFont);
        p.clearLegend();
        //      p.getXAxis().setFixedTickSpacing(fixedMinimumTick, fixedTickSpacing, fixedMaximumFractionarDigits)
        return p;

    }


    public static PlotPanel makeSTAPlot(NeuronViewer nv,
            int rightPadding, int topPadding,
            int verticalAxisThickness, int width, int height,
            Font labelFont, Font letterFont, String letter, Color letterColor,
            double x0, double y0, double viewRadius) {

        PlotPanel p = (PlotPanel) nv.getPlot("sta");
        p.setDoubleBuffered(false);
        p.addBackgroundText(letter, LEFT, TOP, letterFont, letterColor);
        p.setRange(x0 - viewRadius, x0 + viewRadius, y0 - viewRadius / 2,
                y0 + viewRadius / 2);
        p.setAxisVisible(true);
        p.setLabels("x (mm)", "y (mm)");
        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);
        p.axesBorder.getVerticalAxis().setThickness(verticalAxisThickness);
        p.axesBorder.getVerticalAxis().setLabelPosition(1);
        double[] range = p.getRange();
        p.axesBorder.getHorizontalAxis().setManualTicks(new double[] {range[0], (range[0] + range[1])/2, range[1]});
        p.axesBorder.getVerticalAxis().setManualTicks(new double[] {range[2], range[3]});
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
        double[] range = p.getRange();
        range[0] = -660;
        range[1] = 0;
        p.axesBorder.getHorizontalAxis().setManualTicks(new double[] {range[0], (range[0] + range[1])/2, range[1]});
        p.axesBorder.getVerticalAxis().setManualTicks(new double[] {range[2], 0, range[3]});
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
        double[] range = p.getRange();
        range[0] = -660;
        range[1] = 0;
        p.setRange(range);
        p.axesBorder.getHorizontalAxis().setManualTicks(new double[] {range[0], (range[0] + range[1])/2, range[1]});
        p.axesBorder.getVerticalAxis().setManualTicks(new double[] {range[2], range[3]});
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
        double[] range = p.getRange();
        range[0] = 0;
        range[1] = 100;
        p.setRange(range);
        p.axesBorder.getHorizontalAxis().setManualTicks(new double[] {range[0], (range[0] + range[1])/2, range[1]});
        p.axesBorder.getVerticalAxis().setManualTicks(new double[] {range[2], range[3]});
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
        p.setRange(p.getRange()[0]-.02, p.getRange()[1]+1, 0,
                p.getRange()[3] + p.getRange()[3] * .2);
        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);
        p.axesBorder.getVerticalAxis().setThickness(verticalAxisThickness);
        p.axesBorder.getVerticalAxis().setLabelPosition(1);
        double[] range = p.getRange();
        p.axesBorder.getVerticalAxis().setManualTicks(new double[] {range[2], range[3]});
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
        double[] range = p.getRange();
        range[0] = 0;
        range[1] = 850;
        range[2] = 0;
        range[3]*= 1.2;
        p.setRange(range);
        p.axesBorder.getHorizontalAxis().setManualTicks(new double[] {range[0], (range[0] + range[1])/2,range[1]});
        p.axesBorder.getVerticalAxis().setManualTicks(new double[] {range[2],  range[3]});
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
        double[] range = p.getRange();
        p.axesBorder.getHorizontalAxis().setManualTicks(new double[] {range[0], (range[0] + range[1])/2, range[1]});
        p.axesBorder.getVerticalAxis().setManualTicks(new double[] {range[2], range[3]});
        p.setPreferredSize(new Dimension(width, height));
        p.setLabelFont(labelFont);
        return p;
    }
    
    public static PlotPanel makeACFWPlot(NeuronViewer nv, int rightPadding,
            int topPadding,
            int verticalAxisThickness, int width,
            int height, Font labelFont, Font letterFont,
            String letter, int[] neurons) {
        PlotPanel p = null;
        try {
        p = NeuronViewer.makeHistogram(nv, "acfMean",
                2, 0, 100,
                null, neurons);
        } catch (CannotEvaluateException e) {
            e.printStackTrace();
        }

        p.setDoubleBuffered(false);
        p.setLegendVisible(false);
        p.setLabels("Time(ms)", "Count");
        p.addBackgroundText(letter, LEFT, TOP, letterFont, Color.red);
        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);
        p.axesBorder.getVerticalAxis().setThickness(verticalAxisThickness);
        p.axesBorder.getVerticalAxis().setLabelPosition(1);
        double[] range = p.getRange();
        p.axesBorder.getHorizontalAxis().setManualTicks(new double[] {range[0], (range[0] + range[1])/2, range[1]});
        p.axesBorder.getVerticalAxis().setManualTicks(new double[] {range[2], range[3]});
        p.setPreferredSize(new Dimension(width, height));
        p.setLabelFont(labelFont);
        return p;
    }
    
    public static PlotPanel makeSizePlot(NeuronViewer nv, int rightPadding,
            int topPadding,
            int verticalAxisThickness, int width,
            int height, Font labelFont, Font letterFont,
            String letter, int[] neurons) {
        PlotPanel p = null;
        try {
        p = NeuronViewer.makeHistogram(nv, "5.8*2*((simpleContourArea/3.14159)^.5)",
                2, 0, 120,
                null, neurons);
        } catch (CannotEvaluateException e) {
            e.printStackTrace();
        }

        p.setDoubleBuffered(false);
        p.setLegendVisible(false);
        p.setLabels("Simplified Countour Diameter (microns)", "Count");
        p.addBackgroundText(letter, LEFT, TOP, letterFont, Color.red);
        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);
        p.axesBorder.getVerticalAxis().setThickness(verticalAxisThickness);
        p.axesBorder.getVerticalAxis().setLabelPosition(1);
        double[] range = p.getRange();
        p.axesBorder.getHorizontalAxis().setManualTicks(new double[] {range[0], (range[0] + range[1])/2, range[1]});
        p.axesBorder.getVerticalAxis().setManualTicks(new double[] {range[2], range[3]});
        p.setPreferredSize(new Dimension(width, height));
        p.setLabelFont(labelFont);
        return p;
    }
    
    public static PlotPanel makeSpatTempPlot(NeuronViewer nv, int rightPadding,
            int topPadding,
            int verticalAxisThickness, int width,
            int height, Font labelFont, Font letterFont,
            String letter, int[] neurons) {
        PlotPanel p = null;
        try {
        p = nv.makeScatterPlot("spatialMean","temporalMean", "spatial mean", "temporal mean",
                null, neurons);
         p.setRange(-.77, 1.63, -.63, 1.48);
        } catch (CannotEvaluateException e) {
            e.printStackTrace();
        }

        p.setDoubleBuffered(false);
        p.setLegendVisible(false);
        p.setLabels("spatial mean", "temporal mean");
        p.addBackgroundText(letter, LEFT, TOP, letterFont, Color.red);
        p.axesBorder.setRightPadding(rightPadding);
        p.axesBorder.setTopPadding(topPadding);
        p.axesBorder.getVerticalAxis().setThickness(verticalAxisThickness);
        p.axesBorder.getVerticalAxis().setLabelPosition(1);
        double[] range = p.getRange();
        p.axesBorder.getHorizontalAxis().setManualTicks(new double[] {range[0], (range[0] + range[1])/2, range[1]});
        p.axesBorder.getVerticalAxis().setManualTicks(new double[] {range[2], range[3]});
        p.setPreferredSize(new Dimension(width, height));
        p.setLabelFont(labelFont);
        return p;
    }

    public static void makeClassFigure(int nv, String className, String classNumber) throws Exception {
        String neuronToShow = null;
        int[] neuronsInClass = parameterFiles[nv].getNeuronsInClass(className);
        if(neuronsInClass != null && neuronsInClass.length > 0) {


            neuronToShow = new Integer(neuronsInClass[0]).toString();

            makeClassFigure(nv, className, classNumber, neuronToShow);
        }

    }

    public static void makeClassFigure(int nv, String className, String classNumber,
            String neuronToShow) throws
            Exception {

        JFrame frame = new JFrame();

        frame.setLayout(new GridBagLayout());
        GridBagConstraints c = new GridBagConstraints();

        double dpi = 72.0;
        int width = (int) (2.41 * (dpi));
        int height = (int) (1.360 /*1.7*/ * (dpi));
        int staHeight = (int) (1.360 /*1.6*/ * (dpi));
        int verticalAxisThickness = 43;
        int leftPadding = 0;
        int rightPadding = 15;
        int topPadding = 4;
        int bottomPadding = 0;
        Font labelFont = Font.decode("Arial ITALIC 8");
        //Font letterFont = Font.decode("Arial Bold 10");
        Font letterFont = Font.decode("Arial Bold 7");
        frame.setSize(width * 3 + 8, height * 5 + 34);

        whitenedNVS[nv].showNeuronsInSubclasses = true;
        whitenedNVS[nv].select(className);
        nvs[nv].showNeuronsInSubclasses = true;
        nvs[nv].select(className);
        c.gridx = 0;
        c.gridy = 0;

//		frame.add(makeMosaicPlot(nvs[nv], leftPadding, rightPadding, topPadding, bottomPadding,
//				verticalAxisThickness,
//				width, staHeight, labelFont, letterFont, "Piece: " + (nv + 1) + " Class: " + classNumber, Color.red), c);

        frame.add(makeMosaicPlot(nvs[nv], leftPadding, rightPadding, topPadding, bottomPadding,
        verticalAxisThickness,
        width, staHeight, labelFont, letterFont, "A", Color.red), c);

        whitenedNVS[nv].showNeuronsInSubclasses = false;

        whitenedNVS[nv].select(className + "/" + neuronToShow);

        nvs[nv].showNeuronsInSubclasses = false;

        nvs[nv].select(className + "/" + neuronToShow);

        c.gridx = 1;
        c.gridy = 0;
        GlobalsFile.RunTimeMovieParams runParams = globalsFiles[nv].getRunTimeMovieParams();
        double x0 = (( (Double) parameterFiles[nv].evaluate("x0", new Integer(neuronToShow))) +
                runParams.xOffset/runParams.micronsPerStixelX)*runParams.micronsPerStixelX/1000;
        double y0 = (( (Double) parameterFiles[nv].evaluate("y0", new Integer(neuronToShow))) +
                runParams.yOffset/runParams.micronsPerStixelY)*runParams.micronsPerStixelY/1000;

        frame.add(makeSTAPlot(whitenedNVS[nv], rightPadding, topPadding,
                verticalAxisThickness,
                width, staHeight, labelFont, letterFont, "B", Color.white, x0, y0, .8),
                c);

//		whitenedNVS[nv].select(className);
//
//		c.gridx = 2;
//		c.gridy = 0;
//		frame.add(makeTCPlot(whitenedNVS[nv], rightPadding, topPadding,
//				verticalAxisThickness,
//				width, height, labelFont, letterFont, "C"), c);

        nvs[nv].select(className);

        c.gridx = 2;
        c.gridy = 0;
        frame.add(makeTCPlot(nvs[nv], rightPadding, topPadding,
                verticalAxisThickness,
                width, height, labelFont, letterFont, "C"), c);

        c.gridx = 0;
        c.gridy = 1;
        frame.add(makeTCVPlot(nvs[nv], rightPadding, topPadding,
                verticalAxisThickness,
                width, height, labelFont, letterFont, "D"), c);

        c.gridx = 1;
        c.gridy = 1;
        frame.add(makeAverageGrating2DPlot(nvs[nv], rightPadding,
                topPadding, verticalAxisThickness, width,
                height, labelFont, letterFont, "E"), c);

        c.gridx = 2;
        c.gridy = 1;
        frame.add(makeACFPlot(nvs[nv], rightPadding, topPadding,
                verticalAxisThickness,
                width, height, labelFont, letterFont, "F"), c);

        c.gridx = 0;
        c.gridy = 2;
        frame.add(makeF1F2Plot(nvs[nv], rightPadding, topPadding,
                verticalAxisThickness,
                width, height, labelFont, letterFont, "G"), c);

        c.gridx = 1;
        c.gridy = 2;
        frame.add(makeOMSPlot(nvs[nv], rightPadding, topPadding,
                verticalAxisThickness,
                width, height, labelFont, letterFont, "H"), c);



        c.gridx = 2;
        c.gridy = 2;
        frame.add(makeFlashesPlot(nvs[nv], rightPadding, topPadding,
                verticalAxisThickness,
                width, height, labelFont, letterFont, "I"), c);
        
//		c.gridx = 1;
//		c.gridy = 3;
//		frame.add(makeACFWPlot(nvs[nv], rightPadding, topPadding,
//				verticalAxisThickness,
//				width, height, labelFont, letterFont, "K", parameterFiles[nv].getNeuronsInClass(className))
//				, c);
//		
//		c.gridx = 2;
//		c.gridy = 3;
//		frame.add(makeSpatTempPlot(nvs[nv], rightPadding, topPadding,
//				verticalAxisThickness,
//				width, height, labelFont, letterFont, "L", parameterFiles[nv].getNeuronsInClass(className))
//				, c);
//		
//		c.gridx = 0;
//		c.gridy = 4;
//		frame.add(makeSizePlot(nvs[nv], rightPadding, topPadding,
//				verticalAxisThickness,
//				width, height, labelFont, letterFont, "M", parameterFiles[nv].getNeuronsInClass(className))
//				, c);



        frame.setVisible(true);

        RepaintManager.currentManager(frame).setDoubleBufferingEnabled(false);
        EpsGraphics2D g = new EpsGraphics2D(
                "", Font.decode("Arial ITALIC 8"), 0, 0, frame.getWidth(), frame.getHeight());

        g.scale(72.0 / dpi, 72.0 / dpi);
        g.setAccurateTextMode(true);
        g.setColorDepth(g.RGB);
        frame.paint(g);

        RepaintManager.currentManager(frame).setDoubleBufferingEnabled(true);

        FileWriter writer = new FileWriter(
                "e:\\Publications\\GuineaPaper\\figures\\classes\\" + (nv + 1) +className.replace("/", "") +
        ".eps");

        writer.write(g.toString());
        writer.flush();
        writer.close();

        frame.setVisible(false);
        frame.dispose();

    }




    public static void makeWhitenerFigure(String className, String neuronToShow) throws
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
        int width = (int) (2.33 * (dpi));
        int height = (int) (2.5 / 2.0 /*1.7*/ * (dpi));
//		int staHeight = (int) (1.2 * (dpi));
        int verticalAxisThickness = 39;
        int leftPadding = 0;
        int rightPadding = 9;
        int topPadding = 6;
        int bottomPadding = 0;
        Font labelFont = Font.decode("Arial ITALIC 8");
        Font letterFont = Font.decode("Arial BOLD 10");
        frame.setSize(width * 3 + 8, height * 3 + 34);

        //Make mosaic plot
        wnvs[0].showNeuronsInSubclasses = false;
        wnvs[0].select(className);

        c.gridx = 0;
        c.gridy = 0;
        frame.add(makeMosaicPlot(wnvs[0], leftPadding, rightPadding, topPadding, bottomPadding,
                verticalAxisThickness,
                width, height, labelFont, letterFont, "A", Color.red)
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
                width, height, labelFont, letterFont, "B", Color.white, x0, y0, .9),
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
        frame.add(makeMosaicPlot(wnvs[1], leftPadding, rightPadding, topPadding, bottomPadding,
                verticalAxisThickness,
                width, height, labelFont, letterFont, "D", Color.red)
                , c);

        wnvs[1].showNeuronsInSubclasses = false;

        //Make STA plot


        wnvs[1].select(className + "/" + neuronToShow);

        c.gridx = 1;
        c.gridy = 1;
        frame.add(makeSTAPlot(wnvs[1], rightPadding, topPadding,
                verticalAxisThickness,
                width, height, labelFont, letterFont, "E", Color.white, x0, y0, .9),
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
        frame.add(makeMosaicPlot(wnvs[2], leftPadding, rightPadding, topPadding, bottomPadding,
                verticalAxisThickness,
                width, height, labelFont, letterFont, "G", Color.red)
                , c);

        wnvs[2].showNeuronsInSubclasses = false;

        //Make STA plot

        wnvs[2].select(className + "/" + neuronToShow);
        c.gridx = 1;
        c.gridy = 2;
        frame.add(makeSTAPlot(wnvs[2], rightPadding, topPadding,
                verticalAxisThickness,
                width, height, labelFont, letterFont, "H", Color.white, x0, y0, .900),
                c);

        wnvs[2].select(className);
        c.gridx = 2;
        c.gridy = 2;
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
                "e:\\Publications\\GuineaPaper\\figures\\" + className.replace("/", "") +
        ".eps");

        writer.write(g.toString());
        writer.flush();
        writer.close();
//		frame.setVisible(true);
//		frame.dispose();

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
//		p.axesBorder.getVerticalAxis().setThickness(verticalAxisThickness);
//		p.axesBorder.getVerticalAxis().setLabelPosition(1);
        p.setPreferredSize(new Dimension(width, height));
        p.setLabelFont(labelFont);
        p.autoscale();

        p.setRange(0, 10, .5 /*p.getRange()[2]*/, 1.1);

        PlotUtil.showData("", p);

        p.saveAsEPS(new File(
                "e:\\Publications\\thesis\\figures\\results\\separation" +
                preparationNumber + ".eps"), 3, 3, false);

    }


    public static void makeDSMosaicsPlot() {
        JFrame frame = new JFrame();

        double dpi = 72.0;
        int width = (int) (3.2 * (dpi));
        int height = (int) (2.3 /*1.7*/ * (dpi));
        int staHeight = (int) (1.6 /*1.6*/ * (dpi));
        int verticalAxisThickness = 0;
        int rightPadding = 2;
        int topPadding = 1;
        int bottomPadding = 1;
        int leftPadding = 0;
        Font labelFont = Font.decode("Arial ITALIC 10");
        Font letterFont = Font.decode("Arial Bold 10");
        frame.setSize(7 * 72, 7 * 72 * 3 / 4);

//		whitenedNVS[0].showNeuronsInSubclasses = true;
//		whitenedNVS[1].showNeuronsInSubclasses = true;
//		whitenedNVS[2].showNeuronsInSubclasses = true;
        nvs[0].showNeuronsInSubclasses = true;
        nvs[1].showNeuronsInSubclasses = true;
        nvs[2].showNeuronsInSubclasses = true;
        JPanel panel = new JPanel(new PlaceLayout());
        int cpNumber = 1;
        int currentNumber = 0;
        for (int mosaicNumber = 0; mosaicNumber < dsClassList.length; mosaicNumber++) {
            if (mosaicNumber == 0) {
                cpNumber = 4;
            } else if (mosaicNumber == 3) {
                cpNumber = 13;

            } else if (mosaicNumber == 4) {
                cpNumber = 21;
            }
//			whitenedNVS[dsPrepList[mosaicNumber]].select(dsClassList[mosaicNumber]);
            nvs[dsPrepList[mosaicNumber]].select(dsClassList[mosaicNumber]);

//			PlotPanel p = makeMosaicPlot(whitenedNVS[dsPrepList[mosaicNumber]],
            PlotPanel p = makeMosaicPlot(nvs[dsPrepList[mosaicNumber]],
                    leftPadding,
                    rightPadding,
                    topPadding, bottomPadding,
                    verticalAxisThickness,
                    width, staHeight, labelFont, letterFont,
                    "CP" + cpNumber, Color.black);
            p.setAxisVisible(false);
            if (currentNumber == 0) {
                p.setScaleBar(.500, 1, 0.05, 0.1);
            }

            panel.add(p,
                    new PlaceC( (currentNumber % 4) / 4.0,
                            (1.0 / 6.0) * (currentNumber / 4)
                            , 0, 0, 1 / 4.0, 1 / 6.0));

            frame.add(panel);
            cpNumber++;
            currentNumber++;
        }
        frame.setVisible(true);

        RepaintManager.currentManager(frame).setDoubleBufferingEnabled(false);
        EpsGraphics2D g = new EpsGraphics2D(
                "", Font.decode("Arial ITALIC 8"), 0, 0, frame.getWidth(), frame.getHeight());

        g.scale(72.0 / dpi, 72.0 / dpi);
        g.setAccurateTextMode(true);
        frame.paint(g);

        RepaintManager.currentManager(frame).setDoubleBufferingEnabled(true);
        try {
            FileWriter writer = new FileWriter(
            "e:\\Publications\\GuineaPaper\\figures\\classesDS.eps");
            writer.write(g.toString());
            writer.flush();
            writer.close();

        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public static void makeNonDSMosaicsPlot() {
        JFrame frame = new JFrame();

        double dpi = 72.0;
        int width = (int) (3.2 * (dpi));
        int height = (int) (2.3 /*1.7*/ * (dpi));
        int staHeight = (int) (1.6 /*1.6*/ * (dpi));
        int verticalAxisThickness = 0;
        int rightPadding = 2;
        int topPadding = 1;
        int bottomPadding = 1;
        int leftPadding = 0;
        Font labelFont = Font.decode("Arial ITALIC 10");
        Font letterFont = Font.decode("Arial Bold 10");
        frame.setSize(7 * 72, 7 * 72 * 3 / 4);

//		whitenedNVS[0].showNeuronsInSubclasses = true;
//		whitenedNVS[1].showNeuronsInSubclasses = true;
//		whitenedNVS[2].showNeuronsInSubclasses = true;
        nvs[0].showNeuronsInSubclasses = true;
        nvs[1].showNeuronsInSubclasses = true;
        nvs[2].showNeuronsInSubclasses = true;
        JPanel panel = new JPanel(new PlaceLayout());
        int cpNumber = 1;
        int currentNumber = 0;
        for (int mosaicNumber = 0; mosaicNumber < classList.length; mosaicNumber++) {
            if (mosaicNumber == 3) {
                mosaicNumber++;
                cpNumber += 3;
            } else if (mosaicNumber == 10) {
                mosaicNumber++;
                cpNumber += 1;
            } else if (mosaicNumber == 18) {
                mosaicNumber++;
                cpNumber += 4;
            }
//			whitenedNVS[prepList[mosaicNumber]].select(classList[mosaicNumber]);
            nvs[prepList[mosaicNumber]].select(classList[mosaicNumber]);

//			PlotPanel p = makeMosaicPlot(whitenedNVS[prepList[mosaicNumber]],
            PlotPanel p = makeMosaicPlot(nvs[prepList[mosaicNumber]],
                    leftPadding,
                    rightPadding,
                    topPadding, bottomPadding,
                    verticalAxisThickness,
                    width, staHeight, labelFont, letterFont,
                    "CP" + cpNumber, Color.black);
            p.setAxisVisible(false);
            if (cpNumber == 1) {
                p.setScaleBar(.500, 1, 0.05, 0.1);
            }

            panel.add(p,
                    new PlaceC( (currentNumber % 4) / 4.0,
                            (1.0 / 6.0) * (currentNumber / 4)
                            , 0, 0, 1 / 4.0, 1 / 6.0));

            frame.add(panel);
            cpNumber++;
            currentNumber++;
        }
        frame.setVisible(true);

        RepaintManager.currentManager(frame).setDoubleBufferingEnabled(false);
        EpsGraphics2D g = new EpsGraphics2D(
                "", Font.decode("Arial ITALIC 8"), 0, 0, frame.getWidth(), frame.getHeight());

        g.scale(72.0 / dpi, 72.0 / dpi);
        g.setAccurateTextMode(true);
        frame.paint(g);

        RepaintManager.currentManager(frame).setDoubleBufferingEnabled(true);
        try {
            FileWriter writer = new FileWriter(
            "e:\\Publications\\GuineaPaper\\figures\\classesNonDS.eps");
            writer.write(g.toString());
            writer.flush();
            writer.close();

        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public static void makeTotalMosaicsPlot(int prepNumber) {
        JFrame frame = new JFrame();

        double dpi = 72.0;
        int width = (int) (3.2 * (dpi));
        int height = (int) (2.3 /*1.7*/ * (dpi));
        int staHeight = (int) (1.6 /*1.6*/ * (dpi));
        int verticalAxisThickness = 0;
        int rightPadding = 2;
        int topPadding = 1;
        int bottomPadding = 1;
        int leftPadding = 0;
        Font labelFont = Font.decode("Arial ITALIC 10");
        Font letterFont = Font.decode("Arial Bold 10");
        frame.setSize(7 * 72, 7 * 72 );

//		whitenedNVS[0].showNeuronsInSubclasses = true;
//		whitenedNVS[1].showNeuronsInSubclasses = true;
//		whitenedNVS[2].showNeuronsInSubclasses = true;
        nvs[0].showNeuronsInSubclasses = true;
        nvs[1].showNeuronsInSubclasses = true;
        nvs[2].showNeuronsInSubclasses = true;
        JPanel panel = new JPanel(new PlaceLayout());
        int cpNumber = 1;
        int currentNumber = 0;
        for (int mosaicNumber = 0; mosaicNumber < allClassList.length; mosaicNumber++) {
            prepNumber = allPrepList[mosaicNumber];
            cpNumber = mosaicNumber+1;
//			if (mosaicNumber == 3) {
//			mosaicNumber++;
//			cpNumber += 3;
//			} else if (mosaicNumber == 10) {
//			mosaicNumber++;
//			cpNumber += 1;
//			} else if (mosaicNumber == 18) {
//			mosaicNumber++;
//			cpNumber += 4;
//			}
//			whitenedNVS[prepList[mosaicNumber]].select(classList[mosaicNumber]);


            DefaultMutableTreeNode node = nvs[prepNumber].select(allClassList[mosaicNumber]);

            PlotPanel p;
            if(node!=null) {
//				PlotPanel p = makeMosaicPlot(whitenedNVS[prepList[mosaicNumber]],
                p = makeMosaicPlot(nvs[prepNumber],
                        leftPadding,
                        rightPadding,
                        topPadding, bottomPadding,
                        verticalAxisThickness,
                        width, staHeight, labelFont, letterFont,
                        "CP" + cpNumber, Color.black);
                p.setAxisVisible(false);
                if (cpNumber == 1) {
                    p.setScaleBar(.500, 1, 0.05, 0.1);
                }
            } else {
                System.out.println(prepNumber);
                System.out.println(mosaicNumber);
                System.out.println(allClassList[mosaicNumber]);
                System.out.println();
                p = new PlotPanel();
                p.setAxisVisible(false);
                p.axesBorder.setRightPadding(rightPadding);
                p.axesBorder.setTopPadding(topPadding);
                p.axesBorder.setBottomPadding(bottomPadding);
                p.axesBorder.setLeftPadding(leftPadding);
                p.axesBorder.getVerticalAxis().setThickness(verticalAxisThickness);
                p.axesBorder.getVerticalAxis().setLabelPosition(1);
            }

            panel.add(p,
                    new PlaceC( (currentNumber % 4) / 4.0,
                            (1.0 / 8.0) * (currentNumber / 4)
                            , 0, 0, 1 / 4.0, 1 / 8.0));

            frame.add(panel);
            cpNumber++;
            currentNumber++;
        }
        frame.setVisible(true);

        RepaintManager.currentManager(frame).setDoubleBufferingEnabled(false);
        EpsGraphics2D g = new EpsGraphics2D(
                "", Font.decode("Arial ITALIC 8"), 0, 0, frame.getWidth(), frame.getHeight());

        g.scale(72.0 / dpi, 72.0 / dpi);
        g.setAccurateTextMode(true);
        frame.paint(g);

        RepaintManager.currentManager(frame).setDoubleBufferingEnabled(true);
        try {
            FileWriter writer = new FileWriter(
                    "e:\\Publications\\GuineaPaper\\figures\\mosaics" + (prepNumber+10) + ".eps");
            writer.write(g.toString());
            writer.flush();
            writer.close();

        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public static void makeAllClassFigures() throws Exception {
        for(int i=0; i<allClassList.length; i++) {
    //	for(int j=0; j<3; j++) {
    
        makeClassFigure(2, allClassList[i], "CP" + new Integer(i+1).toString());
    //	makeClassFigure(0, allClassList[0], "CP" + new Integer(0+1).toString());

        }
        //}


    }

    public static void makeClassFigures() throws Exception {
        makeClassFigure(0, "All/On/Small/T2/OMS", "CP1", "67");
//		makeClassFigure(2, "All/On/Small/T3/OMS", "CP2", "264");
//		makeClassFigure(2, "All/On/Large/T1/Triphasic", "CP3", "5");
//		makeClassFigure(0, "All/On/Large/T2/DS/3", "CP4-6", "24");
//		makeClassFigure(0, "All/On/Large/T2/AOMS", "CP7", "12");
//		makeClassFigure(0, "All/On/Large/T2/Composite", "CP8", "253");
//		makeClassFigure(2, "All/On/Large/T2/A", "CP9", "188");
//		makeClassFigure(0, "All/On/Large/T3/UniformityA", "CP10", "175");
//		makeClassFigure(2, "All/On/Large/T3/UniformityB", "CP11", "408");
//		makeClassFigure(1, "All/On/Large/T3/Blue", "CP12", "508");
//		makeClassFigure(1, "All/On/Large/T3/OS", "CP13", "240");
//		makeClassFigure(2, "All/On/Large/T3/Triphasic", "CP14", "491");
//		makeClassFigure(0, "All/On/Large/T3/A", "CP15", "142");
//		makeClassFigure(2, "All/On/Large/T3/B", "CP16", "180");
//		makeClassFigure(0, "All/On/Large/T3/C", "CP17", "370");
//		makeClassFigure(0, "All/On/Large/T3/D", "CP18", "156");
//		makeClassFigure(2, "All/On/Large/T3/E", "CP19", "61");
//		makeClassFigure(0, "All/On/Large/T3/F", "CP20", "258");
//		makeClassFigure(0, "All/Off/Small/T2/DS/2", "CP21-24", "171");
//		makeClassFigure(0, "All/Off/Small/T2/OMS", "CP25", "28");
//		makeClassFigure(2, "All/Off/Small/T3/OMS", "CP26", "229");
//		makeClassFigure(2, "All/Off/Small/T3/A", "CP27", "99");
//		makeClassFigure(0, "All/Off/Small/T3/B", "CP28", "488");
//		makeClassFigure(1, "All/Off/Small/T3/C", "CP29", "489");
//		makeClassFigure(0, "All/Off/Large/T1/AOMS", "CP30", "16");
//		makeClassFigure(0, "All/Off/Large/T1/A", "CP31", "10");
//		makeClassFigure(2, "All/Off/Large/T3/A", "CP32", "463");

    }


    public static void calculateReferenceParameters() throws Exception {
        for (int i = 0; i < parameterFiles.length; i++) {

//			DefaultMutableTreeNode node = whitenedNVS[i].select(referenceType);
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
//			System.out.println("Gauss: " +calculateClassParameter(neuronsInClass, i,

//			"5.8*5*2*((SigmaX*SigmaY)^.5)").x);

            referenceContourArea[i] = calculateClassParameter(neuronsInClass, i,
                    "2*((contourArea/3.14159)^.5)"
            /* "5.8*5*2*((SigmaX*SigmaY)^.5)"*/).x;
            System.out.println("Contour: " + i + " " +
                    calculateClassParameter(neuronsInClass, i,
                    "2*((contourArea/3.14159)^.5)").toString(true));
            referenceRL[i] = calculateClassParameter(neuronsInClass, i, "rl").x;
            System.out.println("rl: " + i + " " +
                    calculateClassParameter(neuronsInClass, i, "rl").toString(true));

            referenceAcfMean[i] = calculateClassParameter(neuronsInClass, i, "acfMean").x;
            System.out.println("acfMean: " + i + " " +
                    calculateClassParameter(neuronsInClass, i, "acfMean").toString(true));

            referenceSpatialMean[i] = calculateClassParameter(neuronsInClass, i, "spatialMean").x;
            System.out.println("spatialMean: " + i + " " +
                    calculateClassParameter(neuronsInClass, i, "spatialMean").toString(true));

            referenceTemporalMean[i] = calculateClassParameter(neuronsInClass, i, "temporalMean").x;
            System.out.println("temporalMean: " + i + " " +
                    calculateClassParameter(neuronsInClass, i, "temporalMean").toString(true));

            referenceRLE[i] = calculateClassParameter(neuronsInClass, i, "rlE").x;
            System.out.println("rlE: " + i + " " +
                    calculateClassParameter(neuronsInClass, i, "rlE").toString(true));
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

//			HashMap values = whitenedParameterFiles[parameterFile].evaluate(expression,
//			neuronsInClass);

            HashMap values = parameterFiles[parameterFile].evaluate(expression,
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


    public static void calculateParametersOld(String[] classNames) throws Exception {

        int rightPadding = 28;
        int topPadding = 4;

        ScatterPlot scats[] = new ScatterPlot[6];
        Gaussian1DFunction funcs[][] = new Gaussian1DFunction[6][27];
        for (int i = 0; i < scats.length; i++) {
            scats[i] = new ScatterPlot();

        }

        String[][] paramNamesAndTypes = new String[8][2];
        paramNamesAndTypes[0]= new String[]{"ID", "Double"};
        paramNamesAndTypes[1]= new String[]{"classID", "String"};
        paramNamesAndTypes[2]= new String[]{"SCED", "Double"};
        paramNamesAndTypes[3]= new String[]{"AS", "Double"};
        paramNamesAndTypes[4]= new String[]{"RL", "Double"};
        paramNamesAndTypes[5]= new String[]{"SRM", "Double"};
        paramNamesAndTypes[6]= new String[]{"F2F1", "Double"};
        paramNamesAndTypes[7]= new String[]{"OMS", "Double"};
        double[] params = new double[6];


        ParametersFile paramsFile = new ParametersFile("e:\\data\\combined\\byClass\\byClass.params",
                paramNamesAndTypes, 10000);


        JPanel jpanel = new JPanel();
        jpanel.setLayout(new GridLayout(2, 3));
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
                                /*"2*((contourArea/3.14159)^.5)"*/
                        "5.8*5*2*((SigmaX*SigmaY)^.5)").div(
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
            params[0] = num.x;

            num = Num.averageRescaleErrors(axonList);
            if (classNames[i].equals(referenceType)) {
                num = new Num(num.x, 0.000);
            }
            scats[1].add(num.x, classNumber, 0, num.err());
            funcs[1][i] = new Gaussian1DFunction(1, num.x, num.err());
            params[1] = num.x;

            num = Num.averageRescaleErrors(rlList);
            if (classNames[i].equals(referenceType)) {
                num = new Num(num.x, 0.000);
            }
            scats[2].add(num.x, classNumber, 0, num.err());
            funcs[2][i] = new Gaussian1DFunction(1, num.x, num.err());
            params[2] = num.x;

            num = Num.averageRescaleErrors(srmList);
            scats[3].add(num.x, classNumber, 0, num.err());
            funcs[3][i] = new Gaussian1DFunction(1, num.x, num.err());
            params[3] = num.x;

            num = Num.averageRescaleErrors(f2f1List);
            scats[4].add(num.x, classNumber, 0, num.err());
            funcs[4][i] = new Gaussian1DFunction(1, num.x, num.err());
            params[4] = num.x;

            num = Num.averageRescaleErrors(omsList);
            scats[5].add(num.x, classNumber, 0, num.err());
            funcs[5][i] = new Gaussian1DFunction(1, num.x, num.err());
            params[5] = num.x;

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

            paramsFile.addRow(new Object[]{
                    new Double(classNumber), 
                    "All",
                    new Double(params[0]),
                    new Double(params[1]),
                    new Double(params[2]),
                    new Double(params[3]),
                    new Double(params[4]),
                    new Double(params[5])});

            paramNamesAndTypes[0]= new String[]{"ID", "Double"};
            paramNamesAndTypes[1]= new String[]{"classID", "String"};
            paramNamesAndTypes[2]= new String[]{"SCED", "Double"};
            paramNamesAndTypes[3]= new String[]{"AS", "Double"};
            paramNamesAndTypes[4]= new String[]{"RL", "Double"};
            paramNamesAndTypes[5]= new String[]{"SRM", "Double"};
            paramNamesAndTypes[6]= new String[]{"F2F1", "Double"};
            paramNamesAndTypes[7]= new String[]{"OMS", "Double"};

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

//		p.addBackgroundText("" + 'A', LEFT, TOP, letterFont, Color.black);
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

//		p.addBackgroundText("" + 'B', LEFT, TOP, letterFont, Color.black);

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

//		PlotPanel axonsOnly = new PlotPanel();
//		axonsOnly.addData(scats[1],
//		new ScatterPlotStyle("", SymbolType.NONE, 0, Color.black));
//		axonsOnly.setLabels("Axon Speed", "Type Number");
//		axonsOnly.setLabelFont(labelFont);
//		axonsOnly.axesBorder.setRightPadding(rightPadding);
//		axonsOnly.axesBorder.setTopPadding(topPadding);
//		axonsOnly.setDoubleBuffered(false);
//		axonsOnly.setRange(.3, 1.3, 0, 33);
//		axonsOnly.saveAsPNG(new File(
//		"e:\\Publications\\thesis\\figures\\results\\axonsScatters.png"), 3, 3);



        p = new PlotPanel();
        p.addData(scats[2],
                new ScatterPlotStyle("Scatter", SymbolType.FILLED_SQUARE, 3,
                        Color.black));
        f = new FunctionSum(funcs[2]);
        p.addData(f, new FunctionStyle("f", Color.red, 1));
//		p.addBackgroundText("" + 'C', LEFT, TOP, letterFont, Color.black);
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
//		p.addBackgroundText("" + 'D', LEFT, TOP, letterFont, Color.black);
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
//		p.addBackgroundText("" + 'D', LEFT, TOP, letterFont, Color.black);

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
//		p.addBackgroundText("" + 'D', LEFT, TOP, letterFont, Color.black);

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
        frame.setSize((int) (72 * 7.244), (int) (72 * 5.0));
        frame.add(jpanel);
        frame.setVisible(true);

        RepaintManager.currentManager(frame).setDoubleBufferingEnabled(false);

        EpsGraphics2D g = new EpsGraphics2D(
                "", Font.decode("Arial ITALIC 8"), 0, 0, frame.getWidth(), frame.getHeight());
        g.scale(1, 1);
        g.setAccurateTextMode(true);
        g.setColorDepth(g.RGB);
        frame.paint(g);

        FileWriter writer = new FileWriter(
        "e:\\Publications\\GuineaPaper\\figures\\paramsScatters.eps");
        writer.write(g.toString());
        writer.flush();
        writer.close();
//		frame.setVisible(false);



//		PlotUtil.showData("", p);
        table2.draw(System.out);

//		try {
//		table2.excelDraw(new PrintStream(new File(
//		"e:\\Publications\\thesis\\params.table")));
//		} catch (IOException e) {
//		e.printStackTrace();
//		}
        paramsFile.close(true);


    }

    public static void calculateParameters(String[] classNames) throws Exception {

        int records = 22;
        int rightPadding = 0;//28;
        int topPadding = 0;//4;

        ScatterPlotStyle scatterStyle = new ScatterPlotStyle(
                "Scatter", SymbolType.CROSS, 2, Color.black,
                false, Color.black,
                4f);
        scatterStyle.setErrorSymbolType(SymbolType.NONE);

        ScatterPlot scats[] = new ScatterPlot[records];
        for (int i = 0; i < scats.length; i++) {
            scats[i] = new ScatterPlot();

        }

        ScatterPlot scats2[] = new ScatterPlot[1];
        for(int i=0; i<scats2.length; i++) {
            scats2[i] = new ScatterPlot();
        }
        Gaussian1DFunction funcs[][] = new Gaussian1DFunction[1][27];

        double[][][] paramsAsMatrix = new double[27][records][2]; //class, variable, value or error
        int[] classNums = new int[27];




        JPanel jpanel = new JPanel();
        jpanel.setLayout(new GridLayout(records, 1));

        int currentRow = 0, currentColumn = 1;
        Table table2 = new Table(150, records +2);
        table2.setTitle(0, "");

        table2.setTitle(currentColumn++, "N");
        table2.setTitle(currentColumn++, "Diameter");
        table2.setTitle(currentColumn++, "Simplified Contour Diameter");
        table2.setTitle(currentColumn++, "Axon Speed");
        table2.setTitle(currentColumn++, "RL");
        table2.setTitle(currentColumn++, "SRM");
        table2.setTitle(currentColumn++, "F2/F1");
        table2.setTitle(currentColumn++, "OMS");
        table2.setTitle(currentColumn++, "blueness");
        table2.setTitle(currentColumn++, "magOS");
        table2.setTitle(currentColumn++,"magOnDS");
        table2.setTitle(currentColumn++,"magOffDS");
        table2.setTitle(currentColumn++, "acfMean");
        table2.setTitle(currentColumn++, "spatialMean");
        table2.setTitle(currentColumn++, "temporalMean");
        table2.setTitle(currentColumn++, "RLV");
        table2.setTitle(currentColumn++, "SRMV");
        table2.setTitle(currentColumn++, "amp1");
        table2.setTitle(currentColumn++, "amp2");
        table2.setTitle(currentColumn++, "amp3");
        table2.setTitle(currentColumn++, "amp1V");
        table2.setTitle(currentColumn++, "amp2V");
        table2.setTitle(currentColumn++, "amp3V");

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
            ArrayList<Num> onDSList = new ArrayList<Num> ();
            ArrayList<Num> offDSList = new ArrayList<Num> ();
            ArrayList<Num> acfMeanList = new ArrayList<Num>();
            ArrayList<Num> spatialMeanList = new ArrayList<Num>();
            ArrayList<Num> temporalMeanList = new ArrayList<Num>();
            ArrayList<Num> rlVList = new ArrayList<Num>();
            ArrayList<Num> srmVList = new ArrayList<Num>();
            ArrayList<Num> amp1List = new ArrayList<Num>();
            ArrayList<Num> amp2List = new ArrayList<Num>();
            ArrayList<Num> amp3List = new ArrayList<Num>();
            ArrayList<Num> amp1VList = new ArrayList<Num>();
            ArrayList<Num> amp2VList = new ArrayList<Num>();
            ArrayList<Num> amp3VList = new ArrayList<Num>();
            

            int nNeurons = 0;
            for (int j = 0; j < parameterFiles.length; j++) {
                table2.setCell(currentRow, 0, classNames[i] + "   " + (j + 1));
                currentColumn = 1;
                //	DefaultMutableTreeNode node = whitenedNVS[j].select(classNames[i]);
                DefaultMutableTreeNode node = nvs[j].select(classNames[i]);
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
                        /*"5.8*5*2*((SigmaX*SigmaY)^.5)"*/).div(
                                referenceContourArea[j]);
                        contourList.add(num);
                        table2.setCell(currentRow, currentColumn++, num.toString(true));

                        num = calculateClassParameter(neuronsInClass, j,
                                /*"2*((simpleContourArea/3.14159)^.5)"*/
                        "5.8*5*2*((SigmaX*SigmaY)^.5)").div(
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
                        table2.setCell(currentRow, currentColumn++, num.toString(true));

                        num = calculateClassParameter(neuronsInClass, j,
                        "magOnDS");
                        onDSList.add(num);
                        table2.setCell(currentRow, currentColumn++, num.toString(true));

                        num = calculateClassParameter(neuronsInClass, j,
                        "magOffDS");
                        offDSList.add(num);
                        table2.setCell(currentRow, currentColumn++, num.toString(true));

                        num = calculateClassParameter(neuronsInClass, j,
                        "acfMean").div(
                                referenceAcfMean[j]);
                        acfMeanList.add(num);

                        table2.setCell(currentRow, currentColumn++, num.toString(true));

                        num = calculateClassParameter(neuronsInClass, j,
                        "spatialMean").sub(
                                referenceSpatialMean[j]);
                        spatialMeanList.add(num);
                        table2.setCell(currentRow, currentColumn++, num.toString(true));

                        num = calculateClassParameter(neuronsInClass, j,
                        "temporalMean").sub(
                                referenceTemporalMean[j]);
                        temporalMeanList.add(num);
                        table2.setCell(currentRow, currentColumn++, num.toString(true));

                        num = calculateClassParameter(neuronsInClass, j,
                        "rlV").div(referenceRLE[j]);
                        rlVList.add(num);
                        table2.setCell(currentRow, currentColumn++, num.toString(true));
                        
                        num = calculateClassParameter(neuronsInClass, j,
                        "srmV");
                        srmVList.add(num);
                        table2.setCell(currentRow, currentColumn++, num.toString(true));
                        
                        num = calculateClassParameter(neuronsInClass, j,
                        "amp1");
                        amp1List.add(num);
                        table2.setCell(currentRow, currentColumn++, num.toString(true));
                        
                        num = calculateClassParameter(neuronsInClass, j,
                        "amp2");
                        amp2List.add(num);
                        table2.setCell(currentRow, currentColumn++, num.toString(true));
                        
                        num = calculateClassParameter(neuronsInClass, j,
                        "amp3");
                        amp3List.add(num);
                        table2.setCell(currentRow, currentColumn++, num.toString(true));
                        
                        num = calculateClassParameter(neuronsInClass, j,
                        "amp1V");
                        amp1VList.add(num);
                        table2.setCell(currentRow, currentColumn++, num.toString(true));
                        
                        num = calculateClassParameter(neuronsInClass, j,
                        "amp2V");
                        amp2VList.add(num);
                        table2.setCell(currentRow, currentColumn++, num.toString(true));
                        
                        num = calculateClassParameter(neuronsInClass, j,
                        "amp3V");
                        amp3VList.add(num);
                        table2.setCell(currentRow++, currentColumn++, num.toString(true));
                        
                        

                    } else {
                        System.err.println("No neurons for class " + classNames[i] +
                                " found in preparation " + j);
                    }
                } else System.err.println("Class " + classNames[i] +
                        " not found in preparation " + (j + 1));
            }


            int classNumber = i + 1;
            if (classNumber > 4) classNumber += 2;
            if (classNumber > 21) classNumber += 3;

            classNums[i] = classNumber;

            //on vs off
            if(classNumber <= 20) {
                scats[6].add(classNumber, 1, 0, 0);
                paramsAsMatrix[i][6][0] = 1;
                paramsAsMatrix[i][6][1] = 0;
            }
            else {
                scats[6].add(classNumber, -1, 0, 0);
                paramsAsMatrix[i][6][0] = -1;
                paramsAsMatrix[i][6][1] = 0;
            }

            Num num = Num.averageRescaleErrors(simplifiedContourList);
            if (classNames[i].equals(referenceType)) {
                num = new Num(num.x, 0.000);
            }
            scats[0].add(classNumber, num.x, num.err(), 0);
            paramsAsMatrix[i][0][0] = num.x;
            paramsAsMatrix[i][0][1] = num.err();


            num = Num.averageRescaleErrors(axonList);
            if (classNames[i].equals(referenceType)) {
                num = new Num(num.x, 0.000);
            }
            scats[1].add(classNumber, num.x, num.err(), 0);
            paramsAsMatrix[i][1][0] = num.x;
            paramsAsMatrix[i][1][1] = num.err();


            num = Num.averageRescaleErrors(rlList);
            if (classNames[i].equals(referenceType)) {
                num = new Num(num.x, 0.000);
            }
            scats[2].add(classNumber, num.x, num.err(), 0);
            paramsAsMatrix[i][2][0] = num.x;
            paramsAsMatrix[i][2][1] = num.err();


            num = Num.averageRescaleErrors(srmList);
            scats[3].add(classNumber, num.x, num.err(), 0);
            paramsAsMatrix[i][3][0] = num.x;
            paramsAsMatrix[i][3][1] = num.err();


            num = Num.averageRescaleErrors(f2f1List);
            scats[4].add(classNumber, num.x, num.err(), 0);
            paramsAsMatrix[i][4][0] = num.x;
            paramsAsMatrix[i][4][1] = num.err();


            num = Num.averageRescaleErrors(omsList);
            scats[5].add(classNumber, num.x, num.err(), 0);
            paramsAsMatrix[i][5][0] = num.x;
            paramsAsMatrix[i][5][1] = num.err();

            num = Num.averageRescaleErrors(blueList);
            scats[7].add(classNumber, num.x, num.err(), 0);
            paramsAsMatrix[i][7][0] = num.x;
            paramsAsMatrix[i][7][1] = num.err();

            num = Num.averageRescaleErrors(osList);
            scats[8].add(classNumber, num.x, num.err(), 0);
            paramsAsMatrix[i][8][0] = num.x;
            paramsAsMatrix[i][8][1] = num.err();

            num = Num.averageRescaleErrors(onDSList);
            scats[9].add(classNumber, num.x, num.err(), 0);
            paramsAsMatrix[i][9][0] = num.x;
            paramsAsMatrix[i][9][1] = num.err();

            num = Num.averageRescaleErrors(offDSList);
            scats[10].add(classNumber, num.x, num.err(), 0);
            paramsAsMatrix[i][10][0] = num.x;
            paramsAsMatrix[i][10][1] = num.err();

            num = Num.averageRescaleErrors(acfMeanList);
            if (classNames[i].equals(referenceType)) {
                num = new Num(num.x, 0.000);
            }
            scats[11].add(classNumber, num.x, num.err(), 0);
            paramsAsMatrix[i][11][0] = num.x;
            paramsAsMatrix[i][11][1] = num.err();
            funcs[0][i] = new Gaussian1DFunction(1, num.x, num.err());
            scats2[0].add(num.x, classNumber, 0, num.err());

            num = Num.averageRescaleErrors(spatialMeanList);
            if (classNames[i].equals(referenceType)) {
                num = new Num(num.x, 0.000);
            }
            scats[12].add(classNumber, num.x, num.err(), 0);
            paramsAsMatrix[i][12][0] = num.x;
            paramsAsMatrix[i][12][1] = num.err();

            num = Num.averageRescaleErrors(temporalMeanList);
            if (classNames[i].equals(referenceType)) {
                num = new Num(num.x, 0.000);
            }
            scats[13].add(classNumber, num.x, num.err(), 0);
            paramsAsMatrix[i][13][0] = num.x;
            paramsAsMatrix[i][13][1] = num.err();

            num = Num.averageRescaleErrors(rlVList);
            if (classNames[i].equals(referenceType)) {
                num = new Num(num.x, 0.000);
            }
            scats[14].add(classNumber, num.x, num.err(), 0);
            paramsAsMatrix[i][14][0] = num.x;
            paramsAsMatrix[i][14][1] = num.err();


            num = Num.averageRescaleErrors(srmVList);
            scats[15].add(classNumber, num.x, num.err(), 0);
            paramsAsMatrix[i][15][0] = num.x;
            paramsAsMatrix[i][15][1] = num.err();
            
            num = Num.averageRescaleErrors(amp1List);
            scats[16].add(classNumber, num.x, num.err(), 0);
            paramsAsMatrix[i][16][0] = num.x;
            paramsAsMatrix[i][16][1] = num.err();
            
            num = Num.averageRescaleErrors(amp2List);
            scats[17].add(classNumber, num.x, num.err(), 0);
            paramsAsMatrix[i][17][0] = num.x;
            paramsAsMatrix[i][17][1] = num.err();

            num = Num.averageRescaleErrors(amp3List);  //FIX.  DUMPING NANS OR INFS WHEN ERROR IS ZERO

            if(Double.isNaN(num.x) || Double.isInfinite(num.x)) num = new Num(0,0);
            if(Double.isNaN(num.err()) || Double.isInfinite(num.err())) num = new Num(num.x, 0.0);
            scats[18].add(classNumber, num.x, num.err(), 0);
            paramsAsMatrix[i][18][0] = num.x;
            paramsAsMatrix[i][18][1] = num.err();
            
            num = Num.averageRescaleErrors(amp1VList);
            scats[19].add(classNumber, num.x, num.err(), 0);
            paramsAsMatrix[i][19][0] = num.x;
            paramsAsMatrix[i][19][1] = num.err();
            
            num = Num.averageRescaleErrors(amp2VList);
            scats[20].add(classNumber, num.x, num.err(), 0);
            paramsAsMatrix[i][20][0] = num.x;
            paramsAsMatrix[i][20][1] = num.err();
            
            num = Num.averageRescaleErrors(amp3VList);
            scats[21].add(classNumber, num.x, num.err(), 0);
            paramsAsMatrix[i][21][0] = num.x;
            paramsAsMatrix[i][21][1] = num.err();

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
            table2.setCell(currentRow, currentColumn++,
                    Num.averageRescaleErrors(osList).toString(true));
            table2.setCell(currentRow, currentColumn++,
                    Num.averageRescaleErrors(onDSList).toString(true));
            table2.setCell(currentRow, currentColumn++,
                    Num.averageRescaleErrors(offDSList).toString(true));
            table2.setCell(currentRow, currentColumn++,
                    Num.averageRescaleErrors(acfMeanList).toString(true));
            table2.setCell(currentRow, currentColumn++,
                    Num.averageRescaleErrors(spatialMeanList).toString(true));
            table2.setCell(currentRow, currentColumn++,
                    Num.averageRescaleErrors(temporalMeanList).toString(true));
            table2.setCell(currentRow, currentColumn++,
                    Num.averageRescaleErrors(rlVList).toString(true));
            table2.setCell(currentRow, currentColumn++,
                    Num.averageRescaleErrors(srmVList).toString(true));
            table2.setCell(currentRow, currentColumn++,
                    Num.averageRescaleErrors(amp1List).toString(true));
            table2.setCell(currentRow, currentColumn++,
                    Num.averageRescaleErrors(amp2List).toString(true));
            table2.setCell(currentRow, currentColumn++,
                    Num.averageRescaleErrors(amp3List).toString(true));
            table2.setCell(currentRow, currentColumn++,
                    Num.averageRescaleErrors(amp1VList).toString(true));
            table2.setCell(currentRow, currentColumn++,
                    Num.averageRescaleErrors(amp2VList).toString(true));
            table2.setCell(currentRow++, currentColumn++,
                    Num.averageRescaleErrors(amp3VList).toString(true));

        }

        table2.draw(System.out);

        scats = similaritySort(paramsAsMatrix, classNums);

        for(int i=0; i<scats.length; i++) {
            PlotPanel p = new PlotPanel();
            p.addData(scats[i], scatterStyle);
            p.setAxisVisible(false);
            p.setLabelFont(labelFont);
            p.axesBorder.setRightPadding(rightPadding);
            p.axesBorder.setTopPadding(topPadding);
            p.setDoubleBuffered(false);
            p.setRange(-1, 27, -.5, 1.5);

            jpanel.add(p);
        }


        JFrame frame = new JFrame();
        frame.setSize((int) (72 * 7.244), (int) (72 * 10.0));
        frame.add(jpanel);
        frame.setVisible(true);

        RepaintManager.currentManager(frame).setDoubleBufferingEnabled(false);

        EpsGraphics2D g = new EpsGraphics2D(
                "", Font.decode("Arial ITALIC 8"), 0, 0, frame.getWidth(), frame.getHeight());
        g.scale(1, 1);
        g.setAccurateTextMode(true);
        g.setColorDepth(g.RGB);
        frame.paint(g);

        FileWriter writer = new FileWriter(
        "e:\\Publications\\GuineaPaper\\figures\\paramsScatters.eps");
        writer.write(g.toString());
        writer.flush();
        writer.close();
//		frame.setVisible(false);



//		PlotUtil.showData("", p);
    //	table2.draw(System.out);

//		try {
//		table2.excelDraw(new PrintStream(new File(
//		"e:\\Publications\\thesis\\params.table")));
//		} catch (IOException e) {
//		e.printStackTrace();
//		}


        PlotPanel p = new PlotPanel();
        p.addData(scats2[0],
                new ScatterPlotStyle("Scatter", SymbolType.FILLED_SQUARE, 3,
                        Color.black));
        FunctionSum f = new FunctionSum(funcs[0]);
        p.addData(f, new FunctionStyle("f", Color.red, 1));
        frame = new JFrame();
        frame.setSize((int) (72 * 7.244), (int) (72 * 5.0));
        frame.add(p);
        frame.setVisible(true);

    }

    //similarity sort.  Start from the one that is most unique (farthest from the mean).  
    //then find the one that is closest from the
    //remaining items.  Continue till done.

    public static ScatterPlot[] similaritySort(double[][][] m, int[] classNums) {
        int records = m.length;
        int fields = m[0].length;


        double min = Double.POSITIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;
        for(int field=0; field<fields; field++) {
            for(int record=0; record<records; record++) {
                if(m[record][field][0] < min) min = m[record][field][0];
                if(m[record][field][0] > max) max = m[record][field][0];
            }

            for(int record=0; record<records; record++) {
                m[record][field][0] = (m[record][field][0] - min)/(max-min);
                m[record][field][1] = (m[record][field][1])/(max-min);
            }

            min = Double.POSITIVE_INFINITY;
            max = Double.NEGATIVE_INFINITY;
        }


        double[][][] out = new double[records][fields][2];

        int []classNumsTemp = new int[classNums.length];

        //find average record
        double average[] = new double[fields];	
        for(int record=0; record<m.length; record++) {
            for(int field=0; field<m[0].length; field++) {
                average[field]+= m[record][field][0];
            }
        }
        for(int field=0; field<m[0].length; field++) {
            average[field]/=m.length;
        }

        //find least average record

        int furthestRecord = -1;
        double furthestDistance = 0.0;
        for(int record=0; record<records; record++) {
            double distance = 0.0;
            for(int field=0; field<fields; field++) {
                distance += (m[record][field][0] - average[field])
                *(m[record][field][0] - average[field]);
            }
            if(distance > furthestDistance) {
                furthestRecord = record;
                furthestDistance = distance;
            }
        }


        boolean[] sorted = new boolean[records];
        java.util.Arrays.fill(sorted, false);

//		java.lang.System.arraycopy(m[furthestRecord], 0, out[0], 0, m[0].length);

        for(int field=0; field<fields; field++) {
            out[0][field][0] = m[furthestRecord][field][0];
            out[0][field][1] = m[furthestRecord][field][1];
        }

        sorted[furthestRecord] = true;
        classNumsTemp[0] = classNums[furthestRecord];

        for(int recordSorted = 1; recordSorted<records; recordSorted++) {
            int closestRecord = -1;
            double closestDistance = Double.POSITIVE_INFINITY;
            for(int record=0; record<records; record++) {
                if(!sorted[record]) {
                    double distance = 0.0;
                    for(int field=0; field<fields; field++) {
                        distance += (m[record][field][0] - out[recordSorted-1][field][0])
                        *(m[record][field][0] - out[recordSorted-1][field][0]);
                    }
                    if(distance < closestDistance) {
                        closestRecord = record;
                        closestDistance = distance;
                    }
                }
            }
            for(int field=0; field<fields; field++) {
                out[recordSorted][field][0] = m[closestRecord][field][0];
                out[recordSorted][field][1] = m[closestRecord][field][1];
            }
            //java.lang.System.arraycopy(m[closestRecord], 0, out[recordSorted], 0, m[0].length);

            sorted[closestRecord] = true;
            classNumsTemp[recordSorted] = classNums[closestRecord];
        }
        ScatterPlotStyle scatterStyle = new ScatterPlotStyle(
                "Scatter", SymbolType.CROSS, 2, Color.black,
                false, Color.black,
                4f);
        scatterStyle.setErrorSymbolType(SymbolType.NONE);

        ScatterPlot scats[] = new ScatterPlot[fields];

        for (int i = 0; i < scats.length; i++) {
            scats[i] = new ScatterPlot();

        }
        for(int record=0; record<records; record++) {
            for(int field = 0; field<fields; field++) {

                scats[field].add(record, out[record][field][0], out[record][field][1],0);
            }
        }

        for(int record = 0; record<records; record++) {
            classNums[record] = classNumsTemp[record];
            System.out.println(classNums[record]);
        }

        return scats;
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
                "e:\\data\\2005-05-02-0\\data002\\data002.params"
        };
        String[] name = {
                "e:\\data\\2005-05-02-0\\data005"
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
                    "e:\\Publications\\thesis\\figures\\stimuli\\contrastReversingExample1.eps",
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


    public static void makeParametersFile(String configName) throws Exception {
        Config config = new Config(configName);
        String cName = "Make Parameters File";
        HashMap<String, String> p= config.getParameterList(cName);

        Vision.getInstance().getCalculationManager().runCalculation(cName, p);
    }
    


    public static void mainAlt(String[] args) {
        String path = "c:\\Documents and Settings\\mgrivich\\eclipse\\workspace\\vision6\\";
        String[] configNames = {"0502.xml", "0504.xml", "0825.xml"};
        try {
            for(int i=0; i<configNames.length; i++) {
                makeParametersFile(path + configNames[i]);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    public static void main(String[] args) throws Exception {
        UIDefaults defaults = UIManager.getDefaults();
        defaults.put("Panel.background", new ColorUIResource(255, 255, 255));
        SwingUtilities.invokeAndWait(new Runnable() {
            public void run() {
                try {
                    parameterFiles = new ParametersFile[paramFilesString.length];
                    nvs = new NeuronViewer[paramFilesString.length];
                    globalsFiles = new GlobalsFile[paramFilesString.length];
                    for (int i = 0; i < parameterFiles.length; i++) {
                        parameterFiles[i] = new ParametersFile(paramFilesString[i]);
                        nvs[i] = new NeuronViewer(paramFilesString[i], false, new Config("config.xml"));
                        globalsFiles[i] = new GlobalsFile(globalsFilesString[i], GlobalsFile.READ);
                    }

                    whitenedParameterFiles = new ParametersFile[whitenedParamFilesString.
                                                                length];
                    whitenedNVS = new NeuronViewer[whitenedParamFilesString.length];
                    for (int i = 0; i < parameterFiles.length; i++) {
                        whitenedParameterFiles[i] = new ParametersFile(
                                whitenedParamFilesString[i]);
                        whitenedNVS[i] = new NeuronViewer(whitenedParamFilesString[i], false, new Config("config.xml"));
                    }

//					makeDSMosaicsPlot();
//					makeNonDSMosaicsPlot();
//					makeTotalMosaicsPlot(4);
//					makeTotalMosaicsPlot(1);
//					makeTotalMosaicsPlot(2);
//					referenceSimplifiedContourArea = new double[parameterFiles.length];
//					referenceContourArea = new double[parameterFiles.length];
//					referenceRL = new double[parameterFiles.length];
//					referenceAxon = new double[parameterFiles.length];
//					referenceAcfMean = new double[parameterFiles.length];
//					referenceSpatialMean = new double[parameterFiles.length];
//					referenceTemporalMean = new double[parameterFiles.length];
//					referenceRLE = new double[parameterFiles.length];					
//					calculateReferenceParameters();
//					calculateParameters(classList);
//					makeClassFigures();
//					makeAllClassFigures();
//					makeMethodsFigures();
                    makeResultsFigures();

//					makeWhitenerFigure("/All/Off/Small/T3/A", "326");
//					makeWhitenerFigure("/All/On/Large/T3/C", "66");

//					makeDetailedF1F2();

                    /*                    makeSeparationQualityFigure(new String[] {
                              "e:\\data\\2005-05-02-0\\data002\\data002.params",
                     "e:\\data\\2005-05-02-0\\data002firsthalf\\data002firsthalf.params",
                              "e:\\data\\2005-05-02-0\\data001\\data001.params"}
                                                                    , 1);
                     */
                    /*                      makeSeparationQualityFigure(new String[] {
                                "e:\\data\\2005-05-04-0\\data002\\data002.params",
                     "e:\\data\\2005-05-04-0\\data002firsthalf\\data002firsthalf.params",
                                "e:\\data\\2005-05-04-0\\data001\\data001.params"}
                                                                      , 2);

                                          makeSeparationQualityFigure(new String[] {
                                "e:\\data\\2005-08-25-0\\data003\\data003.params",
                     "e:\\data\\2005-08-25-0\\data003firsthalf\\data003firsthalf.params",
                                "e:\\data\\2005-08-25-0\\data002\\data002.params"},
                                                                      3);
                     */
                } catch (Exception ex) {
                    ex.printStackTrace();
                }
            }
        });
            //	System.exit(1);
    }
}
