package edu.ucsc.neurobiology.vision.neuronviewer;

import java.io.*;
import java.util.*;

import java.awt.*;
import javax.swing.*;
import javax.swing.tree.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.parameters.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.stimulus.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * @author Alexander Sher, University of California, Santa Cruz
 */
public class MakePixelMaskAction
    extends CalculationAction {

    public MakePixelMaskAction() {
        super("Make Pixel Mask", CalculationAction.CLASS_ACTION);
    }


    public void doAction(final IntegerList list, final TreePath classTreePath) {
//        SwingWorker worker = new SwingWorker() {
        double pixelMaskCutOff;
        int frameUsed, colorUsed, nFramesForBaseline, nFramesForTemplate;
        boolean useTemplate, noAdjacencyRequirement;
        PlotMaker plotMaker;
//        double chiSquared = 5;

        ParametersTable table = viewer.configuration.showDialog(
            "Make Pixel Mask", "Make Pixel Mask", viewer.mainFrame);
        if (table == null) {
            return;
        }

        frameUsed = table.getIntParameter("Frame to use");
        colorUsed = table.getIntParameter("Color to use");
        pixelMaskCutOff = table.getDoubleParameter("Pixel Mask Cutoff");
        noAdjacencyRequirement = table.getBooleanParameter("No Adjacency Requirement");
        nFramesForBaseline = table.getIntParameter("Number of Frames For the Baseline");
        useTemplate = table.getBooleanParameter("Use Template");
        nFramesForTemplate = table.getIntParameter("Use Template.Number of Frames");
        String templateClass = table.getStringParameter("Use Template.Class");
        int templateColor = (int)
                            ( (EnumeratorParameter) table.getParameter(
                                "Use Template.Color")).
                            getValue();

//    int templateColor = (int)
//        ( (EnumeratorParameter) table.getParameter("Use Template.Color")).
//        getTextFor(( (EnumeratorParameter) table.getParameter("Use Template.Color")).
//        getValue());

// String kuku = ( (EnumeratorParameter) table.getParameter("Use Template.Color")).
//                getTextFor(templateColor);



//        public String getStringParameter(String name) {
//        return ( (StringParameter) getParameter(name)).getValue()


        ArrayList<Integer> neuronsInClass = new ArrayList();
//                    IntegerList duplicateNeurons = new IntegerList();
        for (int i = 0; i < list.size(); i++) {
            neuronsInClass.add(list.get(i));
        }

        DoubleHistogram h = new DoubleHistogram("", -10, 10, 0.1);
        DoubleHistogram h1 = new DoubleHistogram("", 0, 100, 1);
        DoubleHistogram h2 = new DoubleHistogram("", -0.5, 10.5, 1);
//        DoubleHistogram h2 = new DoubleHistogram("", -15, 15, 0.1);
        MeanVarianceCalculator mvcSconeNumber = new MeanVarianceCalculator(
            MeanVarianceCalculator.UNBIASED);
        MeanVarianceCalculator mvcPixelToCellCount = new MeanVarianceCalculator(
            MeanVarianceCalculator.UNBIASED);

        double[][] pixelMask = new double[3][staFile.getHeight() * staFile.getWidth()];
        int[][] pixelMask1 = new int[staFile.getWidth()][staFile.getHeight()];
        int[][] pixelMask2 = new int[staFile.getWidth()][staFile.getHeight()];
        double[] template = new double[staFile.getSTADepth()];

        int[][] pixelToCellCount = new int[staFile.getWidth()][staFile.getHeight()];

        Arrays.fill(pixelMask[0], -1);
        Arrays.fill(pixelMask[1], -1);
        Arrays.fill(pixelMask[2], 0);

        for (int i = 0; i < pixelMask.length; i++) {
            Arrays.fill(pixelMask1[i], 0);
            Arrays.fill(pixelMask2[i], 0);
            Arrays.fill(pixelToCellCount[i], 0);
        }

        int npixels = 0;
        int npixelsPerNeuron = 0;
        double intensity, error, baseline, intensityCentral;
        double response[];
//            double averageClusterSignal = 0;
        int width, height, length;
//            int nLess = 0;

        width = staFile.getWidth();
        height = staFile.getHeight();
        length = staFile.getSTADepth();

        if (useTemplate) {

            int templateL[] = paramsFile.getNeuronsInClass(templateClass);

            HashMap<Integer, double[]> v;
            try {
                v = paramsFile.evaluate( ( (EnumeratorParameter) table.getParameter(
                    "Use Template.Color")).getTextFor(templateColor) + "TimeCourse",
                                        templateL);
            } catch (CannotEvaluateException ex) {
                System.out.println("Cannot calculate the template");
                v = null;
            }

            // average

            Arrays.fill(template, 0);
            int n = 0;
            for (double[] x : v.values()) {
                if (x != null) {
                    MathUtil.add(template, x);
                    n++;
                }
            }
            if (n == 0) {
                System.out.println("No neurons selected");
            }
            MathUtil.divide(template, n);
        }

        for (int i = neuronsInClass.size() - 1; i >= 0; i--) {
            Vision.getInstance().setProgress(
                100 * (neuronsInClass.size() - i) / neuronsInClass.size());

            int currentNeuron = neuronsInClass.get(i);
            STA sta = null;
            try {
                sta = staFile.getSTA(currentNeuron);
            } catch (IOException ex) {
                ex.printStackTrace();
                continue;
            }

//            sta.getMainFrame();
            ImageFrame frame = sta.getMainFrame();
//            ImageFrame frame = sta.getFrame(frameUsed);
//                ImageFrame frameL = sta.getFrame(frameUsed - 1);
            Point p = new Point(0, 0);
            frame.getMaxAbsColor(p);
//            frame.getMaxAbsColor(p, colorUsed - 1);
            int x0 = p.x;
            int y0 = p.y;
//            System.out.println("id :" + currentNeuron + " X0 = " + x0 + "  Y0 = " + y0);
            int xmin = Math.max(x0 - 15, 0);
            int ymin = Math.max(y0 - 15, 0);
            int xmax = Math.min(x0 + 15, width - 1);
            int ymax = Math.min(y0 + 15, height - 1);
            for (int x = xmin; x <= xmax; x++) {
                for (int y = ymin; y <= ymax; y++) {

                    response = getPixelResponse(sta, x, y, frameUsed, colorUsed,
                                                nFramesForBaseline, true,
                                                useTemplate, nFramesForTemplate, template,
                                                length);
                    intensity = response[0];
                    error = response[1];

                    h.fill(intensity / error, 1);
                    intensityCentral = intensity;
                    if (intensity / error > pixelMaskCutOff) {
//                       if ( (intensity / error > pixelMaskCutOff && !searchForMinSignal)
//                           || (intensity / error < pixelMaskCutOff && searchForMinSignal)) {
                        pixelMask[0][npixels] = x;
                        pixelMask[1][npixels] = y;
                        pixelMask[2][npixels] = intensity / error;
                        npixels++;
                        npixelsPerNeuron++;
//                            nLess = 0;

                        int xmin1 = Math.max(x - 1, 0);
                        int ymin1 = Math.max(y - 1, 0);
                        int xmax1 = Math.min(x + 1, width - 1);
                        int ymax1 = Math.min(y + 1, height - 1);
//                            double xCenter = 0;
//                            double yCenter = 0;
//                            double errorTotal = 0;
//                            double meanTotal = 0;
//                                            pixelMask[2][npixels - 1] = 0;
//                            double intensityTotal = 0;
//                            double normalization = 0;
                        if (!noAdjacencyRequirement) {
                            for (int x1 = xmin1; x1 <= xmax1; x1++) {
                                for (int y1 = ymin1; y1 <= ymax1; y1++) {
                                    response = getPixelResponse(sta, x1, y1, frameUsed,
                                        colorUsed, nFramesForBaseline, true,
                                        useTemplate, nFramesForTemplate, template, length);
                                    intensity = response[0];
                                    error = response[1];

                                    //                                                    pixelMask[2][npixels-1] = pixelMask[2][npixels-1] + Math.pow(intensity, 2);
                                    //                                                    errorTotal = errorTotal + Math.pow(error * intensity, 2);
                                    //                                                    meanTotal = meanTotal + Math.pow(error, 2);
                                    //                                                    xCenter = xCenter + x1*Math.abs(intensity);
                                    //                                                    yCenter = yCenter + y1*Math.abs(intensity);
                                    //                                                    normalization = normalization + Math.abs(intensity);

                                    if (intensityCentral < intensity) {
                                        pixelMask[0][npixels - 1] = -1;
                                        pixelMask[1][npixels - 1] = -1;
                                    }
                                }
                            }
                        }
                        if (pixelMask[0][npixels - 1] == -1) {
                            pixelMask[2][npixels - 1] = 0;
                            npixels = npixels - 1;
                            npixelsPerNeuron = npixelsPerNeuron - 1;
                        }
                        /*
                                                     else {
                         intensityTotal = Math.sqrt(pixelMask[2][npixels-1]);
                         errorTotal = Math.sqrt(errorTotal) / intensityTotal;
                         h2.fill((intensityTotal - Math.sqrt(meanTotal)) / errorTotal, 1);
                         if ((intensityTotal - Math.sqrt(meanTotal)) / errorTotal > 3) {
                         pixelMask[0][npixels - 1] = xCenter / normalization;
                         pixelMask[1][npixels - 1] = yCenter / normalization;
                         pixelMask[2][npixels - 1] = intensityTotal / errorTotal;
                                                          }
                                                          else {
                         pixelMask[0][npixels - 1] = -1;
                         pixelMask[1][npixels - 1] = -1;
                         pixelMask[2][npixels - 1] = 0;
                         npixelsPerNeuron = npixelsPerNeuron - 1;
                                                          }
                                                     }
                         */
                    }
                }
            }

            h1.fill(npixelsPerNeuron, 1);
            mvcSconeNumber.add(npixelsPerNeuron);
            npixelsPerNeuron = 0;
        }

        for (int i = 0; i < pixelMask[0].length; i++) {
            if (pixelMask[0][i] != -1) {
                pixelMask1[ (int) pixelMask[0][i]][ (int) pixelMask[1][i]] = 1;
                pixelToCellCount[ (int) pixelMask[0][i]][ (int) pixelMask[1][i]]++;
            }
        }

        for (int i = 0; i < staFile.getWidth(); i++) {
            for (int j = 0; j < staFile.getHeight(); j++) {
                if (pixelToCellCount[i][j] > 0) {
                    mvcPixelToCellCount.add(pixelToCellCount[i][j]);
                    h2.fill(pixelToCellCount[i][j], 1);
                }
            }
        }

        System.out.println("Average number of cells connected to one S cone = " +
                           mvcPixelToCellCount.getMean() +
                           " +/- " + mvcPixelToCellCount.getMeanVariance());
        System.out.println("Average number S cone per cell = " + mvcSconeNumber.getMean() +
                           " +/- " + mvcSconeNumber.getMeanVariance());

        if (table.getBooleanParameter("Save Pixel Mask to File")) {
            try {
                PrintWriter w = new PrintWriter
                                (new FileWriter(table.getFileParameter(
                                    "Save Pixel Mask to File.Mask File Path") +
                                                "/mask.txt"));
                for (int i = 0; i < pixelMask[0].length; i++) {
                    if (pixelMask[0][i] == -1) {
                        break;
                    }
                    w.println(pixelMask[0][i] + "  " + pixelMask[1][i]);
                }
                w.close();
            } catch (IOException e) {
                Vision.reportException("Can not write the mask to file", e);
            }
        }

        PlotPanel p = new PlotPanel();
        p.setLabels("Response in pixels", "");
        p.addData(h, new HistogramStyle("Spacing"));
        p.autoscale();

        PlotPanel p1 = new PlotPanel();
        p1.setLabels("Number of S cones per cell", "");
        p1.addData(h1, new HistogramStyle("Spacing"));
        p1.autoscale();
        p1.setLegendVisible(true);
        p1.addToLegend("Mean = " + StringUtil.format(mvcSconeNumber.getMean(), 8) +
                       "+/-" +
                       StringUtil.format(mvcSconeNumber.getMeanVariance(), 8));

        PlotPanel p2 = new PlotPanel();
        p2.setLabels("Number of cells Connected to an S cone", "");
        p2.addData(h2, new HistogramStyle("Spacing"));
        p2.autoscale();

        int nColumns = 3;

        JPanel panel = new JPanel(new GridLayout(0, nColumns));
        panel.add(p);
        panel.add(p1);
        panel.add(p2);
        JFrame frame = new JFrame();
        frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        frame.getContentPane().add(panel);
        frame.setBounds(0, 0, 900, 300);
        frame.setVisible(true);

//                            PlotPanel p2 = new PlotPanel();
//
//                            p2.setLabels("Response in clusters", "");
//                            p2.addData(h2, new HistogramStyle("Spacing"));
//                            p2.autoscale();
//                            PlotUtil.showData("Response in blue clusters", p2);


        // add the masked pixels plot to the mosaic and sta plots
        ScatterPlot sp = new ScatterPlot();
        for (int ii = 0; ii < pixelMask[0].length; ii++) {
            if (pixelMask[0][ii] != -1) {
                sp.add( (pixelMask[0][ii] + 0.5) * staFile.getStixelWidth(),
                       (height - pixelMask[1][ii] - 0.5) * staFile.getStixelHeight());
            }
        }
        ScatterPlotStyle style = new ScatterPlotStyle("Masked Pixels",
            SymbolType.FILLED_SQUARE, 3, Color.blue, false, Color.yellow, 1);

        CompoundPlotMaker mosaic = ( (CompoundPlotMaker) viewer.getPlotMaker("Mosaic"));
        mosaic.removePlot(style.getDescription());
        mosaic.addAdditionalPlot(sp, style);

        CompoundPlotMaker sta = ( (CompoundPlotMaker) viewer.getPlotMaker("STA"));
        sta.removePlot(style.getDescription());
        sta.addAdditionalPlot(sp, style);

        viewer.storeData("pixelMask", pixelMask);
        viewer.storeData("pixelMask1", pixelMask1);
        viewer.storeData("pixelMask2", pixelMask2);
        viewer.storeData("template", template);
    }


    public static double[] getPixelResponse(STA sta, int x, int y, int frameUsed,
                                            int colorUsed,
                                            int nFramesForBaseline, boolean useTwoFrames,
                                            boolean useTemplate,
                                            int nFramesForTemplate, double[] template,
                                            int length) {
        double response[] = new double[2];
        double baseline;

        baseline = 0;
        for (int ifr = 0; ifr < nFramesForBaseline; ifr++) {
            ImageFrame fr = sta.getFrame(ifr);
            baseline = baseline + fr.getPixel(x, y, colorUsed - 1) / nFramesForBaseline;
        }

        ImageFrame frame = sta.getFrame(frameUsed);
        ImageFrame frameL = sta.getFrame(frameUsed - 1);

        if (!useTemplate) {
            if (!useTwoFrames) {
                response[0] = frame.getPixel(x, y, colorUsed - 1) - baseline;
                response[1] = frame.getPixelError(x, y, colorUsed - 1);
            } else {
                if (frame.getPixel(x, y, colorUsed - 1) >
                    frameL.getPixel(x, y, colorUsed - 1)) {
                    response[0] = frame.getPixel(x, y, colorUsed - 1) - baseline;
                    response[1] = frame.getPixelError(x, y, colorUsed - 1);
                } else {
                    response[0] = frameL.getPixel(x, y, colorUsed - 1) - baseline;
                    response[1] = frameL.getPixelError(x, y, colorUsed - 1);
                }
            }
        } else {
            response[0] = 0;
            response[1] = 0;
            for (int ifr = length - 2; ifr > length - 2 - nFramesForTemplate; ifr--) {
                response[0] = response[0] +
                              (sta.getFrame(ifr).getPixel(x, y, colorUsed - 1) - baseline)
                              * template[ifr];
                response[1] = response[1] + Math.pow(sta.getFrame(ifr).getPixelError(x, y,
                    colorUsed - 1) * template[ifr], 2);
            }
            response[1] = Math.sqrt(response[1]);
        }

        return response;
    }


}
