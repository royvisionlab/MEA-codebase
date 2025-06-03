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
public class MakeMaskedPixelColorHistogramAction
    extends CalculationAction {

    public MakeMaskedPixelColorHistogramAction() {
        super("Make Masked Pixel Color Histogram", CalculationAction.CLASS_ACTION);
    }


    public void doAction(final IntegerList list, final TreePath classTreePath) {
        int frameUsed, colorUsed, colorLookedAt, minNSpixels, nFramesForBaseline,
            nFramesForTemplate;
        double significanceCutOff;
        boolean negativeMask, useTemplate;
        PlotMaker plotMaker;
        String name;

        ParametersTable table = viewer.configuration.showDialog(
            "Make Masked Pixel Color Histogram", "Make Masked Pixel Color Histogram",
            viewer.mainFrame);
        if (table == null) {
            return;
        }

        frameUsed = table.getIntParameter("Frame to use");
        colorUsed = table.getIntParameter("Color to use");
        colorLookedAt = table.getIntParameter("Color to look at");
        minNSpixels = table.getIntParameter("Minimum number of masked pixels in RF");
        negativeMask = table.getBooleanParameter("Negative mask");
        significanceCutOff = table.getDoubleParameter("Significance cutoff");
        nFramesForBaseline = table.getIntParameter("Number of Frames For the Baseline");
        useTemplate = table.getBooleanParameter("Use Template");
        nFramesForTemplate = table.getIntParameter("Use Template.Number of Frames");
        String templateClass = table.getStringParameter("Use Template.Class");
        int templateColor = (int)
                            ( (EnumeratorParameter) table.getParameter(
                                "Use Template.Color")).
                            getValue();

        ArrayList<Integer> neuronsInClass = new ArrayList();
        for (int i = 0; i < list.size(); i++) {
            neuronsInClass.add(list.get(i));
        }

//        ScatterPlot sp = new ScatterPlot();
//        DoubleHistogram h[] = new DoubleHistogram[2];//("", -1, 1, 0.005);
        MeanVarianceCalculator mvc = new MeanVarianceCalculator(
            MeanVarianceCalculator.UNBIASED);
        double intensity, error, intensity1, significance;
        double error1, baseline, baseline1, blueint;
        double response[];
        int width, height, length;
        STA sta;
        width = staFile.getWidth();
        height = staFile.getHeight();
        length = staFile.getSTADepth();

        double[][] pixelMask = (double[][]) viewer.getStoredData("pixelMask");
        int[][] pixelMask1 = (int[][]) viewer.getStoredData("pixelMask1");
        int[][] pixelMask2 = (int[][]) viewer.getStoredData("pixelMask2");
        double[] template = (double[]) viewer.getStoredData("template");

        if (table.getBooleanParameter("Load Pixel Mask from File")) {
            Arrays.fill(pixelMask[0], -1);
            Arrays.fill(pixelMask[1], -1);
            Arrays.fill(pixelMask[2], 0);

            for (int i = 0; i < pixelMask.length; i++) {
                Arrays.fill(pixelMask1[i], 0);
                Arrays.fill(pixelMask2[i], 0);
            }

            try {
                FileReader fr = new FileReader(table.getFileParameter(
                    "Load Pixel Mask from File.Mask File Path") + "/mask.txt");
                LineNumberReader r = new LineNumberReader(fr);
                String line;
                int ipix = 0;
                while ( (line = r.readLine()) != null) {
                    line = line.trim();
                    int index = line.indexOf(' ');
                    pixelMask[0][ipix] = Double.parseDouble(line.substring(0, index));
                    pixelMask[1][ipix] = Double.parseDouble(line.substring(index + 2,
                        line.length()));
                    ipix++;
                }
                for (int i = 0; i < pixelMask[0].length; i++) {
                    if (pixelMask[0][i] != -1) {
                        pixelMask1[ (int) pixelMask[0][i]][ (int) pixelMask[1][i]] = 1;
                    }
                }
                fr.close();
            } catch (IOException e) {
                Vision.reportException("Can not read the mask file", e);
            }
        }

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

        try {
            int npixels = 0;
            for (int i = 0; i < pixelMask2.length; i++) {
                Arrays.fill(pixelMask2[i], 0);
            }
            int nColumns = 1;
            JPanel panel = new JPanel(new GridLayout(0, nColumns));

            String classes[] = {"All/OFF/Parasol", "All/ON/Parasol"};

//            negativeMask = true;
//            for (int kuku = 0; kuku < 2; kuku++) {
//                negativeMask = !negativeMask;
//                for (int nClass = 0; nClass < classes.length; nClass++) {
//
//                    int neuronList[] = paramsFile.getNeuronsInClass(classes[nClass]);
            DoubleHistogram h = new DoubleHistogram("", -1, 1, 0.005);
//                    for (int i = neuronList.length - 1; i >= 0; i--) {
            for (int i = neuronsInClass.size() - 1; i >= 0; i--) {

                int currentNeuron = neuronsInClass.get(i);
//                        int currentNeuron = neuronList[i];
                sta = staFile.getSTA(currentNeuron);
                npixels = 0;

                ImageFrame frame = sta.getMainFrame();
//            ImageFrame frame = sta.getFrame(frameUsed);
//                ImageFrame frameL = sta.getFrame(frameUsed - 1);
                Point p = new Point(0, 0);
                frame.getMaxAbsColor(p);
                int x0 = p.x;
                int y0 = p.y;
                int xmin = Math.max(x0 - 10, 0);
                int ymin = Math.max(y0 - 10, 0);
                int xmax = Math.min(x0 + 10, width - 1);
                int ymax = Math.min(y0 + 10, height - 1);
                if (minNSpixels > 0) {
                    for (int x = xmin; x <= xmax; x++) {
                        for (int y = ymin; y <= ymax; y++) {

                            response = MakePixelMaskAction.getPixelResponse(sta,
                                x, y,
                                frameUsed, colorUsed, nFramesForBaseline, false,
                                useTemplate, nFramesForTemplate, template, length);
                            intensity = response[0];
                            error = response[1];
                            significance = intensity / error;
                            if (Math.abs(significance) > significanceCutOff &&
                                pixelMask1[x][y] == 1) {
                                npixels++;
                            }
                        }
                    }
                }
                if (npixels >= minNSpixels) {
                    for (int x = xmin; x <= xmax; x++) {
                        for (int y = ymin; y <= ymax; y++) {

                            response = MakePixelMaskAction.getPixelResponse(sta, x, y,
                                frameUsed, colorUsed, nFramesForBaseline, false,
                                useTemplate, nFramesForTemplate, template, length);
                            intensity = response[0];
                            error = response[1];

                            significance = intensity / error;
                            if (Math.abs(significance) > significanceCutOff) {
                                if ( (!negativeMask && pixelMask1[x][y] == 1)
                                    || (negativeMask && pixelMask1[x][y] == 0)) {

                                    if (pixelMask1[x][y] == 0) {
                                        pixelMask2[x][y] = 1;
                                    }

                                    response = MakePixelMaskAction.
                                               getPixelResponse(sta, x, y, frameUsed,
                                        colorUsed, nFramesForBaseline, false,
                                        useTemplate, nFramesForTemplate, template,
                                        length);
                                    intensity = response[0];
                                    error = response[1];

                                    response = MakePixelMaskAction.
                                               getPixelResponse(sta, x, y, frameUsed,
                                        colorLookedAt, nFramesForBaseline, false,
                                        useTemplate, nFramesForTemplate, template,
                                        length);
                                    intensity1 = response[0];
                                    error1 = response[1];

                                    h.fill(intensity1 / intensity, 1);
                                    mvc.add(intensity1 / intensity);

                                    for (int ii = 0; ii < pixelMask[0].length; ii++) {
                                        if (pixelMask[0][ii] == x &&
                                            pixelMask[1][ii] == y) {
                                            blueint = pixelMask[2][ii];
                                        }
                                    }

                                    //                                           sp.add(currentNeuron, (frame.getPixel(x, y,
                                    //                                               colorLookedAt - 1) - baseline1) / (frame.getPixel(x, y,colorUsed - 1)
                                    //                                               -baseline)
                                    //);

                                    //                                           mvc.add((frame.getPixel(x, y,
                                    //                                               colorLookedAt - 1) - baseline1) / (frame.getPixel(x, y,colorUsed - 1)
                                    //                                               -baseline));
                                    //                                           System.out.println((frame.getPixel(x, y,
                                    //                                               colorLookedAt - 1) - baseline1) / (frame.getPixel(x, y,colorUsed - 1)
                                    //                                               -baseline));

                                    //                                               mvc.add((frame.getPixel(x, y,
                                    //                                                  colorLookedAt - 1) - baseline1));
                                    //                                              System.out.println(" X = "+x+" Y = "+y+" Ratio = "+
                                    //                                                   frame.getPixel(x, y,
                                    //                                                   colorLookedAt - 1) / frame.getPixel(x, y,colorUsed - 1)+" Ratio1 = "+
                                    //                                                   (frame.getPixel(x, y,
                                    //                                                   colorLookedAt - 1) - baseline1) / (frame.getPixel(x, y,colorUsed - 1)
                                    //                                                   -baseline));
                                }
                            }
                        }
                    }
                }
            }
            PlotPanel p = new PlotPanel();

            p.setLabels("Response", "");
            p.addData(h, new HistogramStyle("Spacing"));
            //                        p.addData(sp, new ScatterPlotStyle("Spacing"));
            p.autoscale();
            //                    p.setXRange(-0.04, 0.04);
            p.setLegendVisible(true);
            p.addToLegend("Mean = " + StringUtil.format(mvc.getMean(), 8) + "+/-" +
                          StringUtil.format(mvc.getMeanVariance(), 8));
            p.addToLegend(classTreePath.toString());
//                           PlotUtil.showData("kuku", p);

            panel.add(p);

            mvc.reset();
//                    h.clear();
//                }
//            }
            JFrame frame = new JFrame();
            frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
            frame.getContentPane().add(panel);
            frame.setBounds(0, 300, 900, 600);
            frame.setVisible(true);

        } catch (IOException e) {
            Vision.reportException(e);
        }

///////////////////////////////////////////




//                Vision.getInstance().sendMessage("Removing bad neurons...");
//                try {
//                    return makePixelMask(list, chiSquared);
//                } catch (IOException e) {
//                    Vision.reportException(e);
//                    return null;
//                } finally {
//                    Vision.getInstance().sendMessage("Done.");
//                    Vision.getInstance().setProgress(0);
//                }
//            }


//            public void finished() {
//                int[] neurons = (int[]) getValue();
//                if (neurons != null) {
//                    viewer.classificationHappened(neurons,
//                                                  classTreePath.pathByAddingChild(
//                        new DefaultMutableTreeNode("weak " + chiSquared, true)));
//                }
//            }
//        };
//        worker.start();
    }

}
