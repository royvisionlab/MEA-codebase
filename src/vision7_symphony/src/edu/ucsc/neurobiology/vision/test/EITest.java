package edu.ucsc.neurobiology.vision.test;

import java.io.*;
import static java.lang.Math.*;
import java.util.*;

import static java.awt.Color.*;
import javax.swing.*;

import edu.ucsc.neurobiology.vision.electrodemap.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.math.fitting.*;
import edu.ucsc.neurobiology.vision.neuronviewer.*;
import static edu.ucsc.neurobiology.vision.plot.SymbolType.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 *
 * This calculation is now implemented in the main program.  This is a duplicate, stale copy.
 *
 * It is kept for the guinea pig information at the bottom.
 *
 */
public class EITest {
    PhysiologicalImagingFile f;
    ElectrodeMap map;
    int nElectrodes;
    float[] sigma;
    static ScatterPlotStyle style = new ScatterPlotStyle(SQUARE, 2, black, false, black,
        1);
    float[] template;


    public EITest(String name) throws IOException {
        this.f = new PhysiologicalImagingFile(name);
        System.err.println(f.getArrayID());
        map = ElectrodeMapFactory.getElectrodeMap(f.arrayID);
        sigma = IOUtil.loadFloatArray(StringUtil.removeExtension(name) + ".noise");
        nElectrodes = map.getNumberOfElectrodes();
    }


    public void setAxonTemplate(int id, int electrode) throws IOException {
        final float[][][] image = f.getImage(id);
        float[][] average = image[0];

        // prepare the template
        template = new float[average[0].length];
        System.arraycopy(average[electrode], 0, template, 0, template.length);
    }


    public void pca(final int id) throws IOException {
        final float[][][] img = f.getImage(id);

//        UniformSpline spline = new UniformSpline(f.nSamples);

        int N = 0;
        // do the PCA
        CovarianceMatrix c = new CovarianceMatrix(f.nSamples);
        final HashMap<Integer, float[]> vMap = new HashMap<Integer, float[]>();
//        DoubleHistogram h = new DoubleHistogram("", 0, 2, 0.0001);
        for (int e = 1; e < f.nElectrodes; e++) {
            if (map.isDisconnected(e)) {
                continue;
            }

            if (MathUtil.maxAbs(img[0][e]) < 15) {
                continue;
            }

            // amplitude normalization
            MathUtil.divide(img[0][e], MathUtil.sumAbs(img[0][e]));
//            MathUtil.divide(img[0][e], MathUtil.sum(img[0][e]));
//            MathUtil.divide(img[0][e], Math.abs(MathUtil.min(img[0][e])));

            // simple allignment
            int index = MathUtil.minIndex(img[0][e]);
            MathUtil.rotateArray(img[0][e], - (f.nlPoints - index));
            /*
                        // spline allignment
                        try {
                            spline.reSpline(img[0][e]);
                            double x0 = FunctionMinimization.brentParabolic(
                                spline, index - 1, index + 1, 0.001);
                            double n = x0 - f.nlPoints;
                            System.err.println(n);
                            for (int i = 0; i < f.nSamples - 1; i++) {
                                try {
                                    img[0][e][i] = (float) spline.getValueAt(i + n);
                                }catch (ArrayIndexOutOfBoundsException ex) {
                                    img[0][e][i] = 0;
                                }
                            }
                        } catch (CannotEvaluateException ex) {
                            ex.printStackTrace();
                        }
             */
            c.addData(img[0][e]);
            vMap.put(e, img[0][e]);
            N++;
        }
        System.err.println("N " + N);
//        PlotUtil.showData(h);

        PCA pca = new PCA(c.getCovariance());
        try {
            pca.doPCA();
        } catch (TooManyIterationsException e) {
            System.out.println("Impossible to Analyze");
            return;
        }

        // show the eigenvectors
//        pca.printPercentageEigenvalues(nEigenvectors);
//        pca.showEigenVectors(nEigenvectors);

        int c1 = 0;
        int c2 = 1;
        final HashMap<Integer, Double> xMap = new HashMap<Integer, Double>();
        final HashMap<Integer, Double> yMap = new HashMap<Integer, Double>();
        ScatterPlot sp = new ScatterPlot();
        for (Integer e : vMap.keySet()) {
            double x = pca.project(img[0][e], c1);
            double y = pca.project(img[0][e], c2);
            xMap.put(e, x);
            yMap.put(e, y);
            sp.add(x, y);
        }

        ScatterPlotStyle style = new ScatterPlotStyle();
        style.setSymbolSize(2);
        PlotPanel p = PlotUtil.showData(sp, style);
        p.pad();

        p.addSelectionAction(new SelectionAction("Print Electrodes") {
            public void selectionPerformed(JComponent source, Selection selection) {
                SelectionRegion r = selection.getSelection();
                IntegerList idList = new IntegerList();
                PlotPanel p = new PlotPanel();
                ScatterPlotStyle s = new ScatterPlotStyle();
                s.setConnectingPoints(true);
                for (Integer id : xMap.keySet()) {
                    Double ox = xMap.get(id);
                    Double oy = yMap.get(id);
                    double x = (ox != null) ? ox.doubleValue() : Double.NaN;
                    double y = (oy != null) ? oy.doubleValue() : Double.NaN;
                    if (r.contains(x, y)) {
                        idList.add(id.intValue());
                        p.addData(new ScatterPlot(null, vMap.get(id), null, null), s);
                    }
                }
                p.autoscale();
                PlotUtil.showData("", p);

                try {
                    PhysiologicalImagePanel pan = new PhysiologicalImagePanel(f.getImage(
                        id), null, 2,
                        map, -1);
                    pan.setMarkElectrodes(idList.toArray());
                    pan.showInAWindow("");
                } catch (IOException ex) {
                    ex.printStackTrace();
                }
            }
        });

    }


    public double[] average(int id) throws IOException {
        float[][][] image = f.getImage(id);
        float[][] average = image[0];

        double[] amplitude = new double[nElectrodes];
        for (int electrode = 0; electrode < nElectrodes; electrode++) {
//            amplitude[electrode] = -MathUtil.min(average[electrode]);
            amplitude[electrode] = MathUtil.maxAbs(average[electrode]);
        }
//        double minAmplitude = MathUtil.min(amplitude);
        double maxAmplitude = MathUtil.max(amplitude);
        int electrode = MathUtil.maxIndex(amplitude);

        int n = 30;
        double dMin = 300;
        double dAlpha = 2 * PI / n;

        double[] r = new double[n];

        for (int i = 0; i < n; i++) {
            double alpha1 = i * dAlpha;
            double alpha2 = (i + 1) * dAlpha;
            int N = 0;
            double amp = 0;

            for (int e = 1; e < nElectrodes; e++) {
                double dy = map.getYPosition(e) - map.getYPosition(electrode);
                double dx = map.getXPosition(e) - map.getXPosition(electrode);

                double alpha = atan2(dy, dx);
                if (alpha < 0) {
                    alpha += 2 * PI;
                }

                double d = sqrt(dx * dx + dy * dy);
                if (d > dMin && alpha > alpha1 && alpha <= alpha2) {
                    amp += amplitude[e];
                    N++;
                }
            }
            amp = (N == 0) ? 0 : amp / N;
            amp /= maxAmplitude;

            r[i] = amp;
        }

        int i = MathUtil.maxIndex(r);
//        System.err.println(i);
        MathUtil.rotateArray(r, - (n / 2 - i));

        return r;
    }


    public double getSpeed(int id) throws IOException {
        ScatterPlot sp = new ScatterPlot();

        final float[][][] image = f.getImage(id);
        if (image == null) {
            System.err.println("no image for " + id);
            return -1;
        }
        float[][] average = image[0];

        // find the cell body electrode
        double minAmp = Double.MAX_VALUE;
        int minAmpElectrode = -1;
        for (int electrode = 1; electrode < nElectrodes; electrode++) {
            double a = MathUtil.min(average[electrode]);
            if (a < minAmp) {
                minAmp = a;
                minAmpElectrode = electrode;
            }
        }

        // the time at which the cell body fires
        int T0 = MathUtil.minIndex(average[minAmpElectrode]);

//        DoubleHistogram h = new DoubleHistogram("", 0, 10, 0.05);
        final HashMap<Integer, Double> xMap = new HashMap<Integer, Double>();
        final HashMap<Integer, Double> yMap = new HashMap<Integer, Double>();
        double[][] x = new double[1][nElectrodes];
        double[] y = new double[nElectrodes];
        double[] err = new double[nElectrodes];
        int n = 0;
        for (int electrode = 1; electrode < nElectrodes; electrode++) {
            if (!map.isDisconnected(electrode)) {
                int t = MathUtil.minIndex(average[electrode]);
                double chi2 = compare(template, average[electrode]);
                if (t > T0 && Math.abs(average[electrode][t]) > 5) {
//                    h.fill(chi2, 1);

                    if (chi2 < 0.4) {
                        double d = map.getDistance(electrode, minAmpElectrode) / 1000;
                        double time = (t - T0) / 20.0;
                        sp.add(time, d);

                        xMap.put(electrode, (t - T0) / 20.0);
                        yMap.put(electrode, d);

                        x[0][n] = time;
                        y[n] = d;
                        err[n] = 0.1;

                        n++;
                    }
                }
            }
        }

//        PlotUtil.showData(h);

        if (n < 5) {
            return -1;
        }

        Linear1DFunction line = new Linear1DFunction(1, 0);
        Fitter fitter = new Fitter();
        try {
            fitter.fit(line, x, y, err, n);
        } catch (Exception ex) {
            return -1;
        }
//        System.err.println(line);

        PlotPanel p = /*new PlotPanel();/*/ PlotUtil.showData("" + id, sp, style, 600,
            300, 100, 500);
        p.addSelectionAction(new SelectionAction("Print Electrodes") {
            public void selectionPerformed(JComponent source, Selection selection) {
                SelectionRegion r = selection.getSelection();
                IntegerList list = new IntegerList();
                for (Integer id : xMap.keySet()) {
                    Double ox = xMap.get(id);
                    Double oy = yMap.get(id);
                    double x = (ox != null) ? ox.doubleValue() : Double.NaN;
                    double y = (oy != null) ? oy.doubleValue() : Double.NaN;
                    if (r.contains(x, y)) {
                        System.err.println(id);
                        list.add(id);
                    }
                }

                PhysiologicalImagePanel p = new PhysiologicalImagePanel(image, null, 2,
                    map, -1);
                p.setMarkElectrodes(list.toArray());
                p.showInAWindow("");
            }
        });

        p.setLabels("Time (ms)", "Distance (mm)");

        p.setRange(0, 3, 0, 2.5);
//        p.autoscale();

        p.addData(line, new FunctionStyle(""));
        return line.getSlope();
    }


    public double compare(float[] template, float[] data) {
        int t0 = MathUtil.minIndex(template);
        int t = MathUtil.minIndex(data);
        double s1 = MathUtil.sumAbs(template);
        double s2 = MathUtil.sumAbs(data);

        // i goes over the data
        double sum = 0;
        int n = 0;
        for (int i = 0; i < data.length; i++) {
            int j = i + (t0 - t);
            if (j >= 0 && j < data.length) {
                sum += Math.pow(template[j] / s1 - data[i] / s2, 2);
                n++;
            }
        }
        return 1000 * sum / n;
    }


    public void testCompare(int id) throws IOException {
        final float[][][] image = f.getImage(id);
        float[][] average = image[0];

//        System.err.println(StringUtil.format(
//            compare(average[templateElectrode], average[414]), 2));

        IntegerList list = new IntegerList();
        for (int i = 0; i < nElectrodes; i++) {
            int t = MathUtil.minIndex(average[i]);
            if (Math.abs(average[i][t]) > 5 &&
                compare(template, average[i]) < 0.4) {

                list.add(i);
            }
        }

        PhysiologicalImagePanel p = new PhysiologicalImagePanel(image, null, 2,
            map, -1);
        p.setMarkElectrodes(list.toArray());
        p.showInAWindow("" + id);
    }


    /*
                PhysiologicalImagingFile f = new PhysiologicalImagingFile(
                    "f:\\data\\2005-04-26-0\\data009\\data009.ei");

                int id = 675;
                ElectrodeMap map = ElectrodeMapFactory.getElectrodeMap(f.getArrayID());
     PhysiologicalImagePanel p = new PhysiologicalImagePanel(f.getImage(id), null, 0,
                    map, -1);
                p.setTemplates(153, 166);
                JFrame fr = new JFrame();
                fr.add(p);
                fr.setBounds(100, 100, 800, 400);
                fr.setVisible(true);
     */

    public static void main(String args[]) throws IOException {
        EITest t = new EITest("f:\\Good Data\\2005-04-26-0\\data009\\data009.ei");
        ParametersFile pf = new ParametersFile(
            "f:\\Good Data\\2005-04-26-0\\data009\\data009.params");
        t.setAxonTemplate(773, 165);

        int[] id = pf.getNeuronsInClass("All/OFF/Amacrine");
        DoubleHistogram h = new DoubleHistogram("", 0, 5, 0.05);
        MeanVarianceCalculator mvc = new MeanVarianceCalculator();
        for (int i = 0; i < id.length; i++) {
//                    t.testCompare(id[i]);
            double speed = t.getSpeed(id[i]);
            if (speed > 0) {
                mvc.add(speed);
                h.fill(speed, 1);
            }
        }
        PlotUtil.showData("", h);
        System.err.println(mvc);
    }

    /*** Calculation infrastructure for guinea pig data.  Please do not
       remove without consulting Matthew Grivich

     EITest t = new EITest("g:\\data\\2005-05-04-0\\data002\\data002.ei");
       ParametersFile pf = new ParametersFile(
           "g:\\data\\2005-05-04-0\\data002\\data002.params");
//    t.setAxonTemplate(312, 70);
       t.setAxonTemplate(166, 370);
//
//
//        EITest t = new EITest("g:\\data\\2005-05-02-0\\data002\\data002.ei");
//        ParametersFile pf = new ParametersFile(
//            "g:\\data\\2005-05-02-0\\data002\\data002.params");
//        t.setAxonTemplate(26, 22);
//                t.setAxonTemplate(243, 192); //blue
//t.setAxonTemplate(189,284);//sensitive

//      EITest t = new EITest("g:\\data\\2005-08-25-0\\data003\\data003.ei");
//    ParametersFile pf = new ParametersFile(
//        "g:\\data\\2005-08-25-0\\data003\\data003.params");
//      t.setAxonTemplate(155, 303);
//      t.setAxonTemplate(408,221);

           DoubleHistogram hist = new DoubleHistogram("", .1, 1.9, .0005);
//      for (int j = 0; j < GrivichThesis.classList.length; j++) {

           int[] id = pf.getNeuronsInClass("All/Off/X/T2/A");
           int j = 20;
//          int[] id = pf.getNeuronsInClass(GrivichThesis.classList[j]);

           DoubleHistogram h = new DoubleHistogram("", 0, 5, 0.05);
           MeanVarianceCalculator mvc = new MeanVarianceCalculator();
           for (int i = 0; i < id.length; i++) {
               t.testCompare(id[i]);
               double speed = t.getSpeed(id[i]);
               if (speed > 0 &&
     id[i] != -1  && id[i] !=185 && id[i] != 499 && id[i] != 352 && id[i] != 394) {
                   mvc.add(speed);
                   h.fill(speed, 1);
               }
           }
     PlotUtil.showData(GrivichThesis.classList[j], h, new HistogramStyle(), 400, 400,
                             700, 500);
           System.err.println(GrivichThesis.classList[j] + "   " + mvc.getMean() + " ± " + mvc.getMeanVariance());
           if (mvc.getMeanVariance() < 10 && mvc.getMeanVariance() != 0) {
               hist.fill(mvc.getMean(), 1.0);
           }

//      }
//        PlotUtil.showData("", hist);

//    }


     }
     */

}
