package edu.ucsc.neurobiology.vision.simulation;

import java.io.*;
import static java.lang.Math.*;
import java.util.*;

import java.awt.*;

import edu.ucsc.neurobiology.vision.io.*;
import static edu.ucsc.neurobiology.vision.math.MathUtil.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.math.fitting.*;
import edu.ucsc.neurobiology.vision.plot.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class SpatialSummation {
    DOG1DFunction subunitDOG;
    final int seed = 11111;
    final static double micr = 1 / 1000.0;
    double jitter = 0.0; // percent of spacing
    static double B = 0;

    final double noise;
    final double envelopeSize, s1, s2, vRatio, d;
    String name;
    double[] f1Data, f2Data;
    double[] freqs;


    public SpatialSummation(String name, double s, double s1, double s2,
                            double surrFraction, /* double farSurrFraction,*/
                            double noise, double[] f1Data, double[] f2Data,
                            double[] freqs) {
        this.envelopeSize = s;
        this.s1 = s1;
        this.s2 = s2;
        this.vRatio = surrFraction;
        this.d = 1 * s1;
        this.name = name;
        this.f1Data = f1Data;
        this.f2Data = f2Data;
        this.freqs = freqs;
        this.noise = noise;

        subunitDOG = new DOG1DFunction(1, surrFraction, 0, 1 / s1, 1 / s2);

        try {
            MathUtil.multiply(f1Data, f1.getValueAt(0) / f1Data[f1Data.length - 1]);
            MathUtil.multiply(f2Data, f2.getValueAt(0) / f2Data[f2Data.length - 1]);
        } catch (CannotEvaluateException ex) {
        }
    }


    public int getRange(double size) {
        return (int) Math.round(Math.sqrt(2) * size / d);
//        return 15;
    }


    public double getXi(Random random, int i) {
        return i * d + jitter * (random.nextDouble() - 0.5) * d;
    }


    public static double getWeight(double size, double xi) {
        return B + exp( -0.5 * pow(abs(xi) / size, 2));
    }


    public double getA() {
        Random random = new Random(seed);

        int range = getRange(envelopeSize);
        double v = 0;
        for (int i = 0; i <= range; i++) {
            double xi = getXi(random, i);
            if (i != 0) {
                v += 2 * getWeight(envelopeSize, xi);
            } else {
                v += 1 * getWeight(envelopeSize, xi);
            }
        }
        return v;
    }


    FunctionData f1 = new AbstractFunction() {
        Random random = new Random(seed);

        public double getValueAt(double w) {
            random.setSeed(seed);

            int range = getRange(envelopeSize);
            double v = 0;
            for (int i = -range; i <= range; i++) {
                double xi = getXi(random, i);

                double ratio = vRatio * (1 + 2 * (random.nextDouble() - 0.5) * jitter);
                double sigma1 = s1 * (1 + 2 * (random.nextDouble() - 0.5) * jitter);
                double sigma2 = s2 * (1 + 2 * (random.nextDouble() - 0.5) * jitter);
                double sub =
                    exp( -0.5 * w * w * sigma1 * sigma1) +
                    ratio * exp( -0.5 * w * w * sigma2 * sigma2);

                v += sub * getWeight(envelopeSize, xi) * cos(w * xi);
            }

            return abs(v);
        }
    };


    FunctionData f1Average = new AbstractFunction() {
        public double getValueAt(double w) {
            double p = .30;
            int n = 20;

            double sum = 0;
            for (int i = -n; i <= n; i++) {
                double wi = w * (1 + i * p / n);
                try {
                    sum += f1.getValueAt(wi);
                } catch (CannotEvaluateException ex) {
                }
            }

            return sum / (2 * n + 1);
        }
    };


    FunctionData f2WeightingFunction = new AbstractFunction() {
        Random random = new Random(seed);

        public double getValueAt(double w) {
            random.setSeed(seed);

            int range = getRange(envelopeSize);
            double v = 0;
            for (int i = -range; i <= range; i++) {
                double xi = getXi(random, i);

                // average over all phases
//                v += getWeight(envelopeSize, xi);

                // pi/4 phase
                v += getWeight(envelopeSize, xi) * abs(sin(w * xi + PI / 4));
            }

            return v;
        }
    };


    FunctionData f2 = new AbstractFunction() {
        Random random = new Random(seed);

        public double getValueAt(double w) {
            random.setSeed(seed);

            int range = getRange(envelopeSize);
            double v = 0;
            for (int i = -range; i <= range; i++) {
                double xi = getXi(random, i);

                // average over all phases
                v += getWeight(envelopeSize, xi);

                // pi/4 phase
//                v += getWeight(envelopeSize, xi) * abs(sin(w * xi + PI/4));
            }

            return abs(v) * subunitDOG.getValueAt(w) / 8 + noise;
        }
    };


    FunctionData sf = new FunctionData() {
        Random random = new Random(seed);

        public double getValueAt(double s) {
            random.setSeed(seed);

            int range = getRange(s);
            double sum = 0, sum1 = 0, sum2 = 0;
            for (int i = -range; i <= range; i++) {
                double xi = getXi(random, i);

                sum += getWeight(s, xi);
                sum1 += getWeight(s, xi) * getWeight(s1, xi);
                sum2 += getWeight(s, xi) * getWeight(s2, xi);
            }

            return (s2 + vRatio * s1) * sum / (s2 * sum1 + s1 * vRatio * sum2);
//            return (1 / s1 + vRatio * 1 / s2) /
//                (1 / sqrt(s1 * s1 + s * s) + vRatio / sqrt(s2 * s2 + s * s));
        }


        public String getDescription() {
            return "rf";
        }
    };


    public void showPlots(Rectangle r) throws Exception {
        System.err.println(s1);

        PlotPanel p = new PlotPanel();
        p.addData(f1, new FunctionStyle("F1", Color.black, 1));
//        p.addData(f1Average, new FunctionStyle("F1 Average", Color.black, 2));
        p.addData(new ScatterPlot(freqs, f1Data, null),
                  new ScatterPlotStyle(SymbolType.SQUARE, 4, Color.black, true,
                                       Color.black, 2));

        p.addData(f2, new FunctionStyle("F2", Color.red, 1));
        p.addData(new ScatterPlot(freqs, f2Data, null),
                  new ScatterPlotStyle(SymbolType.SQUARE, 4, Color.red, true,
                                       Color.red, 2));

//        p.addData(f2WeightingFunction, new FunctionStyle("F2 weighting", Color.cyan, 1));

//        p.addData(new DOG1DFourierFunction(
//            getA(), getA() * vRatio,
//            sqrt(s * s + s1 * s1),
//            sqrt(s * s + s2 * s2), 0
//                  ), new FunctionStyle("RF", Color.gray, 10));

//        p.addData(new DOG1DFunction(
//            getA(), getA() * vRatio, 0,
//            1 / sqrt(s * s + s1 * s1),
//            1 / sqrt(s * s + s2 * s2)
//                  ), new FunctionStyle("RF", Color.green, 1));

//        p.addData(new DOG1DFunction(
//            1.8, 0, 0,
//            1/0.100, 0
//                  ), new FunctionStyle("RF1", Color.magenta, 1));


        p.setRange(0.22 * 2 * PI, 28.7 * 2 * PI, 0.03, 10);
//        p.setAxesType(AxisType.LINEAR, AxisType.LINEAR);
        p.setAxesType(AxisType.LOG10, AxisType.LOG10);
        p.setLabels("Spatial Frequency", "Amplitude");
        PlotUtil.showData(name, p, r);
    }


    public void showFS() {
        PlotPanel p1 = new PlotPanel();
        p1.addData(sf, new FunctionStyle("fs", Color.black, 1));
        p1.setRange(0, 0.5, 0, 50);
        p1.setLabels("s (mm)", "f(s)");
        PlotUtil.showData(name, p1);
    }


    public static void main(String[] args) throws Exception {
        ParametersFile pf = new ParametersFile(
            "f:\\data\\2005-09-09-1\\data002\\data002.params");
        double[] freqs = pf.getArrayCell(pf.getIDList()[0], "reversingFrequencies");
        for (int i = 0; i < freqs.length; i++) {
            freqs[i] *= 2 * PI;
        }
        double[][] l = getAverageF1F2(pf.getNeuronsInClass("All/OFF/YFast"), pf);
        double[][] p = getAverageF1F2(pf.getNeuronsInClass("All/OFF/Parasol"), pf);
        double[][] m = getAverageF1F2(pf.getNeuronsInClass("All/OFF/Midget"), pf);

        SpatialSummation L = new SpatialSummation(
            "Large", 162 * micr, 32 * micr, 75 * micr, -0.78, 0.0, l[0], l[1], freqs);
        SpatialSummation P = new SpatialSummation(
            "Parasol", 53 * micr, 21 * micr, 89 * micr, -0.75, 0., p[0], p[1], freqs);
//        Sums M = new Sums("Midget", 20 * micr, 19 * micr, 152 * micr, -0.76, 0.4,
//                          m[0], m[1], freqs);

        int x0 = 50;
        int w = 800;
        L.showPlots(new Rectangle(x0, 50, w, 400));
        P.showPlots(new Rectangle(x0 + w, 50, w, 400));
    }


    public static double[][] getAverageF1F2(int[] id, ParametersFile pFile) throws
        IOException {
        final int nFreq = pFile.getArrayCell(id[0], "reversingFrequencies").length;

        MeanVarianceCalculator[] mvcF1 = new MeanVarianceCalculator[nFreq];
        MeanVarianceCalculator[] mvcF2 = new MeanVarianceCalculator[nFreq];
        for (int periodIndex = 0; periodIndex < nFreq; periodIndex++) {
            mvcF1[periodIndex] = new MeanVarianceCalculator();
            mvcF2[periodIndex] = new MeanVarianceCalculator();
        }

        // calculate overall scaling factor
        double overallsum = 0;
        int nn = 0;
        for (int index = 0; index < id.length; index++) {
            double[] _f1 = pFile.getArrayCell(id[index], "T1reversingF1");
            double[] _f2 = pFile.getArrayCell(id[index], "T1reversingF2");
            if (_f1 != null) {
                overallsum += sum(_f1) + sum(_f2);
                nn++;
            }
        }
        overallsum /= nn;

        for (int index = 0; index < id.length; index++) {
            double[] _f1 = pFile.getArrayCell(id[index], "T1reversingF1");
            double[] _f2 = pFile.getArrayCell(id[index], "T1reversingF2");

            if (_f1 != null) {
                double sum = sum(_f1) + sum(_f2);
                if (sum > 0.0) {
                    multiply(_f1, overallsum / sum);
                    multiply(_f2, overallsum / sum);
                }

                for (int periodIndex = 0; periodIndex < nFreq; periodIndex++) {
                    mvcF1[periodIndex].add(_f1[periodIndex]);
                    mvcF2[periodIndex].add(_f2[periodIndex]);
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
        }

        return new double[][] {f1, f2};
    }

}
