package edu.ucsc.neurobiology.vision.testing;

import static java.lang.Math.*;
import java.util.*;

import cern.colt.matrix.*;
import cern.colt.matrix.impl.*;
import cern.colt.matrix.linalg.*;
import edu.ucsc.neurobiology.vision.io.*;
import static edu.ucsc.neurobiology.vision.math.Matrix.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.math.expressions.*;
import edu.ucsc.neurobiology.vision.math.fitting.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.test.*;


/**
 * @author nobody, anyone can change
 */
class MathTest {

    static SparseDoubleMatrix2D m1 = new SparseDoubleMatrix2D(5, 5);
    static SparseDoubleMatrix2D m2 = new SparseDoubleMatrix2D(new double[][] { {1}
        , {1}
        , {1}
        , {1}
        , {1}
    });


    public static void testMatrix() {
        for (int n = 0; n < 1000; n++) {
            double ai = Math.random() * 10;
            double bi = Math.random() * 10;
            double alphai = Math.random() * 2 * PI;
//
//        double ai = 9.286020284471011;
//        double bi = 8.524991402209688;
//        double alphai = 72.68286519373531 * PI / 180;
            ParametricEllipse ei = new ParametricEllipse(0, 0, ai, bi, alphai);

            for (int i = 0; i < 5; i++) {
                double theta;
                if (i == 4) {
                    theta = PI / 6;
                } else {
                    theta = i * PI / 2;
                }
                double[] p = ei.getPointFor(theta);
                double x = p[0];
                double y = p[1];
                m1.set(i, 0, x * x);
                m1.set(i, 1, x * y);
                m1.set(i, 2, y * y);
                m1.set(i, 3, x);
                m1.set(i, 4, y);
            }

            // solve the system and find the polynomial coefficients
            DoubleMatrix2D s = Algebra.DEFAULT.solve(m1, m2);
            double A = s.get(0, 0);
            double B = s.get(1, 0);
            double C = s.get(2, 0);
            double D = s.get(3, 0);
            double E = s.get(4, 0);

            // find the ellipse parameters
            double t = 1 / (B * B - 4 * A * C);
            double M = (C * D * D - B * D * E + A * E * E) * t;
            double x0 = (2 * D * C - B * E) * t;
            double y0 = (2 * A * E - B * D) * t;
            double w = Math.sqrt(B * B + (A - C) * (A - C));
            double a = 1 / Math.sqrt(0.5 * (A + C + w) / (1 - M));
            double b = 1 / Math.sqrt(0.5 * (A + C - w) / (1 - M));
            double alpha = 0.5 * Math.atan2(B, (A - C));

            while (alpha < 0) {
                alpha += PI;
            }

            // do corrections
            if (!correct(ei, new ParametricEllipse(x0, y0, a, b, alpha))) {
                for (int i = 0; i < 4; i++) {
                    double da = i * PI / 2;
                    if (correct(ei, new ParametricEllipse(x0, y0, a, b, alpha + da))) {
//                    System.out.println("worked: +" + da * 180 / PI);
                        alpha += da;
                        break;
                    }
                    if (correct(ei, new ParametricEllipse(x0, y0, b, a, alpha + da))) {
//                    System.out.println("worked: inversion +" + da * 180 / PI);
                        alpha += da;
                        double temp = a;
                        a = b;
                        b = temp;
                        break;
                    }
                }
            }

            while (alpha > 2 * PI) {
                alpha -= 2 * PI;
            }

            if (Math.abs(ai - a) < err && Math.abs(bi - b) < err &&
                Math.abs(alphai - alpha) < err) {
            } else {
                System.out.println("I: " + ai + ", " + bi + ": " + alphai * 180 / PI);
                System.out.println("F: " + a + ", " + b + ": " + alpha * 180 / PI);
            }
            /*
                    ScatterPlotStyle st = new ScatterPlotStyle();
                    st.setSymbolSize(3);
                    PlotPanel pi = new PlotPanel();
                    pi.addData(ei, new FunctionStyle());
                    pi.addData(ei.getCardinalPoints(), st);
                    pi.setRange( -20, 20, -20, 20);
                    PlotPanel p = new PlotPanel();
                    ParametricEllipse e = new ParametricEllipse( x0, y0, a, b, alpha);
                    p.addData(e, new FunctionStyle());
                    p.addData(e.getCardinalPoints(), st);
                    p.setRange( -20, 20, -20, 20);
                    JFrame f = new JFrame();
                    f.getContentPane().setLayout(new GridLayout(2, 1));
                    f.add(pi);
                    f.add(p);
                    f.setBounds(100, 100, 500, 500 + 19);
                    f.setVisible(true);
             */
        }
    }


    public static void testSpline() {
        Spline s = new Spline(100);

        double[] x = new double[100];
        double[] y = new double[100];
        for (int i = 0; i < x.length; i += 1) {
            x[i] = i;
            y[i] = Math.sin(i / 2.0) * Math.exp( -0.05 * i);
        }
        s.reSpline(x, y);

        PlotPanel p = new PlotPanel();
        p.addData(s, new FunctionStyle("spline"));

        p.addData(new ScatterPlot(x, y, null), "DISK 5 red");
        p.autoscale();
        PlotUtil.showData("", p);
    }


    public static void testLinear1DFunction() throws Exception {
        double[][] x = { {1, 2, 3, 4}
        };
        double[] y = {2.99, 4.89, 5, 8};
        double[] sig = {1, 1, 1, 1};

        Linear1DFunction l = Linear1DFunction.lineFit(x[0], y, sig);
        System.out.println(l);

        Linear1DFunction f = new Linear1DFunction();
        f.setParameters(new double[] {0.5, 0.5});

        Fitter fit = new Fitter();
        fit.fit(f, x, y, sig);

        System.out.println("\n" + f);
    }


    public static void testUniformSpline() {
        UniformSpline s = new UniformSpline(100);
        double[] y = new double[100];
        for (int i = 0; i < y.length; i++) {
            y[i] = i * i;
        }
        s.reSpline(y);

        double t1 = System.currentTimeMillis();
        for (double i = 1; i < y.length - 1; i += 1e-6) {
            s.getValueAt(i);
        }
        double t2 = System.currentTimeMillis();
        System.out.println( (t2 - t1) / 1000.0);
    }


    public static void testRandomMersenne() {
        RandomMersenne rand = new RandomMersenne();
        rand.setSeed(11111);
        DoubleHistogram h = new DoubleHistogram("", -2 * Math.PI,
                                                Math.PI * 2, Math.PI / 1000.0);
        double min = Double.POSITIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < 10000000; i++) {
            double bin = rand.nextGaussian(0, 1);

            min = (min < bin) ? min : bin;
            max = (max > bin) ? max : bin;
            h.fill(bin, 1.0);
        }
        System.out.println("Max: " + max);
        System.out.println("Min: " + min);
        System.out.println("PI: " + Math.PI);

        PlotPanel p = new PlotPanel();
        p.addData(h, new HistogramStyle());
        p.autoscale();
        PlotUtil.showData("", p);
    }


    public static void testArbitraryFunction() throws Exception {
        double[] x = new double[100];
        double[] y = new double[x.length];
        double[] err = new double[x.length];
        {
            Random r = new Random(111);
            FunctionSum f = new FunctionSum(new Gaussian1DFunction(4, -1, 0.2),
                                            new Gaussian1DFunction( -2, -1, 0.5));
            for (int i = 0; i < x.length; i++) {
                err[i] = 0.3;
                x[i] = -3 + i * 0.05;
                y[i] = f.getValueAt(x[i]) + (r.nextDouble() - 0.5) * err[i];
            }
        }

        // ads
        ArbitraryFunction ff = new ArbitraryFunction(
            "A1*exp(-((x-x0)^2)/(2*s1^2)) + A2*exp(-((x-x0)^2)/(2*s2^2))",
            "x",
            "A1 A2 s1 s2 x0",
            4, -1, 0.2, 0.5, -.75);

        Fitter.fit1D(ff, x, y, err, x.length);
        System.err.println(ff);

        // show
        PlotPanel p = new PlotPanel();
        p.addData(new ScatterPlot(x, y, err), new ScatterPlotStyle());
        p.addData(ff, new FunctionStyle(""));
        p.setRange( -3, 2, -5, 5);
        PlotUtil.showData("", p);
    }


    public static void testFFT() throws Exception {
        int N = 1024 * 2;
        int n = 1000;

        double[] re = new double[N];
        double[] im = new double[N];
        double[] pow = new double[N];
        for (int i = 0; i < n; i++) {
            re[i] = Math.sin(i);
        }

        FFT.fft(re, im, -1);

        for (int i = 0; i < pow.length; i++) {
            pow[i] = (2.0 / n) * Math.sqrt(re[i] * re[i] + im[i] * im[i]);
        }

        ScatterPlotStyle s = new ScatterPlotStyle();
        s.setConnectingPoints(true);
        s.setSymbolType(SymbolType.NONE);
        PlotUtil.showData(new ScatterPlot(pow), s);
    }


    public static void testExpession() throws Exception {
        HashMap v = new HashMap();
        v.put("x", new double[] {1, 2, 3, 4, 5});

        Object r = new Expression("(x ^ x)[2]").evaluate(v);
        if (r instanceof double[]) {
            IOUtil.printArray( (double[]) r);
        } else {
            System.err.println(r);
        }
    }
    
    public static void testCovarianceMatrix() throws Exception {
        // size of input vectors
        int size = 4;
        // create the covariance matrix
        CovarianceMatrix cov = new CovarianceMatrix(size);
        // constants to add
        double oA = 0.1; double oB = 0.65; double oC = -12; double oD = 0;
        // add the input vectors
        cov.addData(new double[] {-1+oA,1+oB,2+oC,-5+oD});
        cov.addData(new double[] {-2+oA,3+oB,1+oC,-3+oD});
        cov.addData(new double[] {4+oA,0+oB,3+oC,7+oD});
        // print the output
        float[] output = cov.getCovariance();
        IOUtil.printArray(output);
    }


    public static void main(String[] args) throws Exception {
        testCovarianceMatrix();
    }
}
