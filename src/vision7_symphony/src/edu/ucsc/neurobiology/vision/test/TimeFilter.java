package edu.ucsc.neurobiology.vision.test;

import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.math.fitting.*;
import edu.ucsc.neurobiology.vision.plot.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class TimeFilter
    extends Function implements FunctionMinimization.IterationInfo {

    double[] tc, err, p;
    double T = 8.34 / 1000;
    int N;
    ScatterPlotStyle s = new ScatterPlotStyle();
    boolean debug = false;


    public TimeFilter(double[] _tc, double[] _err) {
        this.tc = new double[64];
        this.err = new double[64];
        for (int j = 0; j < _tc.length; j++) {
            tc[j] = _tc[_tc.length - j - 1];
            err[j] = _err[_tc.length - j - 1];
        }
        for (int j = _tc.length; j < tc.length; j++) {
            tc[j] = _tc[0];
            err[j] = _err[0];
        }
        N = tc.length;
        System.out.println("N " + N);
        s.setConnectingPoints(true);
    }


    /*
        public static double[] getFilterFFT(
            int N, double T, double _tL, double nL,
                                            double _tH,
                                            double nH, double aH, double g) {
            Complex i = new Complex(0, 1);
            Complex one = new Complex(1, 0);
            double dw = 2 * Math.PI / (N * T);
            double[] re = new double[N];
            double[] im = new double[N];
            for (int j = 0; j < N; j++) {
                double w = j * dw;
                Complex lp = Complex.pow(
                    one.over(one.plus(i.times(w * _tH))), nL);
                Complex hp = Complex.pow(
                    one.minus(one.times(aH).over(one.plus(i.times(w * _tL)))), nH);
                Complex f = Complex.times(lp, hp).times(g);
                re[j] = f.real();
                im[j] = f.imag();
            }
            FFT fft = new FFT();
            fft.fft(re, im, -1);
            VisionUtilities.divide(re, N);
            return re;
        }
     */

    public static double[] getFilter(
        int N, double T, double tL, int nL, double tH, int _nH, double aH, double g) {

        double zL = -1 / tL;
        double zH = -1 / tH;
//        System.out.println("zL " + zL);
//        System.out.println("zH " + zH);
//        System.out.println("n " + nL);

        double[] f = new double[N];
        double g1 = g * aH * zL * Math.pow( -zH / (zL - zH), nL);
        double g2 = g * Math.pow( -zH, nL) / MathUtil.fact(nL - 1);
        double g3 = g2 * aH * zL;
//        System.out.println("g1 " + g1);
//        System.out.println("g2 " + g2);
//        System.out.println("g3 " + g3);

        for (int j = 0; j < N; j++) {
            double t = j * T;
            double d = 0;
            for (int i = 0; i <= nL - 1; i++) {
                d += (MathUtil.fact(nL - 1) / MathUtil.fact(i)) *
                    Math.pow( -1, i) * Math.pow(t, i) / Math.pow(zH - zL, nL - i);
            }
            d *= Math.exp(zH * t);
            f[j] = g1 * Math.exp(zL * t) +
                   g2 * Math.pow(t, nL - 1) * Math.exp(zH * t) +
                   g3 * d;
            f[j] *= T;
        }

        return f;
    }


    public double getValue(double[] p) {
        this.p = p;
        double g = p[0];
        double tL = p[1] / 1000;
        double nL = p[2];
        double tH = p[3] / 1000;
        double aH = p[4];
        double nH = 1;

        double[] f = getFilter(N, T, tL, (int) nL, tH, (int) nH, aH, g);
//        double[] f = getFilterFFT(N, T, tL, nL, tH, nH, aH, g);
        if (debug) {
            double[] x = new double[N];
            for (int j = 0; j < N; j++) {
                x[j] = j * T;
            }
            PlotPanel pan = new PlotPanel();
            pan.addData(new ScatterPlot(x, tc, err, null), s);
            pan.addData(new ScatterPlot(x, f, null, null), s);
            PlotUtil.showData("Filter", pan);
            pan.autoscale();
        }

        double chi2 = 0;
        for (int j = 0; j < N; j++) {
            chi2 += Math.pow( (tc[j] - f[j]) / err[j], 2);
        }

        return chi2 / N;
    }


    public void receiveInfo(String info) {
        System.out.println(info);
        for (int k = 0; k < p.length; k++) {
            System.out.println(p[k] + ",");
        }
    }


    /*
        public double getA1() {
            return parameters[0];
        }
        public double getA2() {
            return -parameters[0] * parameters[1] / parameters[2];
        }
        public double getT1() {
            return parameters[1];
        }
        public double getT2() {
            return parameters[2];
        }
     */

    public String getDescription() {
        return "STATimeFunction";
    }

}
