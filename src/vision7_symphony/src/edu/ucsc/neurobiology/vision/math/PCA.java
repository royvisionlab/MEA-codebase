package edu.ucsc.neurobiology.vision.math;

import cern.colt.matrix.impl.*;
import cern.colt.matrix.linalg.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * Principal Component Analysis (PCA) implementation. Works with both complete
 * and reduced covariance matrices.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class PCA {
    private double[][] covariance;
    private double[] d, e;
    private int size;


    /**
     * Initializes PCA with a complete form matrix.
     *
     * @param covariance a complete matrix
     */
    public PCA(double[][] covariance) {
        this.covariance = covariance;
        this.size = covariance.length;
        d = new double[size];
        e = new double[size];
    }


    public double getDeterminant() {
        DenseDoubleMatrix2D m = new DenseDoubleMatrix2D(covariance);
        return Algebra.DEFAULT.trace(m);
    }


    /**
     * Initializes PCA with a reduced form matrix.
     *
     * @param cov
     */
    public PCA(float[] cov) {
        this.size = (int) Math.round( (Math.sqrt(1 + 8 * cov.length) - 1) / 2);
        if (size * (size + 1) / 2 != cov.length) {
            throw new IllegalArgumentException("Wrong length for covariance");
        }
        this.covariance = new double[size][size];
        int index = 0;
        for (int i = 0; i < size; i++) {
            for (int j = i; j < size; j++) {
                covariance[i][j] = covariance[j][i] = cov[index];
                index++;
            }
        }

        d = new double[size];
        e = new double[size];
    }


    public synchronized void doPCA() throws TooManyIterationsException {
        //normalise();
        tred2(covariance, size, d, e);        
        tqli(d, e, size, covariance);

        // exchange raws with columns
        double temp;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < i; j++) {
                temp = covariance[i][j];
                covariance[i][j] = covariance[j][i];
                covariance[j][i] = temp;
            }
        }
        sortEigenValues();
    }


    public double getEigenvaluesPercentage(int n) {
        double sum = 0;
        for (int i = 0; i < size; i++) {
            sum += d[i];
        }

        double nSum = 0;
        for (int i = 0; i < n; i++) {
            double p = 100 * d[i] / sum;
            nSum += p;
        }

        return nSum;
    }


    public void showEigenVectors(int n) {
        if (n == -1) {
            n = size;

        }
        for (int i = 0; i < n; i++) {
            ScatterPlotStyle style = new ScatterPlotStyle();
            style.setSymbolType(SymbolType.NONE);
            style.setConnectingPoints(true);

            PlotPanel p = new PlotPanel();
            p.addData(new ScatterPlot(getEigenVector(i)), style);
            p.setLabels("Index", "Value");
            PlotUtil.showData("Eigenvector " + (i + 1), p, 400, 200);
            p.autoscale();
        }

        ScatterPlot eigSP = new ScatterPlot();
        for (int i = 0; i < n; i++) {
            eigSP.add(i + 1, getEigenValue(i));
        }
        PlotPanel p = new PlotPanel();
        ScatterPlotStyle style = new ScatterPlotStyle();
        style.setSymbolSize(10);
        p.addData(eigSP, style);
        p.setLabels("Index", "Eigenvalue (%)");
        PlotUtil.showData("", p);
        p.setRange(0.9, 5.1, 0, 100);
    }


    public void printPercentageEigenvalues(int n) {
        System.out.print("EigenValues: ");
        double sum = 0;
        for (int i = 0; i < size; i++) {
            sum += d[i];
        }

        double nSum = 0;
        for (int i = 0; i < n; i++) {
            double p = 100 * d[i] / sum;
            nSum += p;
            System.out.print(StringUtil.format(p, 0) + ", ");
        }

        System.out.println("Sum: " + StringUtil.format(nSum, 1) + " %");
    }


    // bubble sort
    private void sortEigenValues() {
//        d[j] = d[j - inc];
//        covariance[j] = covariance[j - inc];
        for (int i = 0; i < size; i++) {
            int max = i;
            for (int j = i + 1; j < size; j++) {
                if (d[j] > d[max]) {
                    max = j;
                }
            }

            double temp = d[i];
            d[i] = d[max];
            d[max] = temp;

            double[] temp2 = covariance[i];
            covariance[i] = covariance[max];
            covariance[max] = temp2;
        }
    }


    // The SHELL sorting method from "Numerical Recipes (Fortran)" pp. 229
//    private void sortEigenValues() {
//        int i, j, inc = 1;
//        double v;
//        double[] vector;
//
//        do {
//            inc *= 3;
//            inc++;
//        } while (inc < size); //} while (inc <= n);
//
//        do {
//            inc /= 3;
//            for (i = inc + 1; i < size; i++) { //for (i = inc + 1; i <= n; i++) {
//                v = d[i];
//                vector = covariance[i];
//                j = i;
//                while (d[j - inc] > v) {
//                    d[j] = d[j - inc];
//                    covariance[j] = covariance[j - inc];
//                    j -= inc;
//                    if (j <= inc) break;
//                }
//                d[j] = v;
//                covariance[j] = vector;
//            }
//        } while (inc > 1);
//    }


    public double getEigenValue(int n) {
        return 100 * d[n] / MathUtil.sum(d);
    }


    public double[][] getEigenVectors() {
        return covariance;
    }


    public double[] getEigenVector(int vectorIndex) {
        return covariance[vectorIndex];
    }


    public final double project(double[] x, int vectorIndex) {
        double v = 0;

        for (int i = 0; i < size; i++) {
            v += x[i] * covariance[vectorIndex][i];
        }

        return v;
    }


    public final double project(float[] x, int vectorIndex) {
        double v = 0;

        for (int i = 0; i < size; i++) {
            v += x[i] * covariance[vectorIndex][i];
        }

        return v;
    }


    public final double project(byte[] x, int vectorIndex) {
        double v = 0;

        for (int i = 0; i < size; i++) {
            v += x[i] * covariance[vectorIndex][i];
        }

        return v;
    }


    public final double project(short[] x, int vectorIndex) {
        double v = 0;

        for (int i = 0; i < size; i++) {
            v += x[i] * covariance[vectorIndex][i];
        }

        return v;
    }


    /**
     * The Householder tridiagonal reduction method from
     * "Numerical Recipes (Fortran)" pp. 354
     */
    private static void tred2(double[][] a, int n, double d[], double e[]) {
        int l, k, j, i;
        double scale, hh, h, g, f;

        for (i = n - 1; i >= 1; i--) {
            l = i - 1;
            h = scale = 0;
            if (l > 0) {
                for (k = 0; k <= l; k++) {
                    scale += Math.abs(a[i][k]);
                }

                if (scale == 0) {
                    e[i] = a[i][l];
                } else {
                    for (k = 0; k <= l; k++) {
                        a[i][k] /= scale;
                        h += a[i][k] * a[i][k];
                    }
                    f = a[i][l];
                    g = (f >= 0) ? - (double) Math.sqrt(h) : (double) Math.sqrt(h);
                    e[i] = scale * g;
                    h -= f * g;
                    a[i][l] = f - g;
                    f = 0;
                    for (j = 0; j <= l; j++) {
                        a[j][i] = a[i][j] / h;
                        g = 0;
                        for (k = 0; k <= j; k++) {
                            g += a[j][k] * a[i][k];
                        }
                        for (k = j + 1; k <= l; k++) {
                            g += a[k][j] * a[i][k];
                        }
                        e[j] = g / h;
                        f += e[j] * a[i][j];
                    }
                    hh = f / (h + h);
                    for (j = 0; j <= l; j++) {
                        f = a[i][j];
                        e[j] = g = e[j] - hh * f;
                        for (k = 0; k <= j; k++) {
                            a[j][k] -= (f * e[k] + g * a[i][k]);
                        }
                    }
                }
            } else {
                e[i] = a[i][l];
            }
            d[i] = h;
        }
        d[0] = 0;
        e[0] = 0;

        /* Contents of this loop can be omitted if eigenvectors not
                 wanted except for statement d[i] = a[i][i]; */
        for (i = 0; i < n; i++) {
            l = i - 1;
            if (d[i] != 0) {
                for (j = 0; j <= l; j++) {
                    g = 0;
                    for (k = 0; k <= l; k++) {
                        g += a[i][k] * a[k][j];
                    }
                    for (k = 0; k <= l; k++) {
                        a[k][j] -= g * a[k][i];
                    }
                }
            }
            d[i] = a[i][i];
            a[i][i] = 1;
            for (j = 0; j <= l; j++) {
                a[j][i] = a[i][j] = 0;
            }
        }
    }
    
    
    // The QL algorithm with implicit shifts from "Numerical Recipes (Fortran)" pp. 362
    private static void tqli(double[] d, double[] e, int n, double[][] z) throws
        TooManyIterationsException {

        int m, l, iter, i, k;
        double s, r, p, g, f, dd, c, b;

        for (i = 1; i < n; i++) {
            e[i - 1] = e[i];
        }
        e[n - 1] = 0;

        for (l = 0; l < n; l++) {
            iter = 0;
            do {
                for (m = l; m < n - 1; m++) {
                    dd = Math.abs(d[m]) + Math.abs(d[m + 1]);
                    if (Math.abs(e[m]) + dd == dd) {
                        break;
                    }
                }
                if (m != l) {
                    iter++;
                    g = (d[l + 1] - d[l]) / (2 * e[l]);
                    r = (double) Math.sqrt(1 + g * g);
                    if (g >= 0) {
                        g += Math.abs(r);
                    } else {
                        g -= Math.abs(r);
                    }
                    g = d[m] - d[l] + e[l] / g;
                    s = c = 1;
                    p = 0;
                    for (i = m - 1; i >= l; i--) {
                        f = s * e[i];
                        b = c * e[i];
                        e[i + 1] = r = Math.sqrt(f * f + g * g);
                        if (r == 0) {
                            d[i + 1] -= p;
                            e[m] = 0;
                            break;
                        }
                        s = f / r;
                        c = g / r;
                        g = d[i + 1] - p;
                        r = (d[i] - g) * s + 2 * c * b;
                        d[i + 1] = g + (p = s * r);
                        g = c * r - b;
                        for (k = 0; k < n; k++) {
                            f = z[k][i + 1];
                            z[k][i + 1] = s * z[k][i] + c * f;
                            z[k][i] = c * z[k][i] - s * f;
                        }
                    }
                    if (r == 0. && i >= l) {
                        continue;
                    }
                    d[l] -= p;
                    e[l] = g;
                    e[m] = 0;
                }
                if (iter > 5000) {
                    int a = 0;
                    throw new TooManyIterationsException("Too many iterations in tqli.");
                }
            } while (m != l);
        }
    }

}
