package edu.ucsc.neurobiology.vision.math;

import edu.ucsc.neurobiology.vision.plot.*;


/**
 * A spline that can be used even if the points on the x axis are NOT evenly spaced.
 * Taken from Numerical Recipes in Fortran.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class Spline
    implements FunctionData {

    private final int n;
    private float[] x, y;
    private double[] secondDerivatives, u;


    public Spline(int n) {
        this.n = n;
        x = new float[n];
        y = new float[n];
        u = new double[n];
        secondDerivatives = new double[n];
    }


    public void reSpline(double xArray[], double yArray[]) {
        if (xArray.length != yArray.length) {
            throw new IllegalArgumentException(
                "x.length = " + xArray.length + ", y.length = " + yArray.length);
        }
        if (xArray.length != n) {
            throw new IllegalArgumentException("x.length != n");
        }

        for (int i = 0; i < n; i++) {
            x[i] = (float) xArray[i];
            y[i] = (float) yArray[i];
        }
        spline();
    }


    public void reSpline(float xArray[], float yArray[]) {
        if (xArray.length != yArray.length) {
            throw new IllegalArgumentException(
                "x.length = " + xArray.length + ", y.length = " + yArray.length);
        }
        if (xArray.length != n) {
            throw new IllegalArgumentException("x.length != n");
        }

        System.arraycopy(xArray, 0, x, 0, n);
        System.arraycopy(yArray, 0, y, 0, n);
        spline();
    }


    public void reSpline(float xArray[], short yArray[]) {
        if (xArray.length != yArray.length) {
            throw new IllegalArgumentException(
                "x.length = " + xArray.length + ", y.length = " + yArray.length);
        }
        if (xArray.length != n) {
            throw new IllegalArgumentException("x.length != n");
        }

        System.arraycopy(xArray, 0, x, 0, n);
        for (int i = 0; i < n; i++) {
            y[i] = yArray[i];
        }

        spline();
    }


    public String getDescription() {
        return "spline";
    }


    private void spline() {
        double p, qn, sig, un;

        for (int i = 0; i < n; i++) {
            secondDerivatives[i] = 0;
            u[i] = 0;
        }

        for (int i = 1; i < n - 1; i++) {
            sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
            p = sig * secondDerivatives[i - 1] + 2;
            secondDerivatives[i] = (sig - 1) / p;
            u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) -
                   (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
            u[i] = (6 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
        }

        // To match NR online, the below should read
        // 	 (un - qn * u[n-2]) / (qn * secondDerivatives[n-2] + 1)
        // But of course with qn = un = 0 this whole section is superfluous
        qn = un = 0;
        secondDerivatives[n -
            1] = (un - qn * u[n - 1]) / (qn * secondDerivatives[n - 1] + 1);

        for (int k = n - 2; k >= 0; k--) {
            secondDerivatives[k] = secondDerivatives[k] * secondDerivatives[k + 1] + u[k];
        }
    }


    public final double getValueAt(double coord) {
        int j1 = 0;
        int j2 = n - 1;
        while (j2 - j1 > 1) {
            int k = (j2 + j1) >> 1;
            if (x[k] > coord) {
                j2 = k;
            } else {
                j1 = k;
            }
        }

//        System.out.println(j1 + " : " + j2);

        final double h = x[j2] - x[j1];
        if (h == 0) {
            throw new IllegalArgumentException("Bad x input to routine splint");
        }
        final double a = (x[j2] - coord) / h;
        final double b = (coord - x[j1]) / h;

        return a * y[j1] + b * y[j2] + ( (h * h) / 6) * (
            (a * a * a - a) * secondDerivatives[j1] +
            (b * b * b - b) * secondDerivatives[j2]);
    }


}
