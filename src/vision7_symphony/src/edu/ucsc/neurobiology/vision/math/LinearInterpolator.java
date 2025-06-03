package edu.ucsc.neurobiology.vision.math;

import edu.ucsc.neurobiology.vision.plot.*;


/**
 * Linearly interpolates a set of data points. Acts as a plottable function.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class LinearInterpolator
    implements FunctionData {

    private final int n;
    private float[] x, y;
    private double[] k, b;


    public LinearInterpolator(int n) {
        this.n = n;
        x = new float[n];
        y = new float[n];
        k = new double[n];
        b = new double[n];
    }


    public void interpolate(double xArray[], double yArray[]) {
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
        recalculateCoefficients();
    }


    public void interpolate(float xArray[], float yArray[]) {
        if (xArray.length != yArray.length) {
            throw new IllegalArgumentException(
                "x.length = " + xArray.length + ", y.length = " + yArray.length);
        }
        if (xArray.length != n) {
            throw new IllegalArgumentException("x.length != n");
        }

        System.arraycopy(xArray, 0, x, 0, n);
        System.arraycopy(yArray, 0, y, 0, n);
        recalculateCoefficients();
    }


    public void interpolate(float xArray[], short yArray[]) {
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

        recalculateCoefficients();
    }


    public String getDescription() {
        return "spline";
    }


    private void recalculateCoefficients() {
        for (int i = 0; i < n; i++) {
            k[i] = 0;
            b[i] = 0;
        }

        for (int i = 0; i < n - 1; i++) {
            k[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
            b[i] = y[i] - k[i] * x[i];
        }
    }


    public final double getValueAt(double coord) {
        // find the segment
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
//        System.out.println("x " + x[j1] + " : " + x[j2]);
//        System.out.println("y " + y[j1] + " : " + y[j2]);
//        System.out.println("kb " + k[j1] + " : " + b[j2]);

        return k[j1] * coord + b[j1];
    }

}
