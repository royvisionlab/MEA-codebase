package edu.ucsc.neurobiology.vision.math;


/**
 * A class that calculates the mean and variance of a set of vectors added by add().
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ArrayMeanVarianceCalculator {
//    public static final int BIASED = 0;
//    public static final int UNBIASED = 0;

    double[] sumX;
    double[] sumXSquare;
    int n;
    final int length;


    public ArrayMeanVarianceCalculator(int length) {
        this.length = length;
        reset();
    }


    public void reset() {
        sumX = new double[length];
        sumXSquare = new double[length];
        n = 0;
    }


    public void add(double[] x) {
        for (int i = 0; i < x.length; i++) {
            sumX[i] += x[i];
            sumXSquare[i] += x[i] * x[i];
        }
        n++;
    }


    public double[] getMean() {
        double[] mean = new double[length];
        for (int i = 0; i < length; i++) {
            mean[i] = sumX[i] / n;
        }
        return mean;
    }


    public double[] getVariance() {
        double[] v = new double[length];

        if (n == 1) {
            return v;
        } else {
            for (int i = 0; i < length; i++) {
                v[i] = Math.sqrt( (sumXSquare[i] - sumX[i] * sumX[i] / n) / (n - 1));
            }
            return v;
        }
    }


    public double[] getMeanVariance() {
        double[] v = getVariance();
        for (int i = 0; i < length; i++) {
            v[i] /= Math.sqrt(n);
        }
        return v;
    }
}
