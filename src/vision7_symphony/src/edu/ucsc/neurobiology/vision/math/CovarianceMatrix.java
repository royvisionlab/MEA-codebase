package edu.ucsc.neurobiology.vision.math;

import java.util.*;


/**
 * This class implements a covariance matrix. The matrix is stored in reduced form.
 * Data samples have to added by any addData() method.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class CovarianceMatrix {
    private double[] covariances;
    private double[] means;
    private final int size;
    public int n;
    private boolean normalized = false;


    /**
     * A covariance matrix in reduced form.
     *
     * @param size
     */
    public CovarianceMatrix(int size) {
        this.size = size;
        covariances = new double[size * (size + 1) / 2];
        means = new double[size];
    }


    public void reset() {
        Arrays.fill(covariances, 0);
        Arrays.fill(means, 0);
        n = 0;
    }


    public final void addData(final float[] x) {
        int i, j, index = 0;

        for (i = 0; i < size; i++) {
            means[i] += x[i];
            for (j = i; j < size; j++) {
                covariances[index] += x[i] * x[j];
                index++;
            }
        }

        n++;
    }


    public final void addData(final byte[] x) {
        int i, j, index = 0;

        for (i = 0; i < size; i++) {
            means[i] += x[i];
            for (j = i; j < size; j++) {
                covariances[index] += x[i] * x[j];
                index++;
            }
        }

        n++;
    }


    public final void addData(final short[] x) {
        int i, j, index = 0;

        for (i = 0; i < size; i++) {
            means[i] += x[i];
            for (j = i; j < size; j++) {
                covariances[index] += x[i] * x[j];
                index++;
            }
        }

        n++;
    }


    public final void addData(final double[] x) {
        int i, j, index = 0;

        for (i = 0; i < size; i++) {
            means[i] += x[i];
            for (j = i; j < size; j++) {
                covariances[index] += x[i] * x[j];
                index++;
            }
        }

        n++;
    }


    public void combine(CovarianceMatrix m) {
        this.n += m.n;

        for (int i = 0; i < means.length; i++) {
            this.means[i] += m.means[i];
        }

        for (int i = 0; i < covariances.length; i++) {
            this.covariances[i] += m.covariances[i];
        }
    }


    private void normalize() {
        if (n < 2) {
            return;
        }

        for (int i = 0, index = 0; i < size; i++) {
            for (int j = i; j < size; j++) {
                covariances[index] = (float) ( (covariances[index] -
                                               means[i] * means[j] / n) / (n - 1));
                index++;
            }
        }

        for (int i = 0; i < size; i++) {
            means[i] /= n;
        }
    }


    public float[] getCovariance() {
        if (!normalized) {
            normalize();
            normalized = true;
        }
        
        float[] floatCov = new float[covariances.length];
        for (int i = 0; i < covariances.length; i++) {
            floatCov[i] = (float) covariances[i];
        }
 

        return floatCov;
    }


    public float[] getSigmas() {
        if (!normalized) {
            normalize();
            normalized = true;
        }

        float[] sigmas = new float[size];
        for (int i = 0, index = 0; i < size; i++) {
            sigmas[i] = (float) Math.sqrt(covariances[index]);
            index += size - i;
        }
        return sigmas;
    }


    public float[] getMeans() {
        if (!normalized) {
            normalize();
            normalized = true;
        }
        
        float[] floatMeans = new float[means.length];
        for(int i=0; i<means.length; i++) {
            floatMeans[i] = (float) means[i];
        }

        return floatMeans;
    }


    public int getNPoints() {
        return n;
    }
}
