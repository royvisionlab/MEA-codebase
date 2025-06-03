package edu.ucsc.neurobiology.vision.analysis;

import edu.ucsc.neurobiology.vision.plot.*;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class CrossCorrelationCalculator {

    private CrossCorrelationCalculator() {
    }


    /**
     * Calculates the Cross-Correlation histogram of two spike trains.
     *
     * @param t1 sample times of the first spike train
     * @param t2 sample times of the second spike train
     * @param binning time binning in samples
     * @param deltaT time range (-deltaT : +deltaT) of the histogrsm in samples.
     * @return the Cross-Correlation histogram
     */
    public static DoubleHistogram getCrossCorrelationHistogram(
        int[] t1, int[] t2, int binning, int deltaT, DoubleHistogram h) {

        if (t1 == null || t2 == null) {
            throw new NullPointerException("t1 or t2 is null");
        }

        if (t1.length == 0) {
            throw new IllegalArgumentException("t1 has length 0");
        }

        if (t2.length == 0) {
            throw new IllegalArgumentException("t2 has length 0");
        }

        if (h == null) {
            h = new DoubleHistogram("", -deltaT, deltaT, binning);
        } else {
            h.clear();
        }

//        for (int i = 0; i < t1.length; i++) {
//            for (int j = 0; j < t2.length; j++) {
//                h.fill(t1[i] - t2[j], 1);
//            }
//        }

        int j0 = 0;
        int N = t2.length - 1;
        for (int i = 0; i < t1.length; i++) {
            int j = j0;
            int n1 = t1[i] - deltaT;
            int n2 = t1[i] + deltaT;

            while (t2[j] < n1 && j < N) {
                j++;
            }

            j0 = j - 1;
            if (j0 < 0) {
                j0 = 0;
            }

            while (t2[j] < n2 && j < N) {
                h.fill(t2[j] - t1[i], 1);
                j++;
            }
        }

        return h;
    }


    public static DoubleHistogram getCrossCorrelationHistogram(
        double[] t1, double[] t2, int binning, int deltaT) {

        return getCrossCorrelationHistogram(t1, t1.length, t2, t2.length, binning, deltaT);
    }


    public static DoubleHistogram getCrossCorrelationHistogram(
        double[] t1, final int n1, double[] t2, final int n2, int binning, int deltaT) {

        DoubleHistogram h = new DoubleHistogram("", -deltaT, deltaT, binning);

        for (int i = 0, j0 = 0; i < n1; i++) {
            int j = j0;
            while (t2[j] < t1[i] - deltaT && j < n2 - 1) {
                j++;
            }

            j0 = j - 1;
            if (j0 < 0) {
                j0 = 0;
            }

            while (t2[j] < t1[i] + deltaT & j < n2 - 1) {
                h.fill(t2[j] - t1[i], 1);
                j++;
            }
        }

        return h;
    }


}
