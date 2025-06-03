package edu.ucsc.neurobiology.vision.plot;

import java.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class HistogramND {
    protected String description;
    private int nDimensions;
    protected double[] min, max;
    protected double binInterval[];
    private int[] nBins;
    protected short[] bins;
    protected short[] underflow, overflow;
    int[] tempBin;


    public HistogramND(double[] min, double[] max, double[] binInterval) {
        if (min.length != max.length || min.length != binInterval.length) {
            throw new IllegalArgumentException("not the same lebgth");
        }

        this.nDimensions = min.length;
        this.min = min;
        this.max = max;
        this.binInterval = binInterval;

        int totalNBins = 1;
        nBins = new int[nDimensions];
        for (int d = 0; d < nDimensions; d++) {
            nBins[d] = (int) Math.ceil( (max[d] - min[d]) / binInterval[d]);
//            System.out.println(min[d] + ":" + max[d] + ":" + binInterval[d]);
            totalNBins *= nBins[d];
        }

        this.tempBin = new int[nDimensions];

        bins = new short[totalNBins];
        underflow = new short[nDimensions];
        overflow = new short[nDimensions];
    }


    public int countNonZero() {
        int n = 0;
        for (int i = 0; i < bins.length; i++) {
            if (bins[i] != 0) {
                n++;
            }
        }
        return n;
    }


    public double getUnderflow() {
        double u = 0;

        for (int d = 0; d < nDimensions; d++) {
            u += underflow[d];
        }

        return u;
    }


    public double getOverflow() {
        double o = 0;

        for (int d = 0; d < nDimensions; d++) {
            o += overflow[d];
        }

        return o;
    }


    public static void minimizeSpace(
        float[][] data, int nPoints, double[] min, double[] max, double[] binInterval,
        int maxLostPoints) {
        /*
                    System.out.println("nPoints " +nPoints);
                    for (int d = 0; d < min.length; d++) {
                        System.out.println(min[d] + ":" + max[d] + ":" + binInterval[d]);
                    }
                    System.out.println("-------------------");
         */

        if (nPoints <= maxLostPoints) {
            return;
        }
        int nDimensions = min.length;

        DoubleHistogram[] h = new DoubleHistogram[nDimensions];
        int[][] limits = new int[nDimensions][2];
        for (int d = 0; d < nDimensions; d++) {
            h[d] = new DoubleHistogram("", min[d], max[d], binInterval[d]);
            for (int i = 0; i < nPoints; i++) {
                h[d].fill(data[d][i], 1);
            }
            limits[d][0] = 0;
            limits[d][1] = h[d].getBinCount() - 1;
//            PlotUtilities.showData("" + d, h[d], new HistogramStyle());
        }

        int nLostSpikes = 0;

        while (nLostSpikes < maxLostPoints) {
            int minLossDimension = -1;
            int minLoss = Integer.MAX_VALUE;
            int minLossSide = -1;

            for (int d = 0; d < nDimensions; d++) {
                for (int side = 0; side < 2; side++) {
                    int n = (int) h[d].getBin(limits[d][side]);
                    if (n < minLoss) {
                        minLoss = n;
                        minLossDimension = d;
                        minLossSide = side;
                    }
                }
            }

            if (minLossSide == 0) {
                // left side
                limits[minLossDimension][minLossSide]++;
            } else {
                // right side
                limits[minLossDimension][minLossSide]--;
            }

            nLostSpikes += minLoss;
        }

        for (int d = 0; d < nDimensions; d++) {
            double minOld = min[d];
            double maxOld = max[d];
            min[d] += limits[d][0] * binInterval[d];
            max[d] -= (h[d].getBinCount() - limits[d][1] - 1) * binInterval[d];

            // Salk, Aug 2005, bug fix, the data is must be very funny;
            if (max[d] <= min[d]) {
                min[d] = minOld;
                max[d] = maxOld;
                return;
            }
        }
    }


    public String getDescription() {
        return description;
    }


    public void setDescription(String description) {
        this.description = description;
    }


    public final int linearIndex(int[] x) {
        int index = x[nDimensions - 1];

        for (int d = nDimensions - 2; d >= 0; d--) {
            index = index * nBins[d] + x[d];
        }

        return index;
    }


    public final boolean containBin(int[] x) {
        for (int d = 0; d < nDimensions; d++) {
            if (x[d] < 0 || x[d] >= nBins[d]) {
                return false;
            }
        }

        return true;
    }


    public final boolean containCoord(double[] x) {
        double c;

        for (int d = 0; d < nDimensions; d++) {
            c = (int) Math.floor( (x[d] - min[d]) / binInterval[d]);
            if (c < 0 || c >= nBins[d]) {
                return false;
            }
        }

        return true;
    }


    public final double getBinValue(int[] x) {
        return bins[linearIndex(x)];
    }


    public final double getBinValue(int linearIndex) {
        return bins[linearIndex];
    }


    public final double getBinValueForCoords(double[] x) {
        return bins[linearIndex(getBinForCoordinate(x, tempBin))];
    }


    public final int[] coords(int linearIndex, int[] x) {
        for (int d = 0; d < nDimensions; d++) {
            x[d] = linearIndex % nBins[d];
            linearIndex /= nBins[d];
        }

        return x;
    }


    public final double[] middlePoint(int[] c, double[] x) {
        for (int d = 0; d < nDimensions; d++) {
            x[d] = min[d] + (c[d] + 0.5) * binInterval[d];
        }

        return x;
    }


    public final double getVolumeElement() {
        double vElement = 1;

        for (int d = 0; d < nDimensions; d++) {
            vElement *= binInterval[d];
        }

        return vElement;
    }


    public double getMaxValueBin(int[] x) {
        double maxValue = bins[0];
        int linearIndex = 0;

        for (int i = 1; i < bins.length; i++) {
            if (bins[i] > maxValue) {
                maxValue = bins[i];
                linearIndex = i;
            }
        }

        for (int d = 0; d < nDimensions; d++) {
            x[d] = linearIndex % nBins[d];
            linearIndex /= nBins[d];
        }

        return maxValue;
    }


    public final void fill(int[] x, int value) {
        bins[linearIndex(x)] += value;
    }


    public final void setBin(int[] x, short value) {
        bins[linearIndex(x)] = value;
    }


    public final void fill(double[] x, int value) {
        for (int d = 0; d < nDimensions; d++) {
            tempBin[d] = (int) Math.floor( (x[d] - min[d]) / binInterval[d]);
            if (tempBin[d] < 0) {
                underflow[d] += value;
                return;
            }
            if (tempBin[d] >= nBins[d]) {
                overflow[d] += value;
                return;
            }
        }

        bins[linearIndex(tempBin)] += value;
    }


    public final int[] getBinForCoordinate(double[] x, int[] bx) {
        for (int d = 0; d < nDimensions; d++) {
            bx[d] = (int) Math.floor( (x[d] - min[d]) / binInterval[d]);
        }

        return bx;
    }


    public int getBinCount() {
        return bins.length;
    }


    public double getBinSum() {
        double s = 0;

        for (int i = 0; i < bins.length; i++) {
            s += bins[i];
        }

        return s;
    }


    public double getBinInterval(int d) {
        return binInterval[d];
    }


    public void applyThreshold(double threshold) {
        for (int i = 0; i < bins.length; i++) {
            if (bins[i] <= threshold) {
                bins[i] = 0;
            }
        }
    }


    public void clear() {
        Arrays.fill(bins, (short) 0);
    }


    public int formDensityCluster(
        int[] x, double[] mean, double[] sigma, int n, double densityThreshold) {

        n++;
        final double previousDensity = getBinValue(x);
        setBin(x, (short) 0);

        int[] xi = new int[nDimensions];
        for (int d = 0; d < nDimensions; d++) {
            mean[d] += x[d];
            sigma[d] += x[d] * x[d];
            xi[d] = x[d] - 1;
        }

        double density;
        final int nCombinations = (int) Math.pow(3, nDimensions);

        for (int i = 0; i < nCombinations; i++) {
            if (containBin(xi) && !equal(xi, x)) {
                density = getBinValue(xi);
                if (density > densityThreshold && density <= previousDensity) {
                    n = formDensityCluster(xi, mean, sigma, n, densityThreshold);
                }
            }

            // n-dimensional increment loop
            for (int d = 0; d < nDimensions; d++) {
                xi[d]++;
                if (xi[d] == x[d] + 2) {
                    xi[d] = x[d] - 1;
                } else {
                    break;
                }
            }
        }

        return n;
    }


    private final boolean equal(int[] x1, int[] x2) {
        for (int d = 0; d < x1.length; d++) {
            if (x1[d] != x2[d]) {
                return false;
            }
        }

        return true;
    }

}
