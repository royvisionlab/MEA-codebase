package edu.ucsc.neurobiology.vision.math;

import java.awt.geom.Point2D;
import java.util.ArrayList;


/**
 * This class is intended to contain miscelaneous math-related utility
 * methods used allover the program.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public final class MathUtil {
    public final static double LN2 = Math.log(2);
    public final static double LN10 = Math.log(10);


    public static boolean areEqual(double x1, double x2) {
        return Math.abs(x1 - x2) < 1e-6;
    }


    public static final double project(double[] x, double[] y) {
        double v = 0;

        for (int i = 0; i < x.length; i++) {
            v += x[i] * y[i];
        }
        return v;
    }


    public static double fastAbs(double a) {
        return (a >= 0) ? a : -a;
    }


    /**
     * Adds the value <i>y</i> to each elecment of array <i>x</i>.
     */
    public static final void add(double[] x, double y) {
        for (int i = 0; i < x.length; i++) {
            x[i] += y;
        }
    }


    public static final void sub(double[] x, double y) {
        for (int i = 0; i < x.length; i++) {
            x[i] -= y;
        }
    }


    /**
     * Adds two arrays element by element. The result is stored in <i>x</i>,
     * <i>y</i> remains unchanged.
     */
    public static final void add(double[] x, double[] y) {
        if (x == null) {
            throw new NullPointerException("x is null");
        }
        if (y == null) {
            throw new NullPointerException("y is null");
        }

        if (x.length != y.length) {
            throw new IllegalArgumentException("Impossible to add, different size");
        }

        for (int i = 0; i < x.length; i++) {
            x[i] += y[i];
        }
    }


    /**
     * Multiplies each elecment of array <i>x</i> with the value <i>y</i>.
     */
    public static double[] multiply(double[] x, double y) {
        for (int i = 0; i < x.length; i++) {
            x[i] *= y;
        }
        return x;
    }


    /**
     * Multiplies each elecment of array <i>x</i> with the value <i>y</i>.
     */
    public static final float[] multiply(float[] x, float y) {
        for (int i = 0; i < x.length; i++) {
            x[i] *= y;
        }
        return x;
    }


    /**
     * Divides each elecment of array <i>x</i> by the value <i>y</i>.
     */
    public static final double[] divide(double[] x, double y) {
        for (int i = 0; i < x.length; i++) {
            x[i] /= y;
        }
        return x;
    }


    /**
     * Divides each elecment of array <i>x</i> by the value <i>y</i>.
     */
    public static final float[] divide(float[] x, double y) {
        for (int i = 0; i < x.length; i++) {
            x[i] /= y;
        }
        return x;
    }


    /**
     * Finds the minimum of the values of the array.
     * @return the minimum value
     */
    public static final double min(double[] x) {
        double min = Double.POSITIVE_INFINITY;

        for (int i = 0; i < x.length; i++) {
            if (x[i] < min) {
                min = x[i];
            }
        }

        return min;
    }


    public static final double extreme(double[] x) {
        double min = x[0];
        double max = x[0];

        for (int i = 1; i < x.length; i++) {
            if (x[i] < min) {
                min = x[i];
            }
            if (x[i] > max) {
                max = x[i];
            }
        }

        return (Math.abs(min) > Math.abs(max)) ? min : max;
    }


    public static final double rms(float[] x) {
        double mean = mean(x);

        double rms = 0;
        for (int i = 0; i < x.length; i++) {
            rms += (x[i] - mean) * (x[i] - mean);
        }
        rms /= x.length;

        return Math.sqrt(rms);
    }


    public static final float min(float[] x) {
        float min = Float.POSITIVE_INFINITY;

        for (int i = 0; i < x.length; i++) {
            if (x[i] < min) {
                min = x[i];
            }
        }

        return min;
    }


    /**
     * Finds the minimum of the values of the array.
     * @return the minimum value
     */
    public static final int min(int[] x) {
        int min = Integer.MAX_VALUE;

        for (int i = 0; i < x.length; i++) {
            if (x[i] < min) {
                min = x[i];
            }
        }

        return min;
    }


    /**
     * Finds the maximum of the values of the array.
     * @return the maximum value
     */
    public static final double max(double[] x) {
        double max = x[0];

        for (int i = 1; i < x.length; i++) {
            if (x[i] > max) {
                max = x[i];
            }
        }

        return max;
    }


    /**
     * Finds the maximum of the values of the array.
     * @return the maximum value
     */
    public static final double max(double[] x, int i1, int i2) {
        double max = Double.NEGATIVE_INFINITY;

        for (int i = i1; i <= i2; i++) {
            if (x[i] > max) {
                max = x[i];
            }
        }

        return max;
    }


    public static final float max(float[] x) {
        float max = x[0];

        for (int i = 1; i < x.length; i++) {
            if (x[i] > max) {
                max = x[i];
            }
        }

        return max;
    }


    /**
     * Finds the maximum of the values of the array.
     * @return the maximum value
     */
    public static final int max(int[] x) {
        int max = x[0];

        for (int i = 1; i < x.length; i++) {
            if (x[i] > max) {
                max = x[i];
            }
        }

        return max;
    }


    /**
     * Finds the maximum absolute value (ignoring the sign) of the values of the array.
     * @return the maximum absolute value
     */
    public static final double maxAbs(double[] x) {
        double max = x[0];

        for (int i = 1; i < x.length; i++) {
            if (Math.abs(x[i]) > max) {
                max = Math.abs(x[i]);
            }
        }

        return max;
    }


    /**
     * Finds the maximum absolute value (ignoring the sign) of the values of the array.
     * @return the maximum absolute value
     */
    public static final float maxAbs(float[] x) {
        float max = x[0];

        for (int i = 1; i < x.length; i++) {
            if (Math.abs(x[i]) > max) {
                max = Math.abs(x[i]);
            }
        }

        return max;
    }


    /**
     * Fixed IndexOutOfBounds bug. maxIndex will never be -1.
     *
     * @param x double[]
     * @return int
     */
    public static final int maxAbsIndex(double[] x) {
        int maxIndex = 0;

        double maxAbs = x[0];
        for (int i = 1; i < x.length; i++) {
            if (Math.abs(x[i]) > maxAbs) {
                maxIndex = i;
                maxAbs = Math.abs(x[maxIndex]);
            }
        }

        return maxIndex;
    }


    public static final int closestMatch(final double v, double[] x) {
        int index = -1;

        double minDifference = Double.MAX_VALUE;
        for (int i = 0; i < x.length; i++) {
            if (Math.abs(v - x[i]) < minDifference) {
                minDifference = Math.abs(v - x[i]);
                index = i;
            }
        }

        return index;
    }


    public static final int maxAbsIndex(float[] x) {
        int maxIndex = -1;

        float maxAbs = Float.NEGATIVE_INFINITY;
        for (int i = 0; i < x.length; i++) {
            if (Math.abs(x[i]) > maxAbs) {
                maxIndex = i;
                maxAbs = Math.abs(x[maxIndex]);
            }
        }

        return maxIndex;
    }


    /**
     * Finds the maximum absolute value (ignoring the sign) of the values of the array.
     * @return the maximum absolute value
     */
    public static final int maxAbs(int[] x) {
        int max = x[0];

        for (int i = 1; i < x.length; i++) {
            if (Math.abs(x[i]) > max) {
                max = Math.abs(x[i]);
            }
        }

        return max;
    }


    /**
     * Finds the index of the maximum element in the array.
     * @param x the array
     * @return the the index of the biggest element
     */
    public static final int maxIndex(double[] x) {
        double max = x[0];
        int index = 0;

        for (int i = 1; i < x.length; i++) {
            if (x[i] > max) {
                max = x[i];
                index = i;
            }
        }

        return index;
    }


    public static final int maxIndex(int[] x) {
        int max = x[0];
        int index = 0;

        for (int i = 1; i < x.length; i++) {
            if (x[i] > max) {
                max = x[i];
                index = i;
            }
        }

        return index;
    }


    /**
     * Finds the index of the maximum element in the array.
     * @param x the array
     * @return the the index of the biggest element
     */
    public static final int maxIndex(float[] x) {
        float max = x[0];
        int index = 0;

        for (int i = 1; i < x.length; i++) {
            if (x[i] > max) {
                max = x[i];
                index = i;
            }
        }

        return index;
    }


    public static final int maxIndex(float[] x, int i1, int i2) {
        float max = x[i1];
        int index = 0;

        for (int i = i1 + 1; i < i2; i++) {
            if (x[i] > max) {
                max = x[i];
                index = i;
            }
        }

        return index;
    }


    /**
     * Finds the index of the minimum element in the array.
     * @param x the array
     * @return the the index of the smallest element
     */
    public static final int minIndex(double[] x) {
        double min = x[0];
        int index = 0;

        for (int i = 1; i < x.length; i++) {
            if (x[i] < min) {
                min = x[i];
                index = i;
            }
        }

        return index;
    }


    public static final int minIndex(float[] x) {
        double min = x[0];
        int index = 0;

        for (int i = 1; i < x.length; i++) {
            if (x[i] < min) {
                min = x[i];
                index = i;
            }
        }

        return index;
    }


    public static final int minIndex(int[] x) {
        int min = x[0];
        int index = 0;

        for (int i = 1; i < x.length; i++) {
            if (x[i] < min) {
                min = x[i];
                index = i;
            }
        }

        return index;
    }


    /**
     * Calculates the sum of the elements of the array.
     * @param x the array
     * @return the sum of elements
     */
    public static final double sum(double[] x) {
        double sum = 0;

        for (int i = 0; i < x.length; i++) {
            sum += x[i];
        }

        return sum;
    }


    public static final long sum(long[] x) {
        long sum = 0;

        for (int i = 0; i < x.length; i++) {
            sum += x[i];
        }

        return sum;
    }


    /**
     * Calculates the sum of the elements of the array.
     * @param x the array
     * @return the sum of elements
     */
    public static final float sum(float[] x) {
        float sum = 0;

        for (int i = 0; i < x.length; i++) {
            sum += x[i];
        }

        return sum;
    }


    public static final int sum(int[] x) {
        int sum = 0;

        for (int i = 0; i < x.length; i++) {
            sum += x[i];
        }

        return sum;
    }


    public static final double sum(double[] x, int i1, int i2) {
        double sum = 0;

        for (int i = i1; i < i2; i++) {
            sum += x[i];
        }

        return sum;
    }


    /**
     * Calculates the sum of the elements of the array.
     * @param x the array
     * @return the sum of elements
     */
    public static final double sumAbs(double[] x) {
        double sum = 0;

        for (int i = 0; i < x.length; i++) {
            sum += Math.abs(x[i]);
        }

        return sum;
    }


    /**
     * Calculates the sum of the elements of the array.
     * @param x the array
     * @return the sum of elements
     */
    public static final float sumAbs(float[] x) {
        float sum = 0;

        for (int i = 0; i < x.length; i++) {
            sum += Math.abs(x[i]);
        }

        return sum;
    }


    public static final double sumSquare(double[] x) {
        double sum = 0;

        for (int i = 0; i < x.length; i++) {
            sum += x[i] * x[i];
        }

        return sum;
    }


    public static final float sumSquare(float[] x) {
        float sum = 0;

        for (int i = 0; i < x.length; i++) {
            sum += x[i] * x[i];
        }

        return sum;
    }


    /**
     *
     * @param array double[] array to be resized
     * @param newSize int length of new array
     * @return double[] array with new size
     */
    public static final double[] resize(double[] array, int newSize) {
        double[] newArray = new double[newSize];
        for (int i = 0; i < newSize && i < array.length; i++) {
            newArray[i] = array[i];
        }
        return newArray;
    }


    /**
     * Does a simple linear (left to right) search for the value <i>v</i>.
     * @param x the array to be searched
     * @return the index of the first value <i>v</i>
     */
    public static final int search(int[] x, int v) {
        for (int i = 0; i < x.length; i++) {
            if (x[i] == v) {
                return i;
            }
        }

        return -1;
    }


    /**
     * Calculates the factorial of the value k
     *
     * @param k int
     * @return double
     */
    public static final double fact(int k) {
        if ( (k == 0) || (k == 1)) {
            return 1;
        } else {
            return k * fact(k - 1);
        }
    }


    public static final double normalizeDegrees(double deg, final double range) {
        while (deg >= range) {
            deg -= range;
        } while (deg < 0) {
            deg += range;

        }
        return deg;
    }


    public static final boolean areEqual(int[] x1, int[] x2) {
        if (x1.length != x2.length) {
            return false;
        } else {
            ArrayList<Integer> l1 = new ArrayList<Integer>();
            for (int i = 0; i < x1.length; i++) {
                l1.add(new Integer(x1[i]));
            }
            for (int i = 0; i < x1.length; i++) {
                Integer n = new Integer(x2[i]);
                if (l1.contains(n)) {
                    l1.remove(n);
                } else {
                    return false;
                }
            }
            return true;
        }
    }


    public static final double log2(double x) {
        return Math.log(x) / LN2;
    }


    public static final double log10(double x) {
        return Math.log(x) / LN10;
    }


    public static final double SQR(double x) {
        return x * x;
    }


    public static final double POW6(double x) {
        return x * x * x * x * x * x;
    }


    /**
     * Calculates the mean (0) ans standard variation (1) of an array of values.
     *
     * @param x the array with values
     * @return {mean, sigma}
     */
    public static final double[] calculateStatistics(double[] x) {
        final double mean = mean(x);

        double sigma = 0.0;
        for (int i = 0; i < x.length; i++) {
            sigma += (x[i] - mean) * (x[i] - mean);
        }

        return new double[] {
            mean, Math.sqrt(sigma / (x.length - 1))};
    }


    /**
     * Counts how many times the value <i>v</i> occurs in the array <i>array</i>.
     */
    public static final int countValues(int v, int[] array) {
        int n = 0;

        for (int i = 0; i < array.length; i++) {
            if (array[i] == v) {
                n++;
            }
        }

        return n;
    }


    /**
     * Counts how many times the value <i>v</i> occurs in the array <i>array</i>.
     */
    public static final int countValues(boolean v, boolean[] array) {
        int n = 0;

        for (int i = 0; i < array.length; i++) {
            if (array[i] == v) {
                n++;
            }
        }

        return n;
    }


    /**
     * Counts how many times the value <i>v</i> occurs in the array <i>array</i>.
     */
    public static final int countValues(boolean v, boolean[][] array) {
        int n = 0;

        for (int i = 0; i < array.length; i++) {
            for (int j = 0; j < array[0].length; j++) {
                if (array[i][j] == v) {
                    n++;
                }
            }
        }

        return n;
    }


    public static final double sqr(double x) {
        return x * x;
    }


    public static final int[] intArray(int n1, int n2) {
        int[] x = new int[n2 - n1];
        for (int i = 0; i < x.length; i++) {
            x[i] = n1 + i;
        }
        return x;
    }


    public static final double[] doubleArray(int n1, int n2) {
        double[] x = new double[n2 - n1];
        for (int i = 0; i < x.length; i++) {
            x[i] = n1 + i;
        }
        return x;
    }


    public static final boolean isPrime(int n) {
        for (int i = 2; i <= n - 1; i++) {
            if (n % i == 0) {
                return false;
            }
        }
        return true;
    }


    public static final double getPointToLineDistance(
        double x, double y, double x1, double y1, double x2, double y2, boolean segment) {

        double u = ( (x - x1) * (x2 - x1) + (y - y1) * (y2 - y1))
                   / ( (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));

        if (segment) {
            if (u < 0) {
                return Math.sqrt( (x - x1) * (x - x1) + (y - y1) * (y - y1));
            } else if (u > 1) {
                return Math.sqrt( (x - x2) * (x - x2) + (y - y2) * (y - y1));
            }
        }

        double xx = x1 + u * (x2 - x1);
        double yy = y1 + u * (y2 - y1);
        return Math.sqrt( (x - xx) * (x - xx) + (y - yy) * (y - yy));
    }


    public static final void getLineIntersection(
        double x1, double y1, double x2, double y2,
        double x3, double y3, double x4, double y4,
        Point2D.Double point, boolean segment) {

        double v = (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1);
        if (v == 0) {
            point.x = Double.NaN;
            point.y = Double.NaN;
            return;
        }
        double ua = ( (x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / v;
        double ub = ( (x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / v;

        if (segment && (ua < 0 || ua > 1 || ub < 0 || ub > 1)) {
            point.x = Double.NaN;
            point.y = Double.NaN;
            return;
        }

        point.x = x1 + ua * (x2 - x1);
        point.y = y1 + ua * (y2 - y1);
    }


    /**
     * Calculates the mean of an array of values.
     *
     * @param x the array with values
     * @return the mean
     */
    public static final double mean(double[] x) {
        double mean = 0;

        for (int i = 0; i < x.length; i++) {
            mean += x[i];
        }

        return mean / x.length;
    }


    public static final double mean(float[] x) {
        double mean = 0;

        for (int i = 0; i < x.length; i++) {
            mean += x[i];
        }

        return mean / x.length;
    }


    public static final double mean(double[] x, int i1, int i2) {
        if (i2 <= i1) {
            throw new IllegalArgumentException("(i2 <= i1)");
        }
        double mean = 0;

        for (int i = i1; i < i2; i++) {
            mean += x[i];
        }

        return mean / (i2 - i1);
    }


    /**
     * Calculates the mean of an array of values.
     *
     * @param x the array with values
     * @return the mean
     */
    public static final double mean(int[] x) {
        double mean = 0;

        for (int i = 0; i < x.length; i++) {
            mean += x[i];
        }

        return mean / x.length;
    }


    public static void rotateArray(float[] a, int k) {
        int c, v;
        float tmp;

        if (a == null || a.length == 0) {
            return;
        }
        if (k < 0 || k >= a.length) {
            k %= a.length;
            if (k < 0) {
                k += a.length;
            }
        }
        if (k == 0) {
            return;
        }

        c = 0;
        for (v = 0; c < a.length; v++) {
            int t = v, tp = v + k;
            tmp = a[v];
            c++;
            while (tp != v) {
                a[t] = a[tp];
                t = tp;
                tp += k;
                if (tp >= a.length) {
                    tp -= a.length;
                }
                c++;
            }
            a[t] = tmp;
        }
    }


    public static void rotateArray(double[] a, int k) {
        int c, v;
        double tmp;

        if (a == null || a.length == 0) {
            return;
        }
        if (k < 0 || k >= a.length) {
            k %= a.length;
            if (k < 0) {
                k += a.length;
            }
        }
        if (k == 0) {
            return;
        }

        c = 0;
        for (v = 0; c < a.length; v++) {
            int t = v, tp = v + k;
            tmp = a[v];
            c++;
            while (tp != v) {
                a[t] = a[tp];
                t = tp;
                tp += k;
                if (tp >= a.length) {
                    tp -= a.length;
                }
                c++;
            }
            a[t] = tmp;
        }
    }

}
