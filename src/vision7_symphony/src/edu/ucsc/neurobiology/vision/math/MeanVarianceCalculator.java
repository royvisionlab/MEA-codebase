package edu.ucsc.neurobiology.vision.math;


/**
 * A class that calculates the mean and variance of a set of data points added by add().
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class MeanVarianceCalculator {
    public static final int BIASED = 0;
    public static final int UNBIASED = 0;

    private double sumX;
    private double sumXSquare;
    private int n;
    private int type;


    public MeanVarianceCalculator(int type) {
        this.type = type;
        reset();
    }


    public MeanVarianceCalculator() {
        this.type = UNBIASED;
        reset();
    }


    public void reset() {
        sumX = 0;
        sumXSquare = 0;
        n = 0;
    }


    public void add(double x) {
        sumX += x;
        sumXSquare += x * x;
        n++;
    }


    public double getMean() {
//        if (n > 0) {
        return sumX / n;
//        } else {
//            return 0;
//        }
    }


    public double getStandardDeviation() {
//        if (n < 2) {
//            return 0;
//        } else {
        return Math.sqrt( (sumXSquare - sumX * sumX / n) / (n - 1));
//        }
    }


    public double getMeanVariance() {
        if (n > 0) {
            return getStandardDeviation() / Math.sqrt(n);
        } else {
            return 0;
        }
    }


    public String toString() {
        return getMean() + " \u00B1 " + getStandardDeviation();
    }

    /*
        public String toString(int fractionarDigits) {
            format.setMaximumFractionDigits(fractionarDigits);
            format.setMinimumFractionDigits(fractionarDigits);
            format.setMinimumIntegerDigits(1);
            format.setMaximumIntegerDigits(Integer.MAX_VALUE);

            return format.format(x) + "\u00B1" + format.format(Math.sqrt(s2));
        }
     */

}
