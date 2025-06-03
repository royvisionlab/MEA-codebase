package edu.ucsc.neurobiology.vision.plot;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class AxisUtil {
    private static final double oneOverLn10 = 1 / Math.log(10);


    public static final double ceil(double d, int exp) {
        double x = 1.0 * Math.pow(10.0, (double) exp);
        return Math.ceil(d / x) * x;
    }


    public static final double floor(double d, int exp) {
        double x = 1.0 * Math.pow(10.0, (double) exp);
        return Math.floor(d / x) * x;
    }
    
    public static final double round(double d, int exp) {
        double x = 1.0 * Math.pow(10.0, (double) exp);
        return Math.round(d / x) * x;    	
    }


    public static final double[] performAutoScale(double min, double max) {
        double[] d = new double[2];

        double diff = max - min;
        d[0] = floor(min, getPowerOf10(diff));
        d[1] = ceil(max, getPowerOf10(diff));

        return d;
    }


    public static final double calculateTickSpacing(double min, double max) {
//        return Math.pow(10, getRoundPowerOf10(max - min) - 1);

        return Math.pow(10, (int) Math.round(Math.log(max - min) / Math.log(10)) - 1);
    }


    private static final double log10(double x) {
        return Math.log(x) * oneOverLn10;
    }


    public static final int getPowerOf10(double x) {
        // should be floor
        return (int) Math.floor(log10(x));
    }


    private static final int getRoundPowerOf10(double x) {
        return (int) Math.round(log10(x));
    }


    public static final int getFractionalDigits(double x) {
        for (int i = 0, mult = 1; i < Integer.MAX_VALUE; i++, mult *= 10) {
            double doubleValue = x * mult;
            int intValue = (int) Math.floor(doubleValue);
            if (doubleValue == intValue) {
                return i;
            }
        }
        throw new IllegalArgumentException(
            "Impossible to get the number of fractional digits for: " + x);
    }

}
