package edu.ucsc.neurobiology.vision.math;

import java.io.*;
import java.text.*;
import java.util.*;


/**
 * A value-error number that performs all operations but also does error
 * propagation on the result.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class Num
    implements Serializable {

    private static DecimalFormat format = new DecimalFormat();

    public static final Num minusOne = new Num( -1, 0);
    public static final Num zero = new Num(0, 0);
    public static final Num one = new Num(1, 0);
    public static final Num two = new Num(2, 0);

    public final double x;
    private final double s2;


    public Num(double x, double s) {
        this.x = x;
        this.s2 = s * s;
    }


    private Num(double x, double s2, int voidp) {
        this.x = x;
        this.s2 = s2;
    }


    public Num(Num n) {
        this.x = n.x;
        this.s2 = n.s2;
    }


    public double err() {
        return Math.sqrt(s2);
    }


    public Num inv() {
        return new Num( -x, s2, 0);
    }


    public Num add(Num n) {
        return new Num(x + n.x, s2 + n.s2, 0);
    }


    public Num add(double nx, double ns) {
        return new Num(x + nx, s2 + ns * ns, 0);
    }


    public Num sub(Num n) {
        return new Num(x - n.x, s2 + n.s2, 0);
    }


    public Num mul(Num n) {
        return new Num(x * n.x, n.x * n.x * s2 + x * x * n.s2, 0);
    }


    public Num mul(double nx, double ns) {
        return new Num(x * nx, nx * nx * s2 + x * x * ns * ns, 0);
    }


    public Num mul(double nx) {
        return mul(nx, 0);
    }


    public Num div(Num n) {
        return new Num(x / n.x, (n.x * n.x * s2 + x * x * n.s2) / (n.x * n.x * n.x * n.x),
                       0);
    }


    public Num div(double nx, double ns) {
        return new Num(x / nx, (nx * nx * s2 + x * x * ns * ns) / (nx * nx * nx * nx), 0);
    }


    public Num div(double nx) {
        return div(nx, 0);
    }
    
    public Num sub(double nx, double ns) {
        return new Num(x - nx,s2 + ns*ns);
    }
    
    public Num sub(double nx) {
        return sub(nx, 0);
    }


    public Num exp() {
        return new Num(Math.exp(x), Math.exp(2 * x) * s2, 0);
    }


    public Num pow(Num y) {
        double err2 = Math.pow(x, 2 * y.x) * (y.x * y.x * s2 / (x * x) +
                                              Math.log(x) * Math.log(x) * y.s2);
        return new Num(Math.pow(x, y.x), err2, 0);
    }


    public Num ln() {
        return new Num(Math.log(x), s2 / (x * x), 0);
    }


    public Num lg() {
        return new Num(Math.log10(x), s2 / Math.pow(x * Math.log(10), 2), 0);
    }


    public Num sqrt() {
        return new Num(Math.sqrt(x), s2 / (4 * x), 0);
    }


    public Num sqr() {
        return new Num(x * x, 4 * x * x * s2, 0);
    }


    public String toString() {
        return x + " \u00B1 " + Math.sqrt(s2);
    }


    //Written by Matthew Grivich
    //Returns Num in significant digits format.
    //Not fully bug tested.  It has not been implemented for <0 fractional digits
    public String toString(boolean dummy) {
        int fractionalDigits = 0;
        if (Math.sqrt(s2) == 0) {
            return x + " \u00B1 " + Math.sqrt(s2);
        }
        if (Math.sqrt(s2) >= 1.0) {
            fractionalDigits = 0;
        } else {
            double temp = Math.sqrt(s2);
            while (temp < 1) {
                temp *= 10;
                fractionalDigits++;
            }
        }
        DecimalFormat format = new DecimalFormat();
        format.setMaximumFractionDigits(fractionalDigits);
        format.setMinimumFractionDigits(fractionalDigits);
        return format.format(x) + " \u00B1 " +
            format.format(Math.sqrt(s2));
//        return x + " \u00B1 " + Math.sqrt(s2);
    }


    public String toString(int fractionarDigits) {
        format.setMaximumFractionDigits(fractionarDigits);
        format.setMinimumFractionDigits(fractionarDigits);
        format.setMinimumIntegerDigits(1);
        format.setMaximumIntegerDigits(Integer.MAX_VALUE);

        return format.format(x) + "\u00B1" + format.format(Math.sqrt(s2));
    }


    public static Num average(ArrayList<Num> values) {
        Num[] valuesArray = new Num[values.size()];
        for (int i = 0; i < values.size(); i++) {
            valuesArray[i] = values.get(i);
        }
        return average(valuesArray);
    }


    public static Num averageRescaleErrors(ArrayList<Num> values) {
        Num[] valuesArray = new Num[values.size()];
        for (int i = 0; i < values.size(); i++) {
            valuesArray[i] = values.get(i);
        }
        return averageRescaleErrors(valuesArray);
    }


    public static Num average(Num ...values) {
        double x = 0, w = 0;
        for (Num v : values) {
            if (v.s2 > 0) {
                x += v.x / v.s2;
                w += 1 / v.s2;
            }
        }
        return new Num(x / w, 1 / w, 0);
    }


    //Review of Particle Physics (long version, 2002) Introduction.  Section 4.2.2
    //Written by Matthew Grivich
    public static Num averageRescaleErrors(Num ...values) {
        if (values.length > 1) {
            int n = 0;
            double x = 0, w = 0;
            //Calculate weighted average
            for (Num v : values) {
                //Do not include points with no error
                if (v.s2 > 0) {
                    x += v.x / v.s2;
                    w += 1 / v.s2;
                    n++;
                }
            }
            double newX = x / w;
            double newS2 = 1 / w;

            //Calculate chi^2 value for this average being reasonable.
            double chi2 = 0.0;
            for (Num v : values) {
                if (v.s2 > 0) {
                    chi2 += (1 / v.s2) * (newX - v.x) * (newX - v.x);
                }
            }

            //If unreasonable, increase the error bars on the initial measurements
            //until it becomes reasonable.
            if (chi2 > n - 1) {
                newS2 = newS2 * chi2 / (n - 1);
            }
            return new Num(newX, newS2, 0);
        } else if (values.length == 1) {
            return values[0];
        } else {
            return Num.zero;
        }
    }


    public Num average(Num v) {
        double w = 1 / s2 + 1 / v.s2;
        return new Num( (x / s2 + v.x / v.s2) / w, 1 / w, 0);
    }


    public static Num rms(Num ...values) {
        MeanVarianceCalculator mvc = new MeanVarianceCalculator();
        for (Num v : values) {
            mvc.add(v.x);
        }
        return new Num(mvc.getMean(), mvc.getMeanVariance());
    }


    public static void main(String[] args) {
        System.err.println(Num.averageRescaleErrors(
            new Num(1, 0.001), new Num(1, 0.001), new Num(1, 0.001)));
    }

}
