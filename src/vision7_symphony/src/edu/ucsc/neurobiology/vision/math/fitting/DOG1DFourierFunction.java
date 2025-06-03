package edu.ucsc.neurobiology.vision.math.fitting;

import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.plot.*;


/**
 * Implements a fittable 1D Gaussian reprezented in frequency space.
 * Parameters:  "A1", "sigma1", "B" (pedestal), "A2", "sigma2".
 * The last 2 params are optional.
 *
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class DOG1DFourierFunction
    extends FittableFunction implements FunctionData {

    boolean surround = false;


    public DOG1DFourierFunction(double A1, double A2, double sigma1, double sigma2,
                                double B) {
        setParameters(new double[] {A1, sigma1, B, A2, sigma2});
        surround = true;
    }


    public DOG1DFourierFunction(double A1, double sigma1, double B) {
        setParameters(new double[] {A1, sigma1, B});
        surround = false;
    }


    public String getDescription() {
        return "";
    }


    public double getValueAndDerivatives(
        final double[] coord, final double[] parameters, final double[] derivatives) {

        double w = coord[0];
        final double a1 = parameters[0];
        final double s1 = parameters[1];
        final double b = parameters[2];
        final double exp1 = Math.exp( -w * w * s1 * s1 / 2);
        derivatives[0] = exp1;
        derivatives[1] = -s1 * a1 * w * w * exp1;
        derivatives[2] = 1;

        if (surround) {
            final double a2 = parameters[3];
            final double s2 = parameters[4];
            final double exp2 = Math.exp( -w * w * s2 * s2 / 2);

            derivatives[3] = exp2;
            derivatives[4] = -s2 * a2 * w * w * exp2;

            return a1 * exp1 + a2 * exp2 + b;
        } else {
            return a1 * exp1 + b;
        }
    }


    public double getValueAt(final double w) {
        final double a1 = parameters[0];
        final double s1 = parameters[1];
        final double b = parameters[2];
        final double exp1 = Math.exp( -w * w * s1 * s1 / 2);

        if (surround) {
            final double a2 = parameters[3];
            final double s2 = parameters[4];
            final double exp2 = Math.exp( -w * w * s2 * s2 / 2);

            return a1 * exp1 + a2 * exp2 + b;
        } else {
            return a1 * exp1 + b;
        }
    }


    public Num getValueAt(Num w) {
        Num w2 = w.sqr();

        //        final double exp1 = Math.exp( -w * w * s1 * s1 / 2);
//        final double exp2 = Math.exp( -w * w * s2 * s2 / 2);

        Num exp1 = new Num( -0.5, 0).mul(w2).mul(s1().sqr()).exp().mul(a1());
        Num exp2 = new Num( -0.5, 0).mul(w2).mul(s2().sqr()).exp().mul(a2());

        return exp1.add(exp2).add(b());
    }


    public String[] getParameterNames() {
        if (surround) {
            return new String[] {
                "A1", "sigma1", "B", "A2", "sigma2"};
        } else {
            return new String[] {
                "A1", "sigma1", "B"};
        }
    }


    public double getA1() {
        return parameters[0];
    }


    public double getS1() {
        return Math.abs(parameters[1]);
    }


    public double getB() {
        return parameters[2];
    }


    public double getA2() {
        return (surround) ? parameters[3] : -1;
    }


    public double getS2() {
        return (surround) ? parameters[4] : -1;
    }


    public Num a1() {
        return param(0);
    }


    public Num s1() {
        return param(1);
    }


    public Num b() {
        return param(2);
    }


    public Num a2() {
        return (surround) ? param(3) : new Num( -1, -1);
    }


    public Num s2() {
        return (surround) ? param(4) : new Num( -1, -1);
    }
}
