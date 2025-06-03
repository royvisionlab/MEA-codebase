package edu.ucsc.neurobiology.vision.math.fitting;

import edu.ucsc.neurobiology.vision.plot.*;


/**
 * Implements a fittable 1D gaussian function.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class Gaussian1DFunction
    extends FittableFunction implements FunctionData {


    private final int type;


    public Gaussian1DFunction(double B, double A, double x0, double sig) {
        setParameters(new double[] {A, x0, sig, B});
        type = 0;
    }


    public Gaussian1DFunction(double A, double x0, double sig) {
        setParameters(new double[] {A, x0, sig});
        type = 1;
    }


    public double getValueAndDerivatives(
        final double[] x, final double[] parameters, final double[] derivatives) {

        if (type == 0) {
            double A = parameters[0];
            double dx = (x[0] - parameters[1]);
            double sig = parameters[2];
            double B = parameters[3];

            derivatives[0] = StrictMath.exp( -dx * dx / (2 * sig * sig));
            derivatives[1] = A * derivatives[0] * dx / (sig * sig);
            derivatives[2] = derivatives[1] * dx / sig;
            derivatives[3] = 1;

            return B + A * derivatives[0];
        } else {
            double A = parameters[0];
            double dx = (x[0] - parameters[1]);
            double sig = parameters[2];

            derivatives[0] = StrictMath.exp( -dx * dx / (2 * sig * sig));
            derivatives[1] = A * derivatives[0] * dx / (sig * sig);
            derivatives[2] = derivatives[1] * dx / sig;

            return A * derivatives[0];
        }
    }


    public double getValueAt(double x) {
        double A = parameters[0];
        double dx = (x - parameters[1]);
        double sig = parameters[2];
        double B;

        if (type == 0) {
            B = parameters[3];
        } else {
            B = 0;
        }

        return B + A * StrictMath.exp( -dx * dx / (2 * sig * sig));
    }


    public String getDescription() {
        return "";
    }


    public double getSigma() {
        return Math.abs(parameters[2]);
    }


    public double getAmplitude() {
        return parameters[0];
    }


    public double getBackground() {
        double B;
        if (type == 0) {
            B = parameters[3];
        } else {
            B = 0;
        }
        return B;
    }


    public double getMean() {
        return parameters[1];
    }
}
