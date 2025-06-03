package edu.ucsc.neurobiology.vision.math.fitting;

import edu.ucsc.neurobiology.vision.plot.*;


/**
 * Implements a fittable 1D Difference of Gaussians (DOG) function.
 * Parameters:  "A1", "A2", "x0", "sigma1", "sigma2"
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class DOG1DFunction
    extends FittableFunction implements FunctionData {


    public DOG1DFunction(double A1, double A2, double x0, double sigma1, double sigma2) {
        setParameters(new double[] {A1, A2, x0, sigma1, sigma2});
    }


    public String getDescription() {
        return "";
    }


    public double getValueAndDerivatives(
        final double[] coord, final double[] parameters, final double[] derivatives) {

        double A1 = parameters[0];
        double A2 = parameters[1];
        double x0 = parameters[2];
        double s1 = parameters[3];
        double s2 = parameters[4];
        final double dx = coord[0] - x0;
        final double exp1 = StrictMath.exp( -0.5 * (dx * dx / (s1 * s1)));
        final double exp2 = StrictMath.exp( -0.5 * (dx * dx / (s2 * s2)));

        derivatives[0] = exp1;
        derivatives[1] = exp2;
        derivatives[2] = A1 * exp1 * dx / (s1 * s1) + A2 * exp2 * dx / (s2 * s2);
        derivatives[3] = A1 * exp1 * dx * dx / (s1 * s1 * s1);
        derivatives[4] = A2 * exp2 * dx * dx / (s2 * s2 * s2);

        return A1 * exp1 + A2 * exp2;
    }


    public double getValueAt(final double x) {
        double A1 = parameters[0];
        double A2 = parameters[1];
        double x0 = parameters[2];
        double s1 = parameters[3];
        double s2 = parameters[4];
        final double dx = x - x0;
        final double exp1 = StrictMath.exp( -0.5 * (dx * dx / (s1 * s1)));
        final double exp2 = StrictMath.exp( -0.5 * (dx * dx / (s2 * s2)));

        return A1 * exp1 + A2 * exp2;
    }


    public String[] getParameterNames() {
        return new String[] {
            "A1", "A2", "x0", "sigma1", "sigma2"};
    }


    public double getA1() {
        return parameters[0];
    }


    public double getA2() {
        return parameters[1];
    }


    public double getX0() {
        return parameters[2];
    }


    public double getSigma1() {
        return parameters[3];
    }


    public double getSigma2() {
        return parameters[4];
    }


}
