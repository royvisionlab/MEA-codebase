package edu.ucsc.neurobiology.vision.math.fitting;

import edu.ucsc.neurobiology.vision.plot.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class DOG1DFourierFunction1
    extends FittableFunction implements FunctionData {


    public DOG1DFourierFunction1(double A1, double sigma1, double B) {
        setParameters(new double[] {A1, sigma1, B});
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

        return a1 * exp1 + b;
    }


    public double getValueAt(final double w) {
        final double a1 = parameters[0];
        final double s1 = parameters[1];
        final double b = parameters[2];

        final double exp1 = Math.exp( -w * w * s1 * s1 / 2);

        return a1 * exp1 + b;
    }


    public String[] getParameterNames() {
        return new String[] {"A1", "sigma1", "B"};
    }


    public double getA1() {
        return parameters[0];
    }


    public double getS1() {
        return parameters[1];
    }


    public double getB() {
        return parameters[2];
    }

}
