package edu.ucsc.neurobiology.vision.math.fitting;


/**
 * Implements a fittable 2D generalized Difference of Gaussians (DOG) function.
 * Parameters:  "x0", "y0", "sigX", "sigY", "r" (surround to center size ratio), "A1", "A2", "theta"
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class DOG2DFunction
    extends FittableFunction {


    public DOG2DFunction() {
        setParameters(new double[8]);
    }


    public double getValueAndDerivatives(
        final double[] coord, final double[] parameters, final double[] derivatives) {

        final double dx = coord[0] - parameters[0];
        final double dy = coord[1] - parameters[1];
        final double s1 = parameters[2];
        final double s2 = parameters[3];
        final double r = parameters[4];
        final double A1 = parameters[5];
        final double A2 = parameters[6];
        final double sin = Math.sin(parameters[7]);
        final double cos = Math.cos(parameters[7]);

        final double s1Sq = s1 * s1;
        final double s2Sq = s2 * s2;

        double sum1 = dx * cos + dy * sin;
        double sum2 = dy * cos - dx * sin;

        final double G1 = A1 * StrictMath.exp(
            -sum1 * sum1 / (2 * s1 * s1) - sum2 * sum2 / (2 * s2 * s2));
        final double G2 = A2 * StrictMath.exp(
            -sum1 * sum1 / (2 * r * s1 * r * s1) - sum2 * sum2 / (2 * r * s2 * r * s2));

        double a = G1 + G2 / (r * r);

        derivatives[0] = (sin * sum2 / s2Sq + cos * sum1 / s1Sq) * a;
        derivatives[1] = (cos * sum2 / s2Sq - sin * sum1 / s1Sq) * a;

        derivatives[2] = (sum1 * sum1 / (s1 * s1 * s1)) * a;
        derivatives[3] = (sum2 * sum2 / (s2 * s2 * s2)) * a;

        derivatives[4] = G2 * (sum1 * sum1 / s1Sq + sum2 * sum2 / s2Sq) / (r * r * r);

        derivatives[5] = G1 / A1;
        derivatives[6] = G2 / A2;

        derivatives[7] = sum1 * sum2 * (1 / s1Sq + 1 / s2Sq) * a;

        return G1 + G2;
    }


    public String[] getParameterNames() {
        return new String[] {
            "x0", "y0", "sigX", "sigY", "r", "A1", "A2", "theta"};
    }


    public double getX0() {
        return parameters[0];
    }


    public double getY0() {
        return parameters[1];
    }


    public double getSigmaX1() {
        return parameters[2];
    }


    public double getSigmaY1() {
        return parameters[3];
    }


    public double getSigmaX2() {
        return parameters[2] * parameters[4];
    }


    public double getSigmaY2() {
        return parameters[3] * parameters[4];
    }


    public double getA1() {
        return parameters[5];
    }


    public double getA2() {
        return parameters[6];
    }


    public double getTheta() {
        return parameters[7];
    }
}
