package edu.ucsc.neurobiology.vision.math.fitting;

import edu.ucsc.neurobiology.vision.plot.*;


/**
 * Implements a fittable exponential function.
 * Parameters: N0-amplitude, tau-the decay constant
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class Exponential
    extends FittableFunction implements FunctionData {


    public Exponential(double N0, double tau) {
        setParameters(new double[] {N0, tau});
    }


    public String getDescription() {
        return "";
    }


    public double getValueAndDerivatives(
        final double[] coord, final double[] parameters, final double[] derivatives) {

        final double N0 = parameters[0];
        final double tau = parameters[1];
        derivatives[0] = Math.exp( -coord[0] / tau);
        //             derivatives[1] = ;


        return N0 * Math.exp( -coord[0] / tau);
    }


    public double getValueAt(final double x) {
        double A1 = parameters[0];
        double A2 = parameters[1];
        double x0 = parameters[2];
        double sigma1 = parameters[3];
        double sigma2 = parameters[4];

        return
            A1 * StrictMath.exp( -0.5 * (x - x0) * (x - x0) / (sigma1 * sigma1)) +
            A2 * StrictMath.exp( -0.5 * (x - x0) * (x - x0) / (sigma2 * sigma2));
    }


    public String[] getParameterNames() {
        return new String[] {
            "A1", "A2", "x0", "sigma1", "sigma2"};
    }

    /*
        public double getTheta() {
            return parameters[5];
        }
        public double getSigmaY() {
            return parameters[4];
        }
        public double getSigmaX() {
            return parameters[3];
        }
        public double getY0() {
            return parameters[2];
        }
        public double getX0() {
            return parameters[1];
        }
        public double getAmplitude() {
            return parameters[0];
        }
        public void setAmplitude(double a) {
            parameters[0] = a;
        }
        public double getPedestal() {
            return parameters[6];
        }
     */
}
