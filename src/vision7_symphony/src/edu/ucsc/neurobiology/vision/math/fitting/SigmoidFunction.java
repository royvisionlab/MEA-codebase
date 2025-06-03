package edu.ucsc.neurobiology.vision.math.fitting;

import edu.ucsc.neurobiology.vision.plot.*;


/**
 * Implements a fittable sigmoidal function.
 * Parameters: "A"-amplitude, "x0"-location, "sigma"-slope.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class SigmoidFunction
    extends FittableFunction implements FunctionData, Cloneable, PlotData {

    private String[] paramNames = {"A", "x0", "sigma"};


    public SigmoidFunction(double A, double x0, double sigma) {
        setParameters(new double[] {A, x0, sigma});
    }


    public String getDescription() {
        return "SpikeRateFunction";
    }


    public double getValueAndDerivatives(
        final double[] x, final double[] parameters, final double[] derivatives) {

        final double A = parameters[0];
        final double dx = x[0] - parameters[1];
        final double sigma = parameters[2];

        final double exp = Math.exp( -dx / sigma);
        final double f = 1 / (1 + exp);

        // derivative w.r.t A
        derivatives[0] = f;

        // derivative w.r.t x0
        derivatives[1] = -A * exp * f * f / sigma;

        // derivatives w.r.t sigma
        derivatives[2] = -A * dx * exp * f * f / (sigma * sigma);

        return A * f;
    }


    public double getValueAt(double x) {
        final double A = parameters[0];
        final double dx = x - parameters[1];
        final double sigma = parameters[2];

        return A / (1 + Math.exp( -dx / sigma));
    }


    public String[] getParameterNames() {
        return paramNames;
    }


    public double getAmplitude() {
        return parameters[0];
    }


    public double getX0() {
        return parameters[1];
    }


    public double getSigma() {
        return parameters[2];
    }


    /*
        public Tag read(int tagId, TaggedInputStream input, int len) throws IOException {
            return super.read(tagId, input, -1);
        }


        public void write(int tagId, TaggedOutputStream output) throws IOException {
            super.write(tagId, output);
        }
     */
}
