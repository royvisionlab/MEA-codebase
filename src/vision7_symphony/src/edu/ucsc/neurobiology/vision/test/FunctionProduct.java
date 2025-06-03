package edu.ucsc.neurobiology.vision.test;

import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.math.fitting.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class FunctionProduct
    extends FittableFunction {

    private FittableFunction f1;
    private FittableFunction f2;
    double[] p1, p2;
    double[] d1, d2;


    public final static double[] concatenate(double[] a1, double[] a2) {
        double[] a = new double[a1.length + a2.length];
        System.arraycopy(a1, 0, a, 0, a1.length);
        System.arraycopy(a2, 0, a, a1.length, a2.length);
        return a;
    }


    public final static double[] subArray(double[] a, int i1, int i2) {
        double[] sub = new double[i2 - i1 + 1];
        System.arraycopy(a, i1, sub, 0, sub.length);
        return sub;
    }


    public FunctionProduct(FittableFunction f1, FittableFunction f2) {
        this.f1 = f1;
        this.f2 = f2;

        this.p1 = new double[f1.getParametersCount()];
        this.p2 = new double[f2.getParametersCount()];
        this.d1 = new double[f1.getParametersCount()];
        this.d2 = new double[f2.getParametersCount()];

        setParameters(concatenate(f1.getParameters(), f2.getParameters()));
    }


    public double getValueAndDerivatives(
        final double[] x, final double[] parameters, final double[] derivatives) throws
        CannotEvaluateException {

        System.arraycopy(parameters, 0, p1, 0, p1.length);
        System.arraycopy(parameters, p1.length, p2, 0, p2.length);

        System.arraycopy(derivatives, 0, d1, 0, d1.length);
        System.arraycopy(derivatives, d1.length, d2, 0, d2.length);

        double v1 = f1.getValueAndDerivatives(x, p1, d1);
        double v2 = f2.getValueAndDerivatives(x, p2, d2);

        return v1 * v2;
    }


    /*
        public final double getValue(final double x, final double y) {
            double x0 = parameters[1];
            double y0 = parameters[2];
            double ct = Math.cos(parameters[5]);
            double st = Math.sin(parameters[5]);

            double px = ( (x - x0) * ct - (y - y0) * st) / parameters[3];
            double py = ( (y - y0) * ct + (x - x0) * st) / parameters[4];

            return parameters[0] * Math.exp( -0.5 * (px * px + py * py)) + parameters[6];
        }
     */

    public double getMahalanobisDistance(double x, double y) {
        double A = parameters[0];
        double x0 = parameters[1];
        double y0 = parameters[2];
        double sx = parameters[3];
        double sy = parameters[4];
        double ct = Math.cos(parameters[5]);
        double st = Math.sin(parameters[5]);

        double px = ( (x - x0) * ct - (y - y0) * st) / sx;
        double py = ( (y - y0) * ct + (x - x0) * st) / sy;

        return Math.sqrt(px * px + py * py);
    }


    public String[] getParameterNames() {
        return new String[] {
            "A", "x0", "y0", "sigX", "sigY", "Theta", "B"};
    }


    public boolean[] getAdjustParameters() {
//                            "A",  "x0", "y0", "sigX", "sigY", "Theta", "B"};
        return new boolean[] {true, true, true, true, true, true, false};
    }

}
