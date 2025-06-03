package edu.ucsc.neurobiology.vision.math.fitting;


/**
 * Implements a fittable 2D generalized gaussian.
 * Parameters: A, x0, y0, sigX, sigY, theta, B (the pedestal).
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class Gaussian2DFunction
    extends FittableFunction {

    public Gaussian2DFunction() {
        setParameters(new double[7]);
    }


    public Gaussian2DFunction(double A, double x0, double y0, double sigX, double sigY,
                              double theta, double B) {

        setParameters(new double[] {A, x0, y0, sigX, sigY, theta, B});

        // added to keep compatibility with old code
        this.setParameterState(6, false);
    }


    /**
     * Scales 
     * @param scaleFactor
     * @return
     */
    public void scaleUp(int scaleFactor) {
        parameters[1] *= scaleFactor;
        parameters[2] *= scaleFactor;
        parameters[3] *= scaleFactor;
        parameters[4] *= scaleFactor;
    }
    

    public void normalizeVariables() {
        // fix sigmas
        parameters[3] = Math.abs(parameters[3]);
        parameters[4] = Math.abs(parameters[4]);

        parameters[5] %= 2 * Math.PI;
        if (parameters[5] < 0) {
            parameters[5] += 2 * Math.PI;
        }
    }


    public double getValueAndDerivatives(
        final double[] _x, final double[] parameters, final double[] derivatives) {

        double x = _x[0];
        double y = _x[1];

        double A1 = parameters[0];
        double x0 = parameters[1];
        double y0 = parameters[2];
        double sx = parameters[3];
        double sy = parameters[4];
        double ct = Math.cos(parameters[5]);
        double st = Math.sin(parameters[5]);
        double B = parameters[6];

        double px = ( (x - x0) * ct - (y - y0) * st) / sx;
        double py = ( (y - y0) * ct + (x - x0) * st) / sy;
        double exp = Math.exp( -0.5 * (px * px + py * py));

        double G = A1 * exp + B;
        double G_A = exp;
        double G_x0 = A1 * exp * (ct * px / sx + st * py / sy);
        double G_y0 = A1 * exp * (ct * py / sy - st * px / sx);
        double G_sx = A1 * exp * px * px / sx;
        double G_sy = A1 * exp * py * py / sy;
        double G_theta = A1 * exp * px * py * (sy / sx - sx / sy);
        double G_B = 1;

        derivatives[0] = G_A;
        derivatives[1] = G_x0;
        derivatives[2] = G_y0;
        derivatives[3] = G_sx;
        derivatives[4] = G_sy;
        derivatives[5] = G_theta;
        derivatives[6] = G_B;

        return G;
    }


    public final double getValue(final double x, final double y) {
        double x0 = parameters[1];
        double y0 = parameters[2];
        double ct = Math.cos(parameters[5]);
        double st = Math.sin(parameters[5]);

        double px = ( (x - x0) * ct - (y - y0) * st) / parameters[3];
        double py = ( (y - y0) * ct + (x - x0) * st) / parameters[4];

        return parameters[0] * Math.exp( -0.5 * (px * px + py * py)) + parameters[6];
    }


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


    /**
     * Removed since this was a crappy way of dealing with params.
     */
    //    public boolean[] getAdjustParameters() {
//                            "A",  "x0", "y0", "sigX", "sigY", "Theta", "B"};
//        return new boolean[] {true, true, true, true, true, true, false};
//    }


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
}
