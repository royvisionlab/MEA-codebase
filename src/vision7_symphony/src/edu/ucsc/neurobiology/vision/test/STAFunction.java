package edu.ucsc.neurobiology.vision.test;

import edu.ucsc.neurobiology.vision.math.fitting.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.stimulus.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class STAFunction
    extends FittableFunction {

    private final double n1, n2;


    public STAFunction(int n1, int n2) {
        this.n1 = n1;
        this.n2 = n2;
    }


    public boolean[] getAdjustParameters() {
        return new boolean[] {
            false, true, true, true, true, true,
            true, true, true, true};
    }


    public STAFunction(
        double A, double x0, double y0, double sigX, double sigY, double theta,
        double a1, double t1, double a2, double t2, int n1, int n2) {

        this(n1, n2);
        setParameters(new double[] {A, x0, y0, sigX, sigY, theta, a1, t1, a2, t2});
    }


    public STATimeFunction1 getTimeFunction() {
        double A = parameters[0];

        double a1 = parameters[6];
        double t1 = parameters[7];
        double a2 = parameters[8];
        double t2 = parameters[9];

        return new STATimeFunction1(A * a1, t1, A * a2, t2, n1, n2);
    }


    public ParametricEllipse getEllipse(double pixelSize) {
        double A = parameters[0];
        double x0 = parameters[1];
        double y0 = parameters[2];
        double sx = parameters[3];
        double sy = parameters[4];
        double ct = Math.cos(parameters[5]);
        double st = Math.sin(parameters[5]);

        return new ParametricEllipse(
            x0 * pixelSize, y0 * pixelSize, sx * pixelSize, sy * pixelSize, -parameters[5]);
    }


    public double getValueAndDerivatives(
        final double[] _x, final double[] parameters, final double[] derivatives) {

        double A = parameters[0];
        double x0 = parameters[1];
        double y0 = parameters[2];
        double sx = parameters[3];
        double sy = parameters[4];
        double ct = Math.cos(parameters[5]);
        double st = Math.sin(parameters[5]);

        double a1 = parameters[6];
        double t1 = parameters[7];
        double a2 = parameters[8];
        double t2 = parameters[9];

        double x = _x[0];
        double y = _x[1];
        double t = _x[2];

        double p1 = Math.pow( (t / t1) * Math.exp(1 - t / t1), n1);
        double p2 = Math.pow( (t / t2) * Math.exp(1 - t / t2), n2);

        double T = a1 * p1 + a2 * p2;
        double T_a1 = p1;
        double T_t1 = a1 * n1 * p1 * (t / t1 - 1) / t1;
        double T_a2 = p2;
        double T_t2 = a2 * n2 * p2 * (t / t2 - 1) / t2;

        double px = ( (x - x0) * ct - (y - y0) * st) / sx;
        double py = ( (y - y0) * ct + (x - x0) * st) / sy;
        double xx = -0.5 * (px * px + py * py);
        double exp = Math.exp(xx);

        double G = A * exp;
        double G_A = exp;
        double G_x0 = A * exp * (ct * px / sx + st * py / sy);
        double G_y0 = A * exp * (ct * py / sy - st * px / sx);
        double G_sx = A * exp * px * px / sx;
        double G_sy = A * exp * py * py / sy;
        double G_theta = A * exp * px * py * (sy / sx - sx / sy);

        derivatives[0] = T * G_A;
        derivatives[1] = T * G_x0;
        derivatives[2] = T * G_y0;
        derivatives[3] = T * G_sx;
        derivatives[4] = T * G_sy;
        derivatives[5] = T * G_theta;

        derivatives[6] = G * T_a1;
        derivatives[7] = G * T_t1;
        derivatives[8] = G * T_a2;
        derivatives[9] = G * T_t2;

        return G * T;
    }


    /*
        public double getA1() {
            return parameters[0];
        }


        public double getA2() {
            return parameters[2];
        }


        public double getT1() {
            return parameters[1];
        }


        public double getT2() {
            return parameters[3];
        }


        public String getDescription() {
            return "STATimeFunction";
        }
     */


    public static int convertFrame(
        STA sta, double[][] x, double[] y, double[] sigma, int colorIndex) {

        final int w = sta.getWidth(), h = sta.getHeight();
        int n = 0;
        int staDepth = sta.size();

        int x1 = 0;
        int x2 = w - 1;
        int y1 = 0;
        int y2 = h - 1;
        double average = 0;
        for (int f = 0; f < staDepth; f++) {
            ImageFrame frame = sta.getFrame(f);
            for (int i = x1; i <= x2; i++) {
                for (int j = y1; j <= y2; j++) {
                    // the x and y coordinates
                    x[0][n] = i + 0.5;
                    x[1][n] = h - j - 1 + 0.5;
                    x[2][n] = - (staDepth - f - 1) * 8.34;

                    // the value, never use abs() here
                    y[n] = frame.getPixel(i, j, colorIndex);

                    average += y[n];

                    // the error
                    sigma[n] = frame.getPixelError(i, j, colorIndex);

                    n++;
                }
            }
        }

        average /= y.length;
        for (int i = 0; i < y.length; i++) {
            y[i] -= average;
        }

        return n;
    }

}
