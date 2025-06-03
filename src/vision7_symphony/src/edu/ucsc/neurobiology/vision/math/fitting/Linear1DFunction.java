package edu.ucsc.neurobiology.vision.math.fitting;

import edu.ucsc.neurobiology.vision.plot.*;


/**
 * Implemets a fittable linear function.
 * Parameters: "slope", "intercept".
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class Linear1DFunction
    extends FittableFunction implements FunctionData {


    public Linear1DFunction(double k, double b) {
        setParameters(new double[] {k, b});
    }


    public Linear1DFunction() {
        setParameters(new double[2]);
    }


    public String[] getParameterNames() {
        return new String[] {
            "slope", "intercept"};
    }


    public double getValueAndDerivatives(
        final double[] x, final double[] parameters, final double[] derivatives) {

        derivatives[0] = x[0];
        derivatives[1] = 1;

        return parameters[0] * x[0] + parameters[1];
    }


    public double getSlope() {
        return parameters[0];
    }


    public double getIntercept() {
        return parameters[1];
    }


    public double getValueAt(double x) {
        return parameters[0] * x + parameters[1];
    }


    public String getDescription() {
        return "Linear 1D Function";
    }


    public static Linear1DFunction lineFit(double[] x, double[] y, double[] sigma) {
        double s = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
        for (int i = 0; i < sigma.length; i++) {
            s += 1 / (sigma[i] * sigma[i]);
            sx += x[i] / (sigma[i] * sigma[i]);
            sy += y[i] / (sigma[i] * sigma[i]);
            sxx += (x[i] * x[i]) / (sigma[i] * sigma[i]);
            sxy += (x[i] * y[i]) / (sigma[i] * sigma[i]);
        }
        double delta = s * sxx - sx * sx;
        double intercept = (sxx * sy - sx * sxy) / delta;
        double slope = (s * sxy - sx * sy) / delta;
        double sigmaIntercept2 = sxx / delta;
        double sigmaSlope2 = s / delta;

        double chi2 = 0;
        for (int i = 0; i < sigma.length; i++) {
            double a = (y[i] - intercept - slope * x[i]) / sigma[i];
            chi2 += a * a;
        }

        Linear1DFunction l = new Linear1DFunction();
        l.setParameters(new double[] {slope, intercept});
        l.setChiSquared(chi2);
        l.setCovarianceMatrix(new double[][] { {sigmaSlope2, 0}
                              , {0, sigmaIntercept2}
        });
        return l;
    }

}
