package edu.ucsc.neurobiology.vision.math.fitting;

import java.util.*;

import edu.ucsc.neurobiology.vision.math.*;


/**
 * The <code>Fitter</code> class implements the Levenberg-Marquardt multidimensional
 * fitting algorithm. Solutions of linear sets of equations are calculated using the
 * Gauss method. The code was ported to Java from the Numerical Recipes (C Version)
 * book. The fitter supports n-dimensional (any n) model functions with any number
 * of free parameters. As the result of the fitting process the best set of parameters
 * is evaluated, the final Chi Squared and the covariance matrix are calculated and
 * stored back in the model function.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class Fitter {
    private static final boolean DEBUG = false;

    private int paramsCount, needAdjustCount, nData;
    private boolean[] needAdjust;
    private double lamda, chi2;
    private double[] y, sigma, params, initialParams, beta, da, derivative, coords;
    private double[][] x, oneDa, covar;
    private FittableFunction function;

    private double deltaChi2 = 1e-4;
    private int maxIterations = 1000;
//    private double deltaChi2 = 1e-8;
//    private int maxIterations = 10000;

    // the Gauss linear equations solver
    private static void gaussj(double[][] a, int n, double[][] b) throws
        FitFailedException {

        final int m = b[0].length;
        int icol = 0, irow = 0;
        int[] indxc = new int[n];
        int[] indxr = new int[n];
        int[] ipiv = new int[n];

        for (int i = 0; i < n; i++) {

            double big = 0, abs;
            for (int j = 0; j < n; j++) {
                if (ipiv[j] != 1) {
                    for (int k = 0; k < n; k++) {
                        if (ipiv[k] == 0) {
                            abs = a[j][k];
                            if (abs < 0) {
                                abs = -abs;
                                //abs = (abs < 0) ? -abs : abs;
                            }
                            if (abs >= big) {
                                big = abs;
                                irow = j;
                                icol = k;
                            }
                        } else if (ipiv[k] > 1) {
                            throw new FitFailedException("gaussj: Singular Matrix-1");
                        }
                    }
                }
            }

            ++ (ipiv[icol]);
            if (irow != icol) {
                for (int l = 0; l < n; l++) {
                    double t = a[irow][l];
                    a[irow][l] = a[icol][l];
                    a[icol][l] = t;
                }
                for (int l = 0; l < m; l++) {
                    double t = b[irow][l];
                    b[irow][l] = b[icol][l];
                    b[icol][l] = t;
                }
            }
            indxr[i] = irow;
            indxc[i] = icol;

            if (a[icol][icol] == 0) {
                throw new FitFailedException("gaussj: Singular Matrix-2");
            }

            double pivinv = 1 / a[icol][icol];
            a[icol][icol] = 1;

            for (int l = 0; l < n; l++) {
                a[icol][l] *= pivinv;
            }
            for (int l = 0; l < m; l++) {
                b[icol][l] *= pivinv;

            }
            for (int ll = 0; ll < n; ll++) {
                if (ll != icol) {
                    double dum = a[ll][icol];
                    a[ll][icol] = 0;
                    for (int l = 0; l < n; l++) {
                        a[ll][l] -= a[icol][l] * dum;
                    }
                    for (int l = 0; l < m; l++) {
                        b[ll][l] -= b[icol][l] * dum;
                    }
                }
            }
        }

        for (int l = n - 1; l >= 0; l--) {
            if (indxr[l] != indxc[l]) {
                for (int k = 0; k < n; k++) {
                    double t = a[k][indxr[l]];
                    a[k][indxr[l]] = a[k][indxc[l]];
                    a[k][indxc[l]] = t;
                }
            }
        }
    }


    /**
     * This method fills the arrays Alpha and Beta and calculates and returns
     * the Chi-Square. Alpha and Beta can then be used to solve the
     * Levenberg-Marquardt set of linear equations w.r.t to da.
     */
    private double prepareMatrices(double[][] alpha, double[] beta) throws
        CannotEvaluateException {
        int i, j, k, l, m;
        double ymod, wt, sig2i, dy, chi2 = 0;

        for (j = 0; j < needAdjustCount; j++) {
            for (k = 0; k <= j; k++) {
                alpha[j][k] = 0;
            }
            beta[j] = 0;
        }

        for (i = 0; i < nData; i++) {
            for (int n = 0; n < x.length; n++) {
                coords[n] = x[n][i];
            }
            ymod = function.getValueAndDerivatives(coords, params, derivative);
            sig2i = 1 / (sigma[i] * sigma[i]);
            dy = y[i] - ymod;
            chi2 += dy * dy * sig2i;
            for (j = -1, l = 0; l < paramsCount; l++) {
                if (needAdjust[l]) {
                    wt = derivative[l] * sig2i;
                    for (j++, k = 0, m = 0; m <= l; m++) {
                        if (needAdjust[m]) {
                            alpha[j][k] += wt * derivative[m];
                            k++;
                        }
                    }
                    beta[j] += dy * wt;
                }
            }
        }

        for (j = 1; j < needAdjustCount; j++) {
            for (k = 0; k < j; k++) {
                alpha[k][j] = alpha[j][k];

            }
        }
        return chi2;
    }


    // the initialization step
    private void initializeFit(double[][] alpha) throws CannotEvaluateException {
        params = new double[paramsCount];
        System.arraycopy(initialParams, 0, params, 0, paramsCount);

        needAdjustCount = 0;
        for (int j = 0; j < paramsCount; j++) {
            if (needAdjust[j]) {
                needAdjustCount++;

            }
        }
        oneDa = new double[needAdjustCount][1];
        beta = new double[paramsCount];
        da = new double[paramsCount];

        lamda = 0.001;
        chi2 = prepareMatrices(alpha, beta);
    }


    // the finalization step
    private void finishFit(double[][] covar, double[][] alpha) throws FitFailedException {

        int j, k, l, m;

        for (j = -1, l = 0; l < paramsCount; l++) {
            if (needAdjust[l]) {
                for (j++, k = 0, m = 0; m < paramsCount; m++) {
                    if (needAdjust[m]) {
                        covar[j][k] = alpha[j][k];
                        k++;
                    }
                }
                covar[j][j] = alpha[j][j] * (1.0 + lamda);
                oneDa[j][0] = beta[j];
            }
        }

        gaussj(covar, needAdjustCount, oneDa);
        for (j = 0; j < needAdjustCount; j++) {
            da[j] = oneDa[j][0];

        }
        covsrt(covar, needAdjust, needAdjustCount);
    }


    // this method performs one iteration
    private void mrqmin(double[][] covar, double[][] alpha) throws FitFailedException,
        CannotEvaluateException {

        int j, k, l, m;

        // Compute Alpha'[i,j] (the left hand side of the linear equation)
        // and fill the right hand side vector. These two will be passed to
        // the linear equations solver to get the parameter variations.
        for (j = -1, l = 0; l < paramsCount; l++) {
            if (needAdjust[l]) {
                for (j++, k = 0, m = 0; m < paramsCount; m++) {
                    if (needAdjust[m]) {
                        covar[j][k] = alpha[j][k];
                        k++;
                    }
                }
                covar[j][j] = alpha[j][j] * (1.0 + lamda);
                oneDa[j][0] = beta[j];
            }
        }

        // Run the linear equation solver to get the parameter variations.
        // Copy the parameter variations (oneDa) into the da array.
        gaussj(covar, needAdjustCount, oneDa);
        for (j = 0; j < needAdjustCount; j++) {
            da[j] = oneDa[j][0];
        }

        // update the values of the running parameters according to da.
        for (j = 0, l = 0; l < paramsCount; l++) {
            if (needAdjust[l]) {
                params[l] = initialParams[l] + da[j];
                j++;
            }
        }

        //
        double newChi2 = prepareMatrices(covar, da);
        if (newChi2 < chi2) {
            lamda *= 0.1;
            chi2 = newChi2;
            for (j = -1, l = 0; l < paramsCount; l++) {
                if (needAdjust[l]) {
                    for (j++, k = 0, m = 0; m < paramsCount; m++) {
                        if (needAdjust[m]) {
                            alpha[j][k] = covar[j][k];
                            k++;
                        }
                    }
                    beta[j] = da[j];
                    initialParams[l] = params[l];
                }
            }
        } else {
            lamda *= 10.0;
        }
    }


    // this rearanges the covariance matrix at the end of the fitting process
    private static void covsrt(double[][] covar, boolean needAdjust[],
                               int needAdjustCount) {
        for (int i = needAdjustCount + 1; i < needAdjust.length; i++) {
            for (int j = 1; j <= i; j++) {
                covar[i][j] = covar[j][i] = 0;

            }
        }
        int k = needAdjustCount - 1;
        for (int j = needAdjust.length - 1; j >= 0; j--) {
            if (needAdjust[j]) {
                for (int i = 0; i < needAdjust.length; i++) {
                    double swap = covar[i][k];
                    covar[i][k] = covar[i][j];
                    covar[i][j] = swap;
                }
                for (int i = 0; i < needAdjust.length; i++) {
                    double swap = covar[k][i];
                    covar[k][i] = covar[j][i];
                    covar[j][i] = swap;
                }
                k--;
            }
        }
    }


    /**
     * Returns the value of the minimized Chi-Squared after the fit.
     */
    /*
         public double getChiSquared() {
        return chi2;
         }
     */

    /**
     * Returns the resulting covarince matrix of the parameters as two dimesional array
     * in row-major format.
     */
    public double[][] getCovarianceMatrix() {
        return covar;
    }


    /**
     * Returns the set of parameters that minimize the Chi-Squared.
     */
    /*
         public double[] getParameters() {
        return (double[])initialParams.clone();
         }
     */

    /**
     * This method allows to set the precision of the Chi-Squared at which the fitting
     * is considered done.
     *
     * @param deltaChi2 the minimum varation between two consecutive iterations at
     *        which the fit stops
     */
    public void setLastChiSquaredVariation(double deltaChi2) {
        if (deltaChi2 < 0) {
            throw new IllegalArgumentException("The Chi2 variation must be positive");
        }

        this.deltaChi2 = deltaChi2;
    }


    /**
     * This method allows to set the total number of iterations made by the fit
     * funtion. If the number of iterations is exceeded a FitFailedException gets
     * thrown by the fit() methods. Default value is 1000.
     *
     * @param maxIterations the minimum number of iterations
     */
    public void setMaxIterations(int maxIterations) {
        if (maxIterations < 1) {
            throw new IllegalArgumentException("maxIterations must be at least 1");
        }

        this.maxIterations = maxIterations;
    }


    /**
     * Performs a Levenberg-Marquardt multiparameter fit of the given multidimensional
     * function to the data set described by the n-dimensional coordinates <code>x</code>
     * and the measured values <code>y</code> with errors <code>sigma</code>. The
     * results of the fit (Chi-Squared, parameters, covariance matrix) are set
     * into the model funtion.
     *
     * @param function The model function to be fitted
     * @param x The n-dimensional coordinates
     * @param y The measured values
     * @param sigma The errors of the measurements
     * @throws FitFailedException In the case the fitting process fails due to
     *         various reasons (is a checked exception).
     * @throws InvalidFitDataException In the case the provided data are invalid
     *         (is NOT a checked exception).
     */
    public double fit(
        FittableFunction function, double[][] x, double[] y, double[] sigma) throws
        FitFailedException, InvalidFitDataException {

        return _fit(function, x, y, sigma, y.length);
    }


    /**
     * Performs a Levenberg-Marquardt multiparameter fit of the given multidimensional
     * function to the data set described by the n-dimensional coordinates <code>x</code>
     * and the measured values <code>y</code> with errors <code>sigma</code>. The
     * results of the fit (Chi-Squared, parameters, covariance matrix) are set
     * into the model funtion.
     *
     * @param function The model function to be fitted
     * @param x The n-dimensional coordinates
     * @param y The measured values
     * @param sigma The errors of the measurements
     * @param nData the number of actual data points contained in <b>y</b> and <b>x</b>.
     * @throws FitFailedException In the case the fitting process fails due to
     *         various reasons (is a checked exception).
     * @throws InvalidFitDataException In the case the provided data are invalid
     *         (is NOT a checked exception).
     */
    public double fit(
        FittableFunction function, double[][] x, double[] y, double[] sigma, int nData) throws
        FitFailedException, InvalidFitDataException {

        return _fit(function, x, y, sigma, nData);
    }


    public static double fit1D(
        FittableFunction function, double[] x, double[] y, double[] sigma, int nData) throws
        FitFailedException, InvalidFitDataException {

        return new Fitter()._fit(function, new double[][] {x}, y, sigma, nData);
    }


    private double _fit(
        FittableFunction function, double[][] x, double[] y, double[] sigma, int nData) throws
        FitFailedException, InvalidFitDataException {

        if (x.length < 1) {
            throw new InvalidFitDataException("x.length < 1");
        }
        for (int n = 0; n < x.length; n++) {
            if ( (x[n].length != y.length) || (x[n].length != sigma.length)) {
                throw new InvalidFitDataException("Inconsistent length of data arrays");
            }
        }
        this.x = x;
        this.y = y;
        this.sigma = sigma;
        this.nData = nData;

        this.function = function;
        this.initialParams = function.getParameters();
        this.paramsCount = initialParams.length;
        derivative = new double[paramsCount];
        coords = new double[x.length];

        covar = new double[paramsCount][paramsCount];
        double[][] alpha = new double[paramsCount][paramsCount];

        needAdjust = function.getAdjustParameters();
        if (needAdjust == null) {
            needAdjust = new boolean[paramsCount];
            Arrays.fill(needAdjust, true);
        }

        double oldChi2;
        try {
            initializeFit(alpha);
        } catch (CannotEvaluateException e) {
            throw new FitFailedException("cannot initialize fit", e);
        }

        for (int iteration = 1; ; iteration++) {
            if (iteration > maxIterations) {
                throw new FitFailedException(
                    "Maximum number of iterations (" + maxIterations + ") exceeded.");
            }

            if (DEBUG) {
                System.out.println("Iteration " + iteration + ", lamda = " + lamda);
                System.out.println("-----------------------------------------");
                function.setParameters(params);
                function.setChiSquared(chi2);
                System.out.println(function);
            }

            oldChi2 = chi2;
//            System.out.println("calling mrqmin");
            try {
                mrqmin(covar, alpha);
            } catch (CannotEvaluateException e) {
                throw new FitFailedException("", e);
            }
//            System.out.println("done mrqmin");

            if (chi2 == oldChi2) {
                continue;
            }
            if (Math.abs(chi2 - oldChi2) < deltaChi2) {
//                System.out.println("iterations  " + iteration);
                break;
            }
        }

        finishFit(covar, alpha);

        function.setParameters(initialParams);
        function.setChiSquared(chi2 / (nData - paramsCount));
        function.setCovarianceMatrix(covar);

        return function.getChiSquared();
    }

}
