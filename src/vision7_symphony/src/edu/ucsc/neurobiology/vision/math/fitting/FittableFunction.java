package edu.ucsc.neurobiology.vision.math.fitting;

import java.io.*;
import java.text.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.math.*;


/**
 * This class represents the model function that needs to be passed to the
 * <code>Fitter</code> class along with the data to be fit. This is an abstract class
 * that implements some basic functionality needed for any function. The concrete
 * implementations must implement the abstract <code>getValueAndDerivatives()</code>
 * function to make the implementation useful.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public abstract class FittableFunction {
    private static DecimalFormat df1 = new DecimalFormat();

    protected double[] parameters;
    protected double chiSquared = Double.NaN;
    protected double[][] covarianceMatix;
    protected int nParameters = -1;
    boolean[] adjustParameters;

    /**
     * Used by the fitter to set the chi squared after a fit.
     * @param chiSquared double
     */
    public void setChiSquared(double chiSquared) {
        this.chiSquared = chiSquared;
    }


    /**
     * Returns the parameters indexed by "i" along with its error.
     * @param i int
     * @return Num
     */
    public Num param(int i) {
        return new Num(parameters[i], Math.sqrt(covarianceMatix[i][i]));
    }


    /**
     * Returns the parameters indexed by "i".
     * @param i int
     * @return double
     */
    public double getParameter(int i) {
        return parameters[i];
    }


    /**
     * Returns the error of the parameter indexed by "i".
     * @param i int
     * @return double
     */
    public double getParameterError(int i) {
        return Math.sqrt(covarianceMatix[i][i]);
    }


    /**
     * Sets the error of the parameter indexed by "i".
     * @param i int
     * @return double
     */
    public void setParameterError(int i, double err) {
        covarianceMatix[i][i] = err * err;
    }


    /**
     * Sets the value and error of the parameter indexed by "i".
     * @param i int
     * @return double
     */
    public void setParameter(int i, double value, double err) {
        parameters[i] = value;
        if (covarianceMatix == null) {
            covarianceMatix = new double[nParameters][nParameters];
        }
        covarianceMatix[i][i] = err * err;
    }


    public void setParameter(int i, double v) {
        parameters[i] = v;
    }


    /**
     * This method conerts a double value to a String keeping only the
     * specified number of fractionar digits.
     */
    protected static String format(double value, int nFractionarDigits) {
        df1.setMinimumFractionDigits(nFractionarDigits);
        df1.setMaximumFractionDigits(nFractionarDigits);
        return df1.format(value);
    }


    /**
     * Returns the Chi-Squared the funtion got after the fit. If called before the fit
     * or if the fit failed this will return Double.NaN.
     */
    public double getChiSquared() {
        return chiSquared;
    }


    /**
     * used by the fitter to write back the error matrix.
     * @param covarianceMatix double[][]
     */
    public void setCovarianceMatrix(double[][] covarianceMatix) {
        double[][] newMatrix = new double[nParameters][nParameters];
        for (int i = 0; i < nParameters; i++) {
            if (covarianceMatix[i].length != nParameters) {
                throw new IllegalArgumentException("Wromg size of Covariance Matrix");
            }
            System.arraycopy(covarianceMatix[i], 0, newMatrix[i], 0, nParameters);
        }

        this.covarianceMatix = newMatrix;
    }


    /**
     * Returns the covarince matrix of the parameters the funtion got after the fit.
     * If called before the fit or if the fit failed this will return <b>null</b>.
     */
    public double[][] getCovarianceMatrix() {
        if (covarianceMatix == null) {
            return null;
        }

        double[][] newMatrix = new double[nParameters][nParameters];
        for (int i = 0; i < nParameters; i++) {
            System.arraycopy(covarianceMatix[i], 0, newMatrix[i], 0, nParameters);
        }

        return newMatrix;
    }


    /**
     * Returns the parameters array. The method will never return null. After a fit this
     * funtion will return the final values of parameters.
     */
    public double[] getParameters() {
        return parameters;
    }


    public int getParametersCount() {
        return nParameters;
    }


    /**
     * If the state of a parameter in "true" it is considered for adjustment by the fitter.
     * It is held contant if the state is "false".
     * @param param int
     * @param adjust boolean
     */
    public void setParameterState(int param, boolean adjust) {
        adjustParameters[param] = adjust;
    }


    /**
     * Used to set the value of parameters. This method MUST be called before the fit
     * happens to properly set the starting values of parameters. If the fit is
     * successful the Fitter will set the final values using this function.
     */
    public void setParameters(double[] p) {
        if (nParameters == -1) {
            // if the first set
            this.parameters = new double[p.length];
            this.nParameters = p.length;
            adjustParameters = new boolean[this.nParameters];
            Arrays.fill(adjustParameters, true);
        } else {
            if (p.length != nParameters) {
                throw new IllegalArgumentException(
                    "Wrong number of parameters: " + p.length +
                    " received but " + nParameters + " required");
            }
        }

        System.arraycopy(p, 0, this.parameters, 0, nParameters);
    }


    /**
     * Returns a boolean array of the same size as the parameters array. If an array
     * element contains <b>true</b> then the coresponding parameter will be changed
     * during the fit, otherwise the value of the parameter will be considered fixed.
     * If the function returns <b>null</b> all the parameters will be changed. This
     * implementation returns <b>null</b>.
     */
    public boolean[] getAdjustParameters() {
        return adjustParameters;
    }


    /**
     * This method MUST be properly implemented by subclasses to provide the value
     * and derivatives (w.r.t parameters) of the model function.
     *
     * @param x the vector of coordinates
     * @param parameters the vector of current parameters than should be used.
     *        Do NOT change the values of parameters
     * @param derivatives the vector in which the calculated derivatives have to be placed
     * @return the value of the function
     */
    public abstract double getValueAndDerivatives(
        final double[] x, final double[] parameters, final double[] derivatives) throws
        CannotEvaluateException;


    /**
     * This function can be overritten by the subclasses to return the parameter names
     * in an array of strings. This implementation return null.
     */
    public String[] getParameterNames() {
        return null;
    }


    /**
     * Returns a clone of this function.
     */
    public Object clone() {
        try {
            return super.clone();
        } catch (CloneNotSupportedException e) {
            return null;
        }
    }


//    public Tag read(int tagId, TaggedInputStream input, int len) throws IOException {
//        this.nParameters = input.readInt();
//        parameters = new double[nParameters];
//
//        for (int i = 0; i < nParameters; i++)
//            parameters[i] = input.readDouble();
//
//        return (Tag)clone();
//    }
//
//
//    public void write(int tagId, TaggedOutputStream output) throws IOException {
//        output.writeInt(nParameters);
//        for (int i = 0; i < parameters.length; i++)
//            output.writeDouble(parameters[i]);
//    }


    /**
     * Returns a string representation of the parameters and their coresponding errors
     * (if they exist). Never returns null or empty strings.
     */
    public String toString() {
        String[] paramNames = getParameterNames();
        String s = this.getClass().getName() + ", Chi2: " + chiSquared + "\n";

        for (int i = 0; i < parameters.length; i++) {
            if (paramNames != null) {
                s += paramNames[i] + " = " + format(parameters[i], 4);
            } else {
                s += "parameters[" + i + "] = " + format(parameters[i], 4);
            }

            if (covarianceMatix != null) {
                double sigma = Math.sqrt(covarianceMatix[i][i]);
                double percent = sigma * 100 / parameters[i];

                s += " +/- " + format(sigma, 4) + " (" + format(percent, 4) + "%)\n";
            } else {
                s += "\n";
            }
        }

        return s;
    }


    /**
     * Returns a string representation of the covariance matrix associated with
     * this function. If called before the fit heppened or if the fit failed then
     * this will return an empty string.
     */
    public String covarianceToString() {
        if (covarianceMatix != null) {
            String s = "        ";
            String[] paramNames = getParameterNames();

            for (int i = 0; i < paramNames.length; i++) {
                s += paramNames[i] + "         ";
            }

            for (int i = 0; i < covarianceMatix.length; i++) {
                s += "\n" + paramNames[i] + ":";
                for (int j = 0; j < covarianceMatix[i].length; j++) {
                    double sigI = Math.sqrt(covarianceMatix[i][i]);
                    double sigJ = Math.sqrt(covarianceMatix[j][j]);
                    double rho = covarianceMatix[i][j] / (sigI * sigJ);
                    s += "  " + format(rho, 3);
                }
            }
            return s;
        } else {
            return "";
        }
    }

}
