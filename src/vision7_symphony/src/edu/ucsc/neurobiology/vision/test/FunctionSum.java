package edu.ucsc.neurobiology.vision.test;

import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.math.fitting.*;
import edu.ucsc.neurobiology.vision.plot.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class FunctionSum
    extends FittableFunction implements FunctionData {

    private FittableFunction[] functions;
    private double[][] params, derivs;
    private double[] values;


    public FunctionSum(FittableFunction ...functions) {
        this.functions = functions;

        params = new double[functions.length][];
        derivs = new double[functions.length][];
        values = new double[functions.length];
        int nParams = 0;
        for (int i = 0; i < functions.length; i++) {
            params[i] = functions[i].getParameters();
            derivs[i] = new double[functions[i].getParameters().length];
            nParams += functions[i].getParameters().length;
        }

        double[] parameters = new double[nParams];
        for (int i = 0, fIndex = 0, pIndex = -1; i < parameters.length; i++) {
            pIndex++;
            if (pIndex == params[fIndex].length) {
                fIndex++;
                pIndex = 0;
            }
            parameters[i] = params[fIndex][pIndex];
        }
        setParameters(parameters);

//        System.out.println("params " + nParams);
    }


    public double getValueAndDerivatives(
        final double[] x, final double[] parameters, final double[] derivatives) throws
        CannotEvaluateException {

        for (int i = 0, fIndex = 0, pIndex = -1; i < parameters.length; i++) {
            pIndex++;
            if (pIndex == params[fIndex].length) {
                fIndex++;
                pIndex = 0;
            }
            params[fIndex][pIndex] = parameters[i];
        }

        for (int i = 0; i < functions.length; i++) {
            values[i] = functions[i].getValueAndDerivatives(x, params[i], derivs[i]);
        }

        for (int i = 0, fIndex = 0, pIndex = -1; i < parameters.length; i++) {
            pIndex++;
            if (pIndex == params[fIndex].length) {
                fIndex++;
                pIndex = 0;
            }
            derivatives[i] = derivs[fIndex][pIndex];
        }

        return MathUtil.sum(values);
    }


    public double getValueAt(double x) throws CannotEvaluateException {
        // copy parameters to individual groups
        for (int i = 0, fIndex = 0, pIndex = -1; i < parameters.length; i++) {
            pIndex++;
            if (pIndex == params[fIndex].length) {
                fIndex++;
                pIndex = 0;
            }
            params[fIndex][pIndex] = parameters[i];
        }

        // sum all functions
        for (int i = 0; i < functions.length; i++) {
            values[i] = functions[i].getValueAndDerivatives(
                new double[] {x}, params[i], derivs[i]);
        }

        return MathUtil.sum(values);
    }


    public String getDescription() {
        return "";
    }
}
