package edu.ucsc.neurobiology.vision.math.expressions;

import java.util.*;

import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.math.fitting.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * A function that can be defined by an expression in the Vision Embedded Language.
 * For an example use see testing.MathTest.testArbitraryFunction().
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ArbitraryFunction
    extends FittableFunction implements FunctionData {

    private String expression;
    private String[] paramNames, variableNames;
    private HashMap<String, Double> parametersTable = new HashMap<String, Double>();


    public ArbitraryFunction(String expression, String vNames, String pNames,
                             double ...paramValues) {

        this.paramNames = StringUtil.decomposeString(pNames, " ");
        this.variableNames = StringUtil.decomposeString(vNames, " ");

        this.expression = expression;

        this.setParameters(paramValues);
    }


    public void setParameters(double[] p) {
        super.setParameters(p);
        for (int i = 0; i < paramNames.length; i++) {
            parametersTable.put(paramNames[i], p[i]);
        }
    }


    public double get(String pName) {
        return parametersTable.get(pName);
    }


//    public void setParameter(String param, double value) {
//        parametersTable.put(param, new Double(value));
//    }


    public double getValueAt(double x) throws CannotEvaluateException {
        parametersTable.put("x", x);
        for (int i = 0; i < paramNames.length; i++) {
            parametersTable.put(paramNames[i], parameters[i]);
        }
        return getValueAt(expression, parametersTable);
    }


    public static double getValueAt(String expression, HashMap<String,Double> params) throws
        CannotEvaluateException {

        try {
            ExpressionScanner scanner = new ExpressionScanner(expression);
            ExpressionParser p = new ExpressionParser(scanner, params);
            p.parse();
            return (Double) p.result;
        } catch (Exception e) {
            throw new CannotEvaluateException(e);
        }
    }


    public double getValueAndDerivatives(
        final double[] x, final double[] parameters, final double[] derivatives) throws
        CannotEvaluateException {

        // copy params
        for (int i = 0; i < x.length; i++) {
            parametersTable.put(variableNames[i], x[i]);
        }
        for (int i = 0; i < paramNames.length; i++) {
            parametersTable.put(paramNames[i], parameters[i]);
        }

        // calculate value
        double f = getValueAt(expression, parametersTable);
        parametersTable.put("f", f);

        // calculate derivatives
        double dx = 0.00001;
        for (int i = 0; i < derivatives.length; i++) {
//            derivatives[i] = getValueAt(derivativeExpression[i], parametersTable);

            parametersTable.put(paramNames[i], parametersTable.get(paramNames[i]) + dx);
            double f2 = getValueAt(expression, parametersTable);
            parametersTable.put(paramNames[i], parametersTable.get(paramNames[i]) - dx);

            derivatives[i] = (f2 - f) / dx;
        }

        return f;
    }


    public String getDescription() {
        return expression;
    }


}
