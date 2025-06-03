package edu.ucsc.neurobiology.vision.math.expressions;

import edu.ucsc.neurobiology.vision.math.*;


/**
 * This utility class defines the semantics of the Vision Embedded Language (VEL).
 * There is a method for almost every operation or function call defined in the VEL.
 * Each methods does type checking and then performs the required operation. Called
 * from MathExpressions.cup
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
class ExpressionParserSemantics {
    private static ExpressionParser parser;


    /**
     * This is a fully static class. The constuctor is hidden.
     */
    private ExpressionParserSemantics() {
    }


    public static final void setParser(ExpressionParser parser) {
        ExpressionParserSemantics.parser = parser;
    }


    public static final Object NORM(Object e) throws Exception {
  
        if (e instanceof double[]) {
            double[] a = (double[]) e;
            double norm = 0;
            for (int i = 0; i < a.length; i++) {
                norm += a[i]*a[i];
            }
            
            norm = Math.sqrt(norm);
            
            for (int i = 0; i < a.length; i++) {
                a[i] /= norm;
            }
            double[] b = new double[a.length];
            System.arraycopy(a, 0, b, 0, a.length);
            return b;
        } else {
            throw new Exception("Wrong method usage: norm");
        }
    }


    public static final Object NORMRMS(Object e) throws Exception {
        

//        if (e instanceof double[]) {
//            double[] a = (double[]) e;
//            MeanVarianceCalculator mvc = new MeanVarianceCalculator(
//                MeanVarianceCalculator.BIASED);
//            for (int i = 0; i < a.length; i++) {
//                mvc.add(a[i]);
//            }
//            MathUtil.divide(a, mvc.getStandardDeviation());
//            return a;
//        } else {
        //implementation makes no sense.  Stub left here as an 
        //example for the future.  --mgivich
            throw new Exception("Method not supported: normrms");
 //       }
    }


    public static final Object FFT(Object array, Object n) throws Exception {
        if (array instanceof double[] && n instanceof Double) {
            double[] a = (double[]) array;
            //int n = (int) Math.pow(2, Math.ceil(Math.log(a.length) / Math.log(2)));
            double[] re = new double[ ( (Double) n).intValue()];
            double[] im = new double[ ( (Double) n).intValue()];
            System.arraycopy(a, 0, re, 0, a.length);
            FFT.fft(re, im, +1);
            return FFT.spectrum(re, im);
        } else {
            throw new Exception("Wrong method usage: fft");
        }
    }


    public static final Object ADD(Object e1, Object e2) throws Exception {
        if (e1 instanceof Double && e2 instanceof Double) {
            return (Double) e1 + (Double) e2;
        } else if (e1 instanceof double[] && e2 instanceof Double) {
            double[] ee1 = (double[]) e1;
            double ee2 = (Double) e2;
            double[] result = new double[ee1.length];
            for (int i = 0; i < result.length; i++) {
                result[i] = ee1[i] + ee2;
            }
            return result;
        } else if (e1 instanceof Double && e2 instanceof double[]) {
            double ee1 = (Double) e1;
            double[] ee2 = (double[]) e2;
            double[] result = new double[ee2.length];
            for (int i = 0; i < result.length; i++) {
                result[i] = ee1 + ee2[i];
            }
            return result;
        } else if (e1 instanceof double[] && e2 instanceof double[]) {
            double[] ee1 = (double[]) e1;
            double[] ee2 = (double[]) e2;
            if (ee1.length != ee2.length) {
                throw new Exception("Different array length");
            }
            double[] result = new double[ee1.length];
            for (int i = 0; i < result.length; i++) {
                result[i] = ee1[i] + ee2[i];
            }
            return result;
        } else if (e1 instanceof String && e2 instanceof String) {
            return (String) e1 + (String) e2;
        } else {
            throw new Exception(
                "Wrong operator usage: " +
                e1.getClass().getName() + " + " + e2.getClass().getName(), null);
        }
    }


    public static final Object JOIN(Object e1, Object e2) throws Exception {
        if (e1 instanceof double[] && e2 instanceof double[]) {
            double[] ee1 = (double[]) e1;
            double[] ee2 = (double[]) e2;

            double[] result = new double[ee1.length + ee2.length];
            System.arraycopy(ee1, 0, result, 0, ee1.length);
            System.arraycopy(ee2, 0, result, ee1.length, ee2.length);
            return result;
        } else {
            throw new Exception(
                "Wrong operator usage: " +
                e1.getClass().getName() + " # " + e2.getClass().getName(), null);
        }
    }


    public static final Object SUB(Object e1, Object e2) throws Exception {
        if (e1 instanceof Double && e2 instanceof Double) {
            return (Double) e1 - (Double) e2;
        } else if (e1 instanceof double[] && e2 instanceof Double) {
            double[] ee1 = (double[]) e1;
            double ee2 = (Double) e2;
            double[] result = new double[ee1.length];
            for (int i = 0; i < result.length; i++) {
                result[i] = ee1[i] - ee2;
            }
            return result;
        } else if (e1 instanceof Double && e2 instanceof double[]) {
            double ee1 = (Double) e1;
            double[] ee2 = (double[]) e2;
            double[] result = new double[ee2.length];
            for (int i = 0; i < result.length; i++) {
                result[i] = ee1 - ee2[i];
            }
            return result;
        } else if (e1 instanceof double[] && e2 instanceof double[]) {
            double[] ee1 = (double[]) e1;
            double[] ee2 = (double[]) e2;
            if (ee1.length != ee2.length) {
                throw new Exception("Different array length");
            }
            double[] result = new double[ee1.length];
            for (int i = 0; i < result.length; i++) {
                result[i] = ee1[i] - ee2[i];
            }
            return result;
        } else {
            throw new Exception(
                "Wrong operator usage: " +
                e1.getClass().getName() + " - " + e2.getClass().getName(), null);
        }
    }


    public static final Object MUL(Object e1, Object e2) throws Exception {
        if (e1 instanceof Double && e2 instanceof Double) {
            return (Double) e1 * (Double) e2;
        } else if (e1 instanceof double[] && e2 instanceof Double) {
            double[] ee1 = (double[]) e1;
            double ee2 = (Double) e2;
            double[] result = new double[ee1.length];
            for (int i = 0; i < result.length; i++) {
                result[i] = ee1[i] * ee2;
            }
            return result;
        } else if (e1 instanceof Double && e2 instanceof double[]) {
            double ee1 = (Double) e1;
            double[] ee2 = (double[]) e2;
            double[] result = new double[ee2.length];
            for (int i = 0; i < result.length; i++) {
                result[i] = ee1 * ee2[i];
            }
            return result;
        } else if (e1 instanceof double[] && e2 instanceof double[]) {
            double[] ee1 = (double[]) e1;
            double[] ee2 = (double[]) e2;
            if (ee1.length != ee2.length) {
                throw new Exception("Different array length");
            }
            double[] result = new double[ee1.length];
            for (int i = 0; i < result.length; i++) {
                result[i] = ee1[i] * ee2[i];
            }
            return result;
        } else {
            throw new Exception(
                "Wrong operator usage: " +
                e1.getClass().getName() + " * " + e2.getClass().getName(), null);
        }
    }


    public static final Object DIV(Object e1, Object e2) throws Exception {
        if (e1 instanceof Double && e2 instanceof Double) {
            return (Double) e1 / (Double) e2;
        } else if (e1 instanceof double[] && e2 instanceof Double) {
            double[] ee1 = (double[]) e1;
            double ee2 = (Double) e2;
            double[] result = new double[ee1.length];
            for (int i = 0; i < result.length; i++) {
                result[i] = ee1[i] / ee2;
            }
            return result;
        } else if (e1 instanceof Double && e2 instanceof double[]) {
            double ee1 = (Double) e1;
            double[] ee2 = (double[]) e2;
            double[] result = new double[ee2.length];
            for (int i = 0; i < result.length; i++) {
                result[i] = ee1 / ee2[i];
            }
            return result;
        } else if (e1 instanceof double[] && e2 instanceof double[]) {
            double[] ee1 = (double[]) e1;
            double[] ee2 = (double[]) e2;
            if (ee1.length != ee2.length) {
                throw new Exception("Different array length");
            }
            double[] result = new double[ee1.length];
            for (int i = 0; i < result.length; i++) {
                result[i] = ee1[i] / ee2[i];
            }
            return result;
        } else {
            throw new Exception(
                "Wrong operator usage: " +
                e1.getClass().getName() + " / " + e2.getClass().getName(), null);
        }
    }


    public static final Object MOD(Object e1, Object e2) throws Exception {
        if (e1 instanceof Double && e2 instanceof Double) {
            return (Double) e1 % (Double) e2;
        } else if (e1 instanceof double[] && e2 instanceof Double) {
            double[] ee1 = (double[]) e1;
            double ee2 = (Double) e2;
            double[] result = new double[ee1.length];
            for (int i = 0; i < result.length; i++) {
                result[i] = ee1[i] % ee2;
            }
            return result;
        } else if (e1 instanceof Double && e2 instanceof double[]) {
            double ee1 = (Double) e1;
            double[] ee2 = (double[]) e2;
            double[] result = new double[ee2.length];
            for (int i = 0; i < result.length; i++) {
                result[i] = ee1 % ee2[i];
            }
            return result;
        } else if (e1 instanceof double[] && e2 instanceof double[]) {
            double[] ee1 = (double[]) e1;
            double[] ee2 = (double[]) e2;
            if (ee1.length != ee2.length) {
                throw new Exception("Different array length");
            }
            double[] result = new double[ee1.length];
            for (int i = 0; i < result.length; i++) {
                result[i] = ee1[i] % ee2[i];
            }
            return result;
        } else {
            throw new Exception(
                "Wrong operator usage: " +
                e1.getClass().getName() + " % " + e2.getClass().getName(), null);
        }
    }


    public static final Object POW(Object e1, Object e2) throws Exception {
        if (e1 instanceof Double && e2 instanceof Double) {
            return Math.pow( (Double) e1, (Double) e2);
        } else if (e1 instanceof double[] && e2 instanceof Double) {
            double[] ee1 = (double[]) e1;
            double ee2 = (Double) e2;
            double[] result = new double[ee1.length];
            for (int i = 0; i < result.length; i++) {
                result[i] = Math.pow(ee1[i], ee2);
            }
            return result;
        } else if (e1 instanceof Double && e2 instanceof double[]) {
            double ee1 = (Double) e1;
            double[] ee2 = (double[]) e2;
            double[] result = new double[ee2.length];
            for (int i = 0; i < result.length; i++) {
                result[i] = Math.pow(ee1, ee2[i]);
            }
            return result;
        } else if (e1 instanceof double[] && e2 instanceof double[]) {
            double[] ee1 = (double[]) e1;
            double[] ee2 = (double[]) e2;
            if (ee1.length != ee2.length) {
                throw new Exception("Different array length");
            }
            double[] result = new double[ee1.length];
            for (int i = 0; i < result.length; i++) {
                result[i] = Math.pow(ee1[i], ee2[i]);
            }
            return result;
        } else {
            throw new Exception(
                "Wrong operator usage: " +
                e1.getClass().getName() + " ^ " + e2.getClass().getName(), null);
        }
    }


    public static final Object IDENTIFIER(String id) throws Exception {
        if (!parser.variables.containsKey(id)) {
            throw new Exception("Variable " + id + " does not exist.");
        }

        Object result = parser.variables.get(id);
        if (result == null) {
            throw new NullPointerException("Variable " + id + " is null.");
        }

        return result;
    }


    public static final Object ARRAY_ACCESS(Object _array, Object _i) throws Exception {
        if (_array instanceof double[] && _i instanceof Double) {
            double[] array = (double[]) _array;
            int i = (int) ( (Double) _i).doubleValue();

            if (i < 0 || i >= array.length) {
                throw new Exception("Wrong array[i] usage: i == " + i);
            }

            return array[i];
        } else {
            throw new Exception("Wrong array[i] usage: object types do not match, " +
                                _array.getClass().getName() + ", " +
                                _i.getClass().getName() + ", "
                );
        }
    }


    public static final Object ARRAY_ACCESS(Object _array, Object _i1, Object _i2) throws
        Exception {

        if (_array instanceof double[] && _i1 instanceof Double && _i2 instanceof Double) {
            double[] array = (double[]) _array;
            int i1 = (int) ( (Double) _i1).doubleValue();
            int i2 = (int) ( (Double) _i2).doubleValue();

            if (i1 < 0 || i1 >= array.length) {
                throw new Exception("Wrong array[i1, i2] usage: i1 == " + i1);
            }

            if (i2 < 0 || i2 >= array.length) {
                throw new Exception("Wrong array[i1, i2] usage: i2 = " + i2);
            }

            double[] newArray = new double[i2 - i1 + 1];
            System.arraycopy(array, i1, newArray, 0, newArray.length);
            return newArray;
        } else {
            throw new Exception("Wrong array[i1, i2] usage: object types do not match, " +
                                _array.getClass().getName() + ", " +
                                _i1.getClass().getName() + ", " +
                                _i2.getClass().getName() + ", "
                );
        }
    }


    public static final Object SUM(Object obj, Object index1, Object index2) throws
        Exception {

        if (obj instanceof double[]) {
            double[] a = (double[]) obj;

            if (! (index1 instanceof Double && index2 instanceof Double)) {
                throw new Exception("Wrong Array Indices");
            }

            int i1 = ( (Double) index1).intValue();
            int i2 = ( (Double) index2).intValue();
            double sum = 0;
            for (int i = i1; i <= i2; i++) {
                sum += a[i];
            }

            return new Double(sum);
        } else {
            throw new Exception("Improper use of sum(). An array sould be passed");
        }
    }


    public static final Object SUM(Object expr) throws
        Exception {

        if (expr instanceof Double) {
            return expr;
        } else if (expr instanceof double[]) {
            double[] a = (double[]) expr;

            double sum = 0;
            for (int i = 0; i < a.length; i++) {
                sum += a[i];
            }

            return new Double(sum);
        } else {
            throw new Exception("Improper use of sum()");
        }
    }


    public static final Object MEAN(Object obj, Object index1, Object index2) throws
        Exception {

        if (obj instanceof double[]) {
            double[] a = (double[]) obj;

            int i1, i2;
            if (index1 == null && index2 == null) {
                i1 = 0;
                i2 = a.length - 1;
            } else if (index1 instanceof Double && index2 instanceof Double) {
                i1 = ( (Double) index1).intValue();
                i2 = ( (Double) index2).intValue();
            } else {
                throw new Exception("Wrong Array Indices");
            }

            double mean = 0;
            for (int i = i1; i <= i2; i++) {
                mean += a[i];
            }
            mean /= (i2 - i1 + 1);

            return new Double(mean);
        } else {
            throw new Exception("Improper use of mean(). An array sould be passed");
        }
    }


    public static final Object RAD_NORM(Object o1, Object o2) throws Exception {
        if (o1 instanceof Double && o2 instanceof Double) {
            double angle = (Double) o1;
            double period = (Double) o2;

            while (angle < 0) {
                angle += period;
            }

            while (angle > period) {
                angle -= period;
            }

            return new Double(angle);
        } else {
            throw new Exception("Wrong function usage");
        }
    }


    public static final Object EQUAL(Object e1, Object e2) throws Exception {
        if (e1 instanceof Boolean && e2 instanceof Boolean) {
            return new Boolean(
                ( (Boolean) e1).booleanValue() == ( (Boolean) e2).booleanValue());
        } else if (e1 instanceof Double && e2 instanceof Double) {
            return new Boolean(
                ( (Double) e1).doubleValue() == ( (Double) e2).doubleValue());
        } else if (e1 instanceof String && e2 instanceof String) {
            return new Boolean(e1.equals(e2));
        } else {
            throw new Exception(
                "Wrong operator usage: " +
                e1.getClass().getName() + " == " + e2.getClass().getName(), null);
        }
    }


    public static final Object NOT_EQUAL(Object e1, Object e2) throws Exception {
        if (e1 instanceof Boolean && e2 instanceof Boolean) {
            return new Boolean( ( (Boolean) e1).booleanValue() !=
                               ( (Boolean) e2).booleanValue());
        } else if (e1 instanceof Double && e2 instanceof Double) {
            return new Boolean( ( (Double) e1).doubleValue() !=
                               ( (Double) e2).doubleValue());
        } else if (e1 instanceof String && e2 instanceof String) {
            return new Boolean(!e1.equals(e2));
        } else {
            throw new Exception(
                "Wrong operator usage: " +
                e1.getClass().getName() + " != " + e2.getClass().getName(), null);
        }
    }


    public static final Object LESS(Object e1, Object e2) throws Exception {
        if (e1 instanceof Double && e2 instanceof Double) {
            return (Double) e1 < (Double) e2;
        } else {
            throw new Exception(
                "Wrong operator usage: " +
                e1.getClass().getName() + " < " + e2.getClass().getName(), null);
        }
    }


    public static final Object LESS_EQUAL(Object e1, Object e2) throws Exception {
        if (e1 instanceof Double && e2 instanceof Double) {
            return (Double) e1 <= (Double) e2;
        } else {
            throw new Exception(
                "Wrong operator usage: " +
                e1.getClass().getName() + " <= " + e2.getClass().getName(), null);
        }
    }


    public static final Object GREATER(Object e1, Object e2) throws Exception {
        if (e1 instanceof Double && e2 instanceof Double) {
            return (Double) e1 > (Double) e2;
        } else {
            throw new Exception(
                "Wrong operator usage: " +
                e1.getClass().getName() + " > " + e2.getClass().getName(), null);
        }
    }


    public static final Object GREATER_EQUAL(Object e1, Object e2) throws Exception {
        if (e1 instanceof Double && e2 instanceof Double) {
            return (Double) e1 >= (Double) e2;
        } else {
            throw new Exception(
                "Wrong operator usage: " +
                e1.getClass().getName() + " >= " + e2.getClass().getName(), null);
        }
    }


    public static final Object NOT(Object e) throws Exception {
        if (e instanceof Boolean) {
            return! (Boolean) e;
        } else {
            throw new Exception(
                "Wrong operator usage: " + "! " + e.getClass().getName(), null);
        }
    }


    public static final Object LOG10(Object e) throws Exception {
        if (e instanceof double[]) {
            double[] a = (double[]) e;
            double[] b = new double[a.length];
            for (int i = 0; i < a.length; i++) {
                b[i] = Math.log(a[i]) / Math.log(10);
            }
            return b;
        } else {
            throw new Exception(
                "Wrong operator usage: " + "! " + e.getClass().getName(), null);
        }
    }


    public static final Object ABS(Object e) throws Exception {
        if (e instanceof Double) {
            return Math.abs( (Double) e);
        } else if (e instanceof double[]) {
            double[] array = (double[]) e;
            for (int i = 0; i < array.length; i++) {
                array[i] = Math.abs(array[i]);
            }
            return array;
        } else {
            throw new Exception("Wrong function usage: exp");
        }
    }


    public static final Object EXP(Object e) throws Exception {
        if (e instanceof Double) {
            return Math.exp( (Double) e);
        } else {
            throw new Exception("Wrong function usage: exp");
        }
    }


    public static final Object AND(Object e1, Object e2) throws Exception {
        if (e1 instanceof Boolean && e2 instanceof Boolean) {
            return (Boolean) e1 && (Boolean) e2;
        } else {
            throw new Exception(
                "Wrong operator usage: " +
                e1.getClass().getName() + " & " + e2.getClass().getName(), null);
        }
    }


    public static final Object OR(Object e1, Object e2) throws Exception {
        if (e1 instanceof Boolean && e2 instanceof Boolean) {
            return (Boolean) e1 || (Boolean) e2;
        } else {
            throw new Exception(
                "Wrong operator usage: " +
                e1.getClass().getName() + " | " + e2.getClass().getName(), null);
        }
    }


    public static final Object UNARY_MINUS(Object e) throws Exception {
        if (e instanceof Double) {
            return - (Double) e;
        }
        if (e instanceof double[]) {
            double[] ee = (double[]) e;
            double[] result = new double[ee.length];
            for (int i = 0; i < ee.length; i++) {
                result[i] = ee[i];
            }
            return result;
        } else {
            throw new Exception(
                "Wrong operator usage: " + "- " + e.getClass().getName(), null);
        }
    }


    public static final Object MIN(Object e1, Object e2) throws Exception {
        if (e1 instanceof Double && e2 instanceof Double) {
            return Math.min( (Double) e1, (Double) e2);
        } else {
            throw new Exception(
                "Wrong function usage: " +
                "min(" + e1.getClass().getName() + ", " + e2.getClass().getName() + ")", null);
        }
    }


    public static final Object MAX(Object e1, Object e2) throws Exception {
        if (e1 instanceof Double && e2 instanceof Double) {
            return Math.max( (Double) e1, (Double) e2);
        } else {
            throw new Exception(
                "Wrong function usage: " +
                "max(" + e1.getClass().getName() + ", " + e2.getClass().getName() + ")", null);
        }
    }


    public static final Object MAX(Object e) throws Exception {
        if (e instanceof double[]) {
            double[] a = (double[]) e;
            double max = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < a.length; i++) {
                if (a[i] > max) {
                    max = a[i];
                }
            }
            return new Double(max);
        } else {
            throw new Exception("Wrong function usage: MAX");
        }
    }


    public static final Object MAX(Object e, Object i1, Object i2) throws Exception {
        if (e instanceof double[] && i1 instanceof Double && i2 instanceof Double) {
            double[] a = (double[]) e;
            double max = Double.NEGATIVE_INFINITY;

            for (int i = ( (Double) i1).intValue(); i <= ( (Double) i2).intValue(); i++) {
                if (a[i] > max) {
                    max = a[i];
                }
            }
            return max;
        } else {
            throw new Exception("Wrong function usage: MAX");
        }
    }


    public static final Object EXTREME(Object e) throws Exception {
        if (e instanceof double[]) {
            double[] a = (double[]) e;
            double max = Double.NEGATIVE_INFINITY;
            double min = Double.POSITIVE_INFINITY;
            for (int i = 0; i < a.length; i++) {
                if (a[i] > max) {
                    max = a[i];
                }
                if (a[i] < min) {
                    min = a[i];
                }
            }
            return (Math.abs(max) > Math.abs(min)) ? new Double(max) : new Double(min);
        } else {
            throw new Exception("Wrong function usage: EXTREME");
        }
    }


    public static final Object D(Object e) throws Exception {
        if (e instanceof double[]) {
            double[] a = (double[]) e;
            double[] b = new double[a.length - 1];

            for (int i = 0; i < a.length - 1; i++) {
                b[i] = a[i + 1] - a[i];
            }

            return b;
        } else {
            throw new Exception("Wrong function usage: D");
        }
    }
}
