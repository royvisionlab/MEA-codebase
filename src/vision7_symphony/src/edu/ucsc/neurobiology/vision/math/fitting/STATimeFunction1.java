package edu.ucsc.neurobiology.vision.math.fitting;

import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.plot.*;


/**
 * Implements a time filter fitting function.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 * Generalized to 1-3 filter by Matthew Grivich
 */
public class STATimeFunction1
    extends FittableFunction implements FunctionData, FunctionError, PlotData {

    public final double n1, n2, n3;
    public double offset;
    private double norm = 1.0;
    


    /**
     * One filter constructor
     *
     * @param a1 double
     * @param t1 double
     * @param n1 double
     */
    public STATimeFunction1(double a1, double t1, double n1) {
        setParameters(new double[] {a1, t1, 0.0, 1.0, 0.0, 1.0});
        adjustParameters = new boolean[] {true, true, false, false, false, false};
        this.n1 = n1;
        this.n2 = 0.0;
        this.n3 = 0.0;
        this.offset = 0.0;
    }
    public STATimeFunction1(
        double a1, double t1, double n1, double offset) {

        setParameters(new double[] {a1, t1, 0.0, 1.0, 0.0, 1.0});
        adjustParameters = new boolean[] {true, true, false, false, false, false};
        this.n1 = n1;
        this.n2 = 0.0;
        this.n3 = 0.0;
        this.offset = offset;
    }


    /**
     * Two filter constructor
     *
     * @param a1 double
     * @param t1 double
     * @param n1 double
     */
    
    public STATimeFunction1(
            double a1, double t1, double a2, double t2, double n1, double n2) {
        setParameters(new double[] {a1, t1, a2, t2, 0.0, 1.0});
        adjustParameters = new boolean[] {true, true, true, true, false, false};
        this.n1 = n1;
        this.n2 = n2;
        this.n3 = 0.0;
        this.offset = 0.0;
        
    }
    public STATimeFunction1(
        double a1, double t1, double a2, double t2, double n1, double n2, double offset) {

        setParameters(new double[] {a1, t1, a2, t2, 0.0, 1.0});
        adjustParameters = new boolean[] {true, true, true, true, false, false};
        this.n1 = n1;
        this.n2 = n2;
        this.n3 = 0.0;
        this.offset = offset;
    }


    /**
     * Three filter constructor
     *
     * @param a1 double
     * @param t1 double
     * @param n1 double
     */
    
    public STATimeFunction1(
            double a1, double t1, double a2, double t2, double a3, double t3, double n1,
            double n2, double n3) {

            setParameters(new double[] {a1, t1, a2, t2, a3, t3});
            adjustParameters = new boolean[] {true, true, true, true, true, true};
            this.n1 = n1;
            this.n2 = n2;
            this.n3 = n3;
            this.offset = 0.0;
        }
    
    public STATimeFunction1(
        double a1, double t1, double a2, double t2, double a3, double t3, double n1,
        double n2, double n3, double offset) {

        setParameters(new double[] {a1, t1, a2, t2, a3, t3});
        adjustParameters = new boolean[] {true, true, true, true, true, true};
        this.n1 = n1;
        this.n2 = n2;
        this.n3 = n3;
        this.offset = offset;
    }


    public Num getValueAndErrorAt(double t) throws CannotEvaluateException {
        Num a1 = param(0);
        Num t1 = param(1);
        Num a2 = param(2);
        Num t2 = param(3);
        Num a3 = param(4);
        Num t3 = param(5);

        Num t_div_t1 = new Num(t, 0).div(t1);
        Num t_div_t2 = new Num(t, 0).div(t2);
        Num t_div_t3 = new Num(t, 0).div(t3);
        Num p1 = t_div_t1.mul(Num.one.sub(t_div_t1).exp()).pow(new Num(n1, 0));
        Num p2 = t_div_t2.mul(Num.one.sub(t_div_t2).exp()).pow(new Num(n2, 0));
        Num p3 = t_div_t3.mul(Num.one.sub(t_div_t3).exp()).pow(new Num(n3, 0));
        Num temp = a1.mul(p1).add(a2.mul(p2).add(a3.mul(p3)));
        temp = temp.div(norm);
        return temp.add(new Num(offset, 0));
    }


    public double getValueAt(double t) {
        double a1 = parameters[0];
        double t1 = parameters[1];
        double a2 = parameters[2];
        double t2 = parameters[3];
        double a3 = parameters[4];
        double t3 = parameters[5];

        double p1 = Math.pow( (t / t1) * Math.exp(1 - t / t1), n1);
        double p2 = Math.pow( (t / t2) * Math.exp(1 - t / t2), n2);
        double p3 = Math.pow( (t / t3) * Math.exp(1 - t / t3), n3);

        return (a1 * p1 + a2 * p2 + a3 * p3)/norm+offset;
    }


    public double pow(double x, double n) {
        if (n - (int) n == 0 && n < 50) {
            double v = 1;
            for (int i = 0; i < n; i++) {
                v *= x;
            }
            return v;
        } else {
            return Math.pow(x, n);
        }
    }


    public double getValueAndDerivatives(
        final double[] x, final double[] p, final double[] derivatives) {

        double a1 = p[0];
        double t1 = p[1];
        double a2 = p[2];
        double t2 = p[3];
        double a3 = p[4];
        double t3 = p[5];
        double t = x[0];

        double p1 = pow( (t / t1) * Math.exp(1 - t / t1), n1);
        double p2 = pow( (t / t2) * Math.exp(1 - t / t2), n2);
        double p3 = pow( (t / t3) * Math.exp(1 - t / t3), n3);

        derivatives[0] = p1/norm;
        derivatives[1] = a1 * n1 * p1 * (t / t1 - 1) / t1/norm;
        derivatives[2] = p2/norm;
        derivatives[3] = a2 * n2 * p2 * (t / t2 - 1) / t2/norm;
        derivatives[4] = p3/norm;
        derivatives[5] = a3 * n3 * p3 * (t / t3 - 1) / t3/norm;
        return (a1 * p1 + a2 * p2 + a3 * p3)/norm+offset;
    }


    public double getA1() {
        return parameters[0];
    }


    public double getA2() {
        return parameters[2];
    }


    public double getA3() {
        return parameters[4];
    }


    public double getT1() {
        return parameters[1];
    }


    public double getT2() {
        return parameters[3];
    }


    public double getT3() {
        return parameters[5];
    }


    public double getN1() {
        return n1;
    }


    public double getN2() {
        return n2;
    }


    public double getN3() {
        return n3;
    }
    
    public double getOffset() {
        return offset;
    }


    public String getDescription() {
        return "STATimeFunction";
    }
    
    
    /**
     * Used to allow a timecourse function and a data timecourse to normalize together.
     * 
     * 
     */
 
    public void setNormalization(double normalization) {
        this.norm = normalization;
    }
    
    public double getNormalization() {
        return norm;
    }

}
