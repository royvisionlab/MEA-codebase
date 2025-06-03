package edu.ucsc.neurobiology.vision.test;

import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.plot.*;


/**
 *
 */
public class LowpassFilter
    extends FunctionDataAdapter {
    int n;
    double timeConstant;
    double c;


    public LowpassFilter(int n, int timeConstant) {
        this.n = n;
        c = 1 / (MathUtil.fact(n - 1) * Math.pow(timeConstant, n));
        this.timeConstant = timeConstant;
    }


    public final double getValueAt(double t) {
        return c * Math.pow(t, n - 1) * Math.exp( -t / timeConstant);
    }
}
