package edu.ucsc.neurobiology.vision.math;

import edu.ucsc.neurobiology.vision.plot.FunctionData;


/**
 * A spline that can be used only if the points on the x axis are evenly spaced.
 * Taken from Numerical Recipes in Fortran: 3.3 Cubic Spline Interpolation.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class UniformSpline implements FunctionData {

    private final double c = 1.0 / 6.0;
    private final int n;
    private float[] y;
    private final double[] derivativesTemplate, secondDerivatives, u;


    public UniformSpline(int n) {
        this.n = n;
        y = new float[n];
        u = new double[n];
        secondDerivatives = new double[n];
        
        derivativesTemplate = new double[n-1];
        derivativesTemplate[0] = 0;
        for (int i = 1; i < n-1; i++)
            derivativesTemplate[i] = -1 / (derivativesTemplate[i-1] + 4);
    }


    public void reSpline(double yArray[]) {
        if (yArray.length != n)
            throw new IllegalArgumentException("y.length != n");

        for (int i = 0; i < n; i++)
            y[i] = (float) yArray[i];
        spline();
    }


    public void reSpline(float yArray[]) {
        if (yArray.length != n)
            throw new IllegalArgumentException("y.length != n");

        System.arraycopy(yArray, 0, y, 0, n);
        spline();
    }


    public void reSpline(short yArray[]) {
        if (yArray.length != n)
            throw new IllegalArgumentException("y.length != n");
        
        for (int i = 0; i < n; i++)
            y[i] = yArray[i];
        spline();
    }


    public String getDescription() {
        return "spline";
    }

    
    // Old version, optimized below
    private void splineOld() {
        double p, qn, un;

        for (int i = 0; i < n; i++) {
            secondDerivatives[i] = 0;
            u[i] = 0;
        }

        for (int i = 1; i < n-1; i++) {
            p = 0.5 * secondDerivatives[i-1] + 2;
            secondDerivatives[i] = -0.5 / p;
            u[i] = (y[i+1] - y[i]) - (y[i] - y[i-1]);
            u[i] = (3 * u[i] - 0.5 * u[i-1]) / p;
        }

        // To match NR online, the below should read
        // 	 (un - qn * u[n-2]) / (qn * secondDerivatives[n-2] + 1)
        // But of course with qn = un = 0 this is irrelevant
        qn = un = 0;
        secondDerivatives[n-1] = (un - qn * u[n-1]) / (qn * secondDerivatives[n-1] + 1);

        for (int k = n-2; k >= 0; k--) {
            secondDerivatives[k] = secondDerivatives[k] * secondDerivatives[k+1] + u[k];
        }

        // added
        for (int k = 0; k < n; k++) {
            secondDerivatives[k] *= c;
        }
    }
    
    
    /**
     * Optimized version, somewhat faster than original above
     * Floating point rounding differences due to the removal of p, means this doesn't come out identical
     * to above, but it is practically identical.
     * @author Peter H. Li, The Salk Institute
     */
    private void spline() {
        u[0] = 0;
        for (int i = 1; i < n-1; i++) {
            u[i] = (y[i+1] - y[i]) - (y[i] - y[i-1]);
            u[i] = (u[i-1] - 6 * u[i]) * derivativesTemplate[i];
        }
        
        secondDerivatives[n-1] = 0;
        for (int k = n-2; k >= 0; k--)
            secondDerivatives[k] = derivativesTemplate[k] * secondDerivatives[k+1] + u[k];

        for (int k = 0; k < n; k++)
            secondDerivatives[k] *= c;
    }


    public final double getValueAt(final double coord) {
        final int j1 = (int) coord;
        final int j2 = j1 + 1;
        final double a = j2 - coord;
        final double b = coord - j1;

        return
            a * (y[j1] + (a*a - 1) * secondDerivatives[j1]) +
            b * (y[j2] + (b*b - 1) * secondDerivatives[j2]);
    }

}
