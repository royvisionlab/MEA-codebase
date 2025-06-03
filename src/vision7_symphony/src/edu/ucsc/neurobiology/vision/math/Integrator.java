package edu.ucsc.neurobiology.vision.math;

import java.util.*;

import edu.ucsc.neurobiology.vision.plot.*;


/**
 * Various integration functions.
 * These funtions do not work with rapidly variing functions.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class Integrator {
    // variables used by gauss8
    static final double x1 = 1.83434642495649805E-01;
    static final double x2 = 5.25532409916328986E-01;
    static final double x3 = 7.96666477413626740E-01;
    static final double x4 = 9.60289856497536232E-01;
    static final double w1 = 3.62683783378361983E-01;
    static final double w2 = 3.13706645877887287E-01;
    static final double w3 = 2.22381034453374471E-01;
    static final double w4 = 1.01228536290376259E-01;
    static final double sq2 = 1.41421356E0;
    static final int nlmn = 1;
    static final int kmx = 5000;
    static final int kml = 6;
    static double aa[] = new double[61];
    static double hh[] = new double[61];
    static double vl[] = new double[61];
    static double gr[] = new double[61];
    static int lr[] = new int[61];


    /**
     * This method integrates real functions of one variable over finite
     * intervals using an adaptive 8-point Legendre-Gauss algorithm.
     * Gaus8 is intended primarily for high accuracy integration
     * or integration of smooth functions.<p>
     * This method is a Java translation of the FORTRAN
     * routine gaus8 written by R. E. Jones.  The version translated
     * is part of the SLATEC library of numerical routines, and can be
     * downloaded at www.netlib.org.<p>
     * The maximum number of significant digits obtainable in the integral
     * is 16. (This is based on the translator's reading of the dgaus8.f
     * documentation and the translator's understanding of Java numerics.
     * Since the translator is a mathematical statistician rather than
     * a numerical analyst, the 16 figure might not be quite correct.)<p>
     * The translation was performed by Steve Verrill during
     * February of 2002.<p>
     *
     * @param   f   A class that defines a method, getValueAt,
     *                      that returns a value that is to be integrated.
     *                      The class must implement the FunctionData
     *                      interface (see the definition in
     *                      FunctionData.java).
     *                      See Gaus8Test.java for an example of such
     *                      a class.  getValueAt must have one double
     *                      valued argument.
     * @param   a            The lower limit of the integral.
     * @param   b            The upper limit of the integral (may be less
     *                      than a).
     * @param   err          On input, err[0] is a requested pseudorelative
     *                      error tolerance.  Normally pick a value
     *                      of abs(err[0]) so that dtol < abs(err[0]) <= .001
     *                      where dtol is the double precision unit
     *                      roundoff (2.22e-16 for IEEE 754 numbers
     *                      according to the SLATEC documentation).  The
     *                      answer returned by gaus8 will normally have
     *                      no more error than Math.abs(err[0]) times the
     *                      integral of the absolute value of
     *                      getValueAt.  Usually, smaller values for
     *                      err[0] yield more accuracy and require more
     *                      function evaluations.
     *                      A negative value for err[0] causes an estimate
     *                      of the absolute error in the integral to be
     *                      returned in err[0].
     */
    synchronized public static double gauss8(
        FunctionData f, double a, double b, double err) throws CannotEvaluateException {

        int k, l, lmn, lmx, mxl, nbits,
            nib, nlmx;

        double ae, anib, area, c, ce, ee, ef, eps, est,
            gl, glr, tol, vr, x, h, g8xh, ans;

        for (int i = 0; i < 61; i++) {
            aa[i] = 0;
            hh[i] = 0;
            vl[i] = 0;
            gr[i] = 0;
            lr[i] = 0;
        }

        nbits = 53;
        nlmx = Math.min(60, (nbits * 5) / 8);
        ans = 0.0;
        ce = 0.0;

        if (a == b) {
            if (err < 0.0) {
                err = ce;
            }
            return ans;
        }

        lmx = nlmx;
        lmn = nlmn;

        if (b != 0.0) {
            if (sign(1.0, b) * a > 0.0) {
                c = Math.abs(1.0 - a / b);
                if (c <= 0.1) {
                    if (c <= 0.0) {
                        if (err < 0.0) {
                            err = ce;
                        }
                        return ans;
                    }

                    anib = 0.5 - Math.log(c) / 0.69314718E0;
                    nib = (int) (anib);
                    lmx = Math.min(nlmx, nbits - nib - 7);
                    if (lmx < 1) {
                        throw new IllegalArgumentException(
                            "\n\nGaus8 --- a and b are too nearly"
                            + " equal to allow normal integration.\n"
                            + "ans is set to 0 and ierr[0] is set to -1.\n\n");
                    }
                    lmn = Math.min(lmn, lmx);
                }
            }
        }

        tol = Math.max(Math.abs(err), Math.pow(2.0, 5 - nbits)) / 2.0;
        // According to the SLATEC documentation 2.22e-16 is the
        // appropriate value for IEEE 754 double precision unit roundoff.
        if (err == 0.0) {
            tol = Math.sqrt(2.22e-16);
        }
        eps = tol;
        hh[1] = (b - a) / 4.0;
        aa[1] = a;
        lr[1] = 1;
        l = 1;
        x = aa[l] + 2.0 * hh[l];
        h = 2.0 * hh[l];
        g8xh = h * (w1 * (f.getValueAt(x - x1 * h) + f.getValueAt(x + x1 * h)) +
                    w2 * (f.getValueAt(x - x2 * h) + f.getValueAt(x + x2 * h)) +
                    w3 * (f.getValueAt(x - x3 * h) + f.getValueAt(x + x3 * h)) +
                    w4 * (f.getValueAt(x - x4 * h) + f.getValueAt(x + x4 * h)));
        est = g8xh;
        k = 8;
        area = Math.abs(est);
        ef = 0.5;
        mxl = 0;

        // Compute refined estimates, estimate the error, etc.
        while (true) {
            x = aa[l] + hh[l];
            h = hh[l];
            g8xh = h * (w1 * (f.getValueAt(x - x1 * h) + f.getValueAt(x + x1 * h)) +
                        w2 * (f.getValueAt(x - x2 * h) + f.getValueAt(x + x2 * h)) +
                        w3 * (f.getValueAt(x - x3 * h) + f.getValueAt(x + x3 * h)) +
                        w4 * (f.getValueAt(x - x4 * h) + f.getValueAt(x + x4 * h)));
            gl = g8xh;
            x = aa[l] + 3.0 * hh[l];
            h = hh[l];
            g8xh = h * (w1 * (f.getValueAt(x - x1 * h) + f.getValueAt(x + x1 * h)) +
                        w2 * (f.getValueAt(x - x2 * h) + f.getValueAt(x + x2 * h)) +
                        w3 * (f.getValueAt(x - x3 * h) + f.getValueAt(x + x3 * h)) +
                        w4 * (f.getValueAt(x - x4 * h) + f.getValueAt(x + x4 * h)));
            gr[l] = g8xh;
            k += 16;
            area += (Math.abs(gl) + Math.abs(gr[l]) - Math.abs(est));
            glr = gl + gr[l];
            ee = Math.abs(est - glr) * ef;
            ae = Math.max(eps * area, tol * Math.abs(glr));
            if (ee - ae > 0.0) {
                if (k > kmx) {
                    lmx = kml;
                }
                if (l < lmx) {
                    l++;
                    eps *= 0.5;
                    ef /= sq2;
                    hh[l] = hh[l - 1] * 0.5;
                    lr[l] = -1;
                    aa[l] = aa[l - 1];
                    est = gl;
                } else {
                    mxl = 1;
                    ce += (est - glr);
                    if (lr[l] <= 0) {
                        vl[l] = glr;
                        est = gr[l - 1];
                        lr[l] = 1;
                        aa[l] += 4.0 * hh[l];
                    } else {
                        vr = glr;
                        while (true) {
                            if (l <= 1) {
                                ans = vr;
                                if ( (mxl != 0) && (Math.abs(ce) > 2.0 * tol * area)) {
                                    throw new IllegalArgumentException(
                                        "nans is probably insufficiently accurate.");
                                }
                                if (err < 0.0) {
                                    err = ce;
                                }
                                return ans;
                            }
                            l--;
                            eps *= 2.0;
                            ef *= sq2;
                            if (lr[l] <= 0.0) {
                                break;
                            }
                            vr += vl[l + 1];
                        }
                        vl[l] = vl[l + 1] + vr;
                        est = gr[l - 1];
                        lr[l] = 1;
                        aa[l] += 4.0 * hh[l];
                    }
                }
            } else {
                ce += (est - glr);

                if (lr[l] <= 0) {
                    vl[l] = glr;
                    est = gr[l - 1];
                    lr[l] = 1;
                    aa[l] += 4.0 * hh[l];

                } else {
                    vr = glr;
                    while (true) {
                        if (l <= 1) {
                            ans = vr;
                            if ( (mxl != 0) && (Math.abs(ce) > 2.0 * tol * area)) {
                                throw new IllegalArgumentException(
                                    "nans is probably insufficiently accurate");
                            }
                            if (err < 0.0) {
                                err = ce;
                            }
                            return ans;
                        }
                        l--;
                        eps *= 2.0;
                        ef *= sq2;
                        if (lr[l] <= 0.0) {
                            break;
                        }
                        vr += vl[l + 1];
                    }
                    vl[l] = vl[l + 1] + vr;
                    est = gr[l - 1];
                    lr[l] = 1;
                    aa[l] += 4.0 * hh[l];
                }
            }
        }
    }


    private final static double sign(double a, double b) {
        if (b < 0.0) {
            return -Math.abs(a);
        } else {
            return Math.abs(a);
        }
    }


    private static double trapzd(FunctionData f, double a, double b, int n, double s) throws
        CannotEvaluateException {

        double x, tnm, sum, del;
        int it, j;

        if (n == 1) {
            return (s = 0.5 * (b - a) * (f.getValueAt(a) + f.getValueAt(b)));
        } else {
            for (it = 1, j = 1; j < n - 1; j++) {
                it <<= 1;
            }
            tnm = it;
            del = (b - a) / tnm;
            x = a + 0.5 * del;
            for (sum = 0.0, j = 1; j <= it; j++, x += del) {
                sum += f.getValueAt(x);
            }
            s = 0.5 * (s + (b - a) * sum / tnm);

            return s;
        }
    }


    public static double integrateSimple(FunctionData f, double a, double b, int nSteps) throws
        CannotEvaluateException {

        double dx = (b - a) / nSteps;
//        double s = f.getValueAt(a) + f.getValueAt(b);\
        double s = 0;
        double x = a + dx / 2;

        for (int j = 0; j < nSteps; j++, x += dx) {
            s += f.getValueAt(x);
        }

        return s * dx;
    }


    static int JMAX = 20;


    public static double qtrap(FunctionData f, double a, double b, double eps) throws
        CannotEvaluateException {
        int j;
        double s = 0, olds;

        olds = -1.0e30;
        for (j = 0; j < JMAX; j++) {
            s = trapzd(f, a, b, j, s);
            if (Math.abs(s - olds) <= eps * Math.abs(olds)) {
                return s;
            }
            olds = s;
        }

        throw new TooManyIterationsException("Too many steps in routine qtrap");
    }


    public static double qsimp(FunctionData f, double a, double b, double eps) throws
        TooManyIterationsException, CannotEvaluateException {

        int j;
        double s, st = 0, ost, os;

        ost = os = -1.0e30;
        for (j = 0; j < JMAX; j++) {
            st = trapzd(f, a, b, j, st);
            s = (4.0 * st - ost) / 3.0;
            if (Math.abs(s - os) <= eps * Math.abs(os)) {
                return s;
            }
            os = s;
            ost = st;
        }

        throw new TooManyIterationsException("a = " + a + ", b = " + b);
    }



    private static FunctionData gauss = new FunctionDataAdapter() {
        public double getValueAt(double x) {
            return Math.exp( -x * x);
        }
    };


//    public static double erf(double x1, double x2) throws CannotEvaluateException {
//        return gauss8(gauss, x1, x2, 1e-6);
//    }
//

    /*It is a bit more accurate to call this min(abs(f))*/
    public static double findZero(FunctionData f, double x1, double x2, double dx) throws
        CannotEvaluateException {

        double yMin = Double.POSITIVE_INFINITY;
        double xMin = Double.NaN;
        for (double x = x1; x < x2; x += dx) {
            double y = Math.abs(f.getValueAt(x));
            if (y < yMin) {
                yMin = y;
                xMin = x;
            }
        }
        return xMin;
    }


    /*written by Matthew Grivich, returns all zeros between x1 and x2 on f to
     the accuracy of dx*/
    public static double[] findZeros(FunctionData f, double x1, double x2, double dx) throws
        CannotEvaluateException {
        ArrayList<Double> zerosList = new ArrayList<Double> ();
        double oldY = 0.0;
        for (double x = x1; x < x2; x += dx) {
            double y = f.getValueAt(x);
            //if y equals zero, zero found
            if (y == 0.0) {
                zerosList.add(new Double(x));
            }

            //if y and old y have a different sign, zero found
            if (y * oldY < 0.0) {
                zerosList.add(new Double(x - dx / 2));
            }
            oldY = y;
        }

        double[] zeros = new double[zerosList.size()];
        for (int i = 0; i < zeros.length; i++) {
            zeros[i] = ( (Double) zerosList.get(i)).doubleValue();
        }
        return zeros;
    }
    
 
}
