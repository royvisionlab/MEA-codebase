package edu.ucsc.neurobiology.vision.math;

import static edu.ucsc.neurobiology.vision.math.MathUtil.*;
import edu.ucsc.neurobiology.vision.math.fitting.*;
import edu.ucsc.neurobiology.vision.plot.*;


/**
 * A collectin of various minimization techniques.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class FunctionMinimization {
    private static final double GOLD = 1.61803399;
    private static final double R = GOLD - 1;
    private static final double C = 1 - R;
    private static final double ZEPS = 1e-10;
    private static final double TINY = 1e-20;
    private static final int GLIMIT = 100;
    private static final int ITMAX = 200;


    private static double SIGN(double a, double b) {
        return (b >= 0) ? fastAbs(a) : -fastAbs(a);
    }


    public static double goldenSection(
        FunctionData func, double ax, double bx, double precision) throws
        CannotEvaluateException {

        double f1, f2, x0, x1, x2, x3;
        double cx = (ax + bx) / 2;
        int iterations = 0;

        x0 = ax;
        x3 = cx;
        if (fastAbs(cx - bx) > fastAbs(bx - ax)) {
            x1 = bx;
            x2 = bx + C * (cx - bx);
        } else {
            x2 = bx;
            x1 = bx - C * (bx - ax);
        }

        f1 = func.getValueAt(x1);
        f2 = func.getValueAt(x2);

        while (fastAbs(x3 - x0) > precision * (fastAbs(x1) + fastAbs(x2))) {
            iterations++;

            if (f2 < f1) {
                x0 = x1;
                x1 = x2;
                x2 = R * x1 + C * x3;
                f1 = f2;
                f2 = func.getValueAt(x2);
            } else {
                x3 = x2;
                x2 = x1;
                x1 = R * x2 + C * x0;
                f2 = f1;
                f1 = func.getValueAt(x1);
            }
        }
        if (f1 < f2) {
            return x1;
        } else {
            return x2;
        }
    }


    static final class MinimumBracket {
        public double ax, bx, cx, fa, fb, fc;

        public void arrange() {
            if (ax > bx && bx > cx) {
                double t = 0;

                t = ax;
                ax = cx;
                cx = t;

                t = fa;
                fa = fc;
                fc = t;
            }
        }


        public String toString() {
            return "\n" +
                ax + " : " + fa + "\n" +
                bx + " : " + fb + "\n" +
                cx + " : " + fc + "\n";
        }
    }


    private static void mnbrak(MinimumBracket b, FunctionData func) throws
        CannotEvaluateException {
        double ulim, u, r, q, fu, dum;

        b.fa = func.getValueAt(b.ax);
        b.fb = func.getValueAt(b.bx);

        if (b.fb > b.fa) {
            dum = b.ax;
            b.ax = b.bx;
            b.bx = dum;

            dum = b.fb;
            b.fb = b.fa;
            b.fa = dum;
        }
        b.cx = b.bx + GOLD * (b.bx - b.ax);
        b.fc = func.getValueAt(b.cx);
        out:while (b.fb > b.fc) {
            r = (b.bx - b.ax) * (b.fb - b.fc);
            q = (b.bx - b.cx) * (b.fb - b.fa);
            u = b.bx - ( (b.bx - b.cx) * q - (b.bx - b.ax) * r) /
                (2.0 * SIGN(Math.max(fastAbs(q - r), TINY), q - r));
            ulim = b.bx + GLIMIT * (b.cx - b.bx);
            if ( (b.bx - u) * (u - b.cx) > 0) {
                fu = func.getValueAt(u);
                if (fu < b.fc) {
                    b.ax = b.bx;
                    b.bx = u;
                    b.fa = b.fb;
                    b.fb = fu;
                    break out;
                } else if (fu > b.fb) {
                    b.cx = u;
                    b.fc = fu;
                    break out;
                }
                u = b.cx + GOLD * (b.cx - b.bx);
                fu = func.getValueAt(u);
            } else if ( (b.cx - u) * (u - ulim) > 0) {
                fu = func.getValueAt(u);
                if (fu < b.fc) {
                    b.bx = b.cx;
                    b.cx = u;
                    u = b.cx + GOLD * (b.cx - b.bx);
                    b.fb = b.fc;
                    b.fc = fu;
                    fu = func.getValueAt(u);
                }
            } else if ( (u - ulim) * (ulim - b.cx) >= 0) {
                u = ulim;
                fu = func.getValueAt(u);
            } else {
                u = b.cx + GOLD * (b.cx - b.bx);
                fu = func.getValueAt(u);
            }
            b.ax = b.bx;
            b.bx = b.cx;
            b.cx = u;
            b.fa = b.fb;
            b.fb = b.fc;
            b.fc = fu;
        }

        b.arrange();
    }


    /**
     * Minimizes a one dimensional function by using the parabolic interpolation method.
     *
     * @param func the function the minimum of which has to be found
     * @param ax the low limit
     * @param bx the high limit
     * @param error the error
     * @return the position of the minimum
     */
    public static double[] _brentParabolic(
        FunctionData func, double ax, double bx, double error) throws
        CannotEvaluateException {

        double a, b, d = 0, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
        double e = 0, cx = (ax + bx) / 2;

        a = (ax < cx) ? ax : cx;
        b = (ax > cx) ? ax : cx;
        x = w = v = bx;
        fw = fv = fx = func.getValueAt(x);

        for (int iter = 1; iter <= 100; iter++) {
            xm = 0.5 * (a + b);
            tol1 = error * fastAbs(x) + ZEPS;
            tol2 = 2 * tol1;
            if (fastAbs(x - xm) <= (tol2 - 0.5 * (b - a))) {
                return new double[] {x, fx};
            }
            if (fastAbs(e) > tol1) {
                r = (x - w) * (fx - fv);
                q = (x - v) * (fx - fw);
                p = (x - v) * q - (x - w) * r;
                q = 2 * (q - r);
                if (q > 0) {
                    p = -p;
                }
                q = fastAbs(q);
                etemp = e;
                e = d;
                if (fastAbs(p) >= fastAbs(0.5 * q * etemp) || p <= q * (a - x) ||
                    p >= q * (b - x)) {
                    d = C * (e = (x >= xm ? a - x : b - x));
                } else {
                    d = p / q;
                    u = x + d;
                    if (u - a < tol2 || b - u < tol2) {
                        d = SIGN(tol1, xm - x);
                    }
                }
            } else {
                e = (x >= xm) ? a - x : b - x;
                d = C * e;
            }
            u = (fastAbs(d) >= tol1) ? x + d : x + SIGN(tol1, d);
            fu = func.getValueAt(u);
            if (fu <= fx) {
                if (u >= x) {
                    a = x;
                } else {
                    b = x;
                }
                v = w;
                w = x;
                x = u;
                fv = fw;
                fw = fx;
                fx = fu;
            } else {
                if (u < x) {
                    a = u;
                } else {
                    b = u;
                }
//                if (fu <= fw || w == x) {
                if (fu <= fw || w == x) {
                    v = w;
                    w = u;
                    fv = fw;
                    fw = fu;
                } else if (fu <= fv || v == x || v == w) {
                    v = u;
                    fv = fu;
                }
            }
        }

        System.out.println("Too many iterations in brent");
        return null;
    }


    public static double brentParabolic(
        FunctionData func, double ax, double bx, double tol) throws
        CannotEvaluateException {

        return _brentParabolic(func, ax, bx, tol)[0];
    }


    public interface IterationInfo {
        public void receiveInfo(String info);
    }


    public static double powell(double[] params, double[][] xi, double ftol,
                                Function func, IterationInfo info) throws
        CannotEvaluateException {

        int i, ibig, j, n = xi.length, iter;
        double del, fp, fptt, fret;
        double[] pt = new double[n];
        double[] ptt = new double[n];
        double[] xit = new double[n];

        fret = func.getValue(params);
        for (j = 0; j < n; j++) {
            pt[j] = params[j];
        }

        for (iter = 1; ; ++iter) {
            fp = fret;
            ibig = 0;
            del = 0;

            for (i = 0; i < n; i++) {
                for (j = 0; j < n; j++) {
                    xit[j] = xi[j][i];
                }
                fptt = fret;

                fret = linmin(params, xit, func);
                if (fastAbs(fptt - fret) > del) {
                    del = fastAbs(fptt - fret);
                    ibig = i;
                }

                if (info != null) {
                    info.receiveInfo("====> Iter " + iter + ", " + i + " : " + fret);
                }
            }

            if (2.0 * fastAbs(fp - fret) <= ftol * (fastAbs(fp) + fastAbs(fret))) {
                System.out.println("Iterations: " + iter);
                return fret;
            }

            if (iter == ITMAX) {
                System.out.println("powell exceeding maximum iterations.");
            }

            for (j = 0; j < n; j++) {
                ptt[j] = 2.0 * params[j] - pt[j];
                xit[j] = params[j] - pt[j];
                pt[j] = params[j];
            }

            fptt = func.getValue(ptt);
            if (fptt < fp) {
                double t = 2.0 * (fp - 2.0 * fret + fptt) *
                           MathUtil.sqr(fp - fret - del) - del * MathUtil.sqr(fp - fptt);
                if (t < 0) {
                    fret = linmin(params, xit, func);
                    if (info != null) {
                        info.receiveInfo("====> Sp Iter " + iter + ", " + i + " : " +
                                         fret);
                    }
                    for (j = 0; j < n; j++) {
                        xi[j][ibig] = xi[j][n - 1];
                        xi[j][n - 1] = xit[j];
                    }
                }
            }
        } // for (iter)
    }


    private static double linmin(
        final double p[], final double xi[], final Function func) throws
        CannotEvaluateException {

        final int n = xi.length;
        FunctionData f = new FunctionData() {
            double[] xt = new double[n];

            public double getValueAt(double x) throws CannotEvaluateException {
                for (int j = 0; j < n; j++) {
                    xt[j] = p[j] + x * xi[j];
                }
                return func.getValue(xt);
            }


            public String getDescription() {
                return "";
            }
        };

        MinimumBracket b = new MinimumBracket();
        // FIXME b.bx was 1 inteasd of 0.1, changed to solve local minima problem
        b.ax = 0;
        b.bx = 0.1;
        mnbrak(b, f);

        double[] min = _brentParabolic(f, b.ax, b.cx, 1e-4);
        for (int j = 0; j < n; j++) {
            xi[j] *= min[0];
            p[j] += xi[j];
        }

        return min[1];
    }

}
