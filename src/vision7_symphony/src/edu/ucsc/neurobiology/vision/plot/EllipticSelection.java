package edu.ucsc.neurobiology.vision.plot;

import java.awt.*;
import javax.swing.*;

import cern.colt.matrix.*;
import cern.colt.matrix.impl.*;
import cern.colt.matrix.linalg.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class EllipticSelection
    extends Selection implements Cloneable {

    private boolean visibleGuides = true;
    private double x0, y0, a, b, alpha; // screen
    private double _x0, _y0, _a, _b, _alpha; // plot


    public EllipticSelection(JComponent parentComponent, Axis hAxis, Axis vAxis) {
        super(parentComponent, hAxis, vAxis);
    }


    public EllipticSelection(JComponent parentComponent, Axis hAxis, Axis vAxis,
                             double _x0, double _y0, double _a, double _b, double _alpha) {
        super(parentComponent, hAxis, vAxis);

        this._x0 = _x0;
        this._y0 = _y0;
        this._a = _a;
        this._b = _b;
        this._alpha = _alpha;

        plot2screen();
        params2points();
    }


    public EllipticSelection(JComponent parentComponent, Axis hAxis, Axis vAxis,
                             double[] p) {
        this(parentComponent, hAxis, vAxis, p[1], p[2], p[3], p[4], p[5]);
    }


    public SelectionRegion getSelection() {
        return new ParametricEllipse(_x0, _y0, _a, _b, _alpha);
    }


    public int getNumberOfControlPoints() {
        return 5;
    }


    public void initializeControlPoints(int x, int y) {
        activeCtrlPt = 1;

        a = 0;
        b = 20;
        x0 = x;
        y0 = y;
        alpha = 0 * Math.PI / 180;

        params2points();
        parentComponent.repaint();
    }


    public void updateActiveControlPoint(int x, int y) {
        switch (activeCtrlPt) {
            case 0:
                x0 = x;
                y0 = y;
                break;

            case 1:
                a = getDistance(x, y, 0);
                alpha = Math.atan2(y - yCtrlPts[0], x - xCtrlPts[0]);
                break;

            case 2:
                a = getDistance(x, y, 0);
                alpha = Math.atan2(y - yCtrlPts[0], x - xCtrlPts[0]) + Math.PI;
                break;

            case 3:
                b = getDistance(x, y, 0);
                alpha = Math.atan2(y - yCtrlPts[0], x - xCtrlPts[0]) + 3 * Math.PI / 2;
                break;

            case 4:
                b = getDistance(x, y, 0);
                alpha = Math.atan2(y - yCtrlPts[0], x - xCtrlPts[0]) + Math.PI / 2;
                break;
        }

        if (isValidSelection()) {
            params2points();
            screen2plot();
            parentComponent.repaint();
        }
    }


    static SparseDoubleMatrix2D m2 = new SparseDoubleMatrix2D(new double[][] { {1}
        , {1}
        , {1}
        , {1}
        , {1}
    });
    static SparseDoubleMatrix2D m1 = new SparseDoubleMatrix2D(5, 5);

    private void screen2plot() {
        ParametricEllipse ei = new ParametricEllipse(x0, y0, a, b, alpha);
        for (int i = 0; i < 5; i++) {
            double theta;
            if (i == 4) {
                theta = Math.PI / 6;
            } else {
                theta = i * Math.PI / 2;
            }
            double[] p = ei.getPointFor(theta);
            double x = hAxis.getPlotCoord(p[0]);
            double y = vAxis.getPlotCoord(p[1]);
            m1.set(i, 0, x * x);
            m1.set(i, 1, x * y);
            m1.set(i, 2, y * y);
            m1.set(i, 3, x);
            m1.set(i, 4, y);
        }

        // solve the system and find the polynomial coefficients
        DoubleMatrix2D s = Algebra.DEFAULT.solve(m1, m2);
        double A = s.get(0, 0);
        double B = s.get(1, 0);
        double C = s.get(2, 0);
        double D = s.get(3, 0);
        double E = s.get(4, 0);

        // find the ellipse parameters
        double t = 1 / (B * B - 4 * A * C);
        double M = (C * D * D - B * D * E + A * E * E) * t;
        double w = Math.sqrt(B * B + (A - C) * (A - C));
        _x0 = (2 * D * C - B * E) * t;
        _y0 = (2 * A * E - B * D) * t;
        _a = 1 / Math.sqrt(0.5 * (A + C + w) / (1 - M));
        _b = 1 / Math.sqrt(0.5 * (A + C - w) / (1 - M));
        _alpha = 0.5 * Math.atan2(B, (A - C));
        while (_alpha < 0) {
            _alpha += Math.PI;
        }
        if (!Matrix.correct(ei, new ParametricEllipse(_x0, _y0, _a, _b, _alpha))) {
            for (int i = 0; i < 4; i++) {
                double da = i * Math.PI / 2;
                if (Matrix.correct(ei,
                                   new ParametricEllipse(_x0, _y0, _a, _b,
                    _alpha + da))) {
                    alpha += da;
                    break;
                }
                if (Matrix.correct(ei,
                                   new ParametricEllipse(_x0, _y0, _b, _a,
                    _alpha + da))) {
                    alpha += da;
                    double temp = a;
                    a = b;
                    b = temp;
                    break;
                }
            }
        } while (_alpha > 2 * Math.PI) {
            _alpha -= 2 * Math.PI;
        }

//        System.out.println("--plot---------------------");
//        System.out.println("x0    = " + VisionUtilities.format(_x0, 3));
//        System.out.println("y0    = " + VisionUtilities.format(_y0, 3));
//        System.out.println("a     = " + VisionUtilities.format(_a, 3));
//        System.out.println("b     = " + VisionUtilities.format(_b, 3));
//        System.out.println("alpha = " + VisionUtilities.format(_alpha * 180 / Math.PI, 3));
    }


    private void plot2screen() {
        ParametricEllipse ei = new ParametricEllipse(_x0, _y0, _a, _b, _alpha);
        for (int i = 0; i < 5; i++) {
            double theta;
            if (i == 4) {
                theta = Math.PI / 6;
            } else {
                theta = i * Math.PI / 2;
            }
            double[] p = ei.getPointFor(theta);
            double x = hAxis.getScreenCoord(p[0]);
            double y = vAxis.getScreenCoord(p[1]);
            m1.set(i, 0, x * x);
            m1.set(i, 1, x * y);
            m1.set(i, 2, y * y);
            m1.set(i, 3, x);
            m1.set(i, 4, y);
        }

        // solve the system and find the polynomial coefficients
        DoubleMatrix2D s = Algebra.DEFAULT.solve(m1, m2);
        double A = s.get(0, 0);
        double B = s.get(1, 0);
        double C = s.get(2, 0);
        double D = s.get(3, 0);
        double E = s.get(4, 0);

        // find the ellipse parameters
        double t = 1 / (B * B - 4 * A * C);
        double M = (C * D * D - B * D * E + A * E * E) * t;
        double w = Math.sqrt(B * B + (A - C) * (A - C));
        x0 = (2 * D * C - B * E) * t;
        y0 = (2 * A * E - B * D) * t;
        a = 1 / Math.sqrt(0.5 * (A + C + w) / (1 - M));
        b = 1 / Math.sqrt(0.5 * (A + C - w) / (1 - M));
        alpha = 0.5 * Math.atan2(B, (A - C));
        while (alpha < 0) {
            alpha += Math.PI;
        }
        if (!Matrix.correct(ei, new ParametricEllipse(x0, y0, a, b, alpha))) {
            for (int i = 0; i < 4; i++) {
                double da = i * Math.PI / 2;
                if (Matrix.correct(ei, new ParametricEllipse(x0, y0, a, b, alpha + da))) {
                    alpha += da;
                    break;
                }
                if (Matrix.correct(ei, new ParametricEllipse(x0, y0, b, a, alpha + da))) {
                    alpha += da;
                    double temp = a;
                    a = b;
                    b = temp;
                    break;
                }
            }
        } while (alpha > 2 * Math.PI) {
            alpha -= 2 * Math.PI;
        }

//        System.out.println("--screen---------------------");
//        System.out.println("x0    = " + VisionUtilities.format(x0, 3));
//        System.out.println("y0    = " + VisionUtilities.format(y0, 3));
//        System.out.println("a     = " + VisionUtilities.format(a, 3));
//        System.out.println("b     = " + VisionUtilities.format(b, 3));
//        System.out.println("alpha = " + VisionUtilities.format(alpha * 180 / Math.PI, 3));
    }


    private void params2points() {
        ParametricEllipse e = new ParametricEllipse(x0, y0, a, b, alpha);
        double[] p;

        xCtrlPts[0] = (int) Math.round(x0);
        yCtrlPts[0] = (int) Math.round(y0);

        p = e.getPointFor(0 * Math.PI / 2);
        xCtrlPts[1] = (int) Math.round(p[0]);
        yCtrlPts[1] = (int) Math.round(p[1]);

        p = e.getPointFor(2 * Math.PI / 2);
        xCtrlPts[2] = (int) Math.round(p[0]);
        yCtrlPts[2] = (int) Math.round(p[1]);

        p = e.getPointFor(1 * Math.PI / 2);
        xCtrlPts[3] = (int) Math.round(p[0]);
        yCtrlPts[3] = (int) Math.round(p[1]);

        p = e.getPointFor(3 * Math.PI / 2);
        xCtrlPts[4] = (int) Math.round(p[0]);
        yCtrlPts[4] = (int) Math.round(p[1]);
    }


    public void componentResized() {
//        if (isValidSelection()) {
        plot2screen();
        params2points();
//        }
    }


    private double getDistance(double x, double y, int ctrlPt) {
        double dx = x - xCtrlPts[ctrlPt];
        double dy = y - yCtrlPts[ctrlPt];
        return Math.sqrt(dx * dx + dy * dy);
    }


    public static void drawEllipse(double x1, double y1,
                                   double x2, double y2,
                                   double x3, double y3,
                                   double x4, double y4,
                                   Graphics2D g) {

        double alpha = Math.atan( (y1 - y2) / (x1 - x2));
        double a = Math.sqrt( (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)) / 2;
        double b = Math.sqrt( (x3 - x4) * (x3 - x4) + (y3 - y4) * (y3 - y4)) / 2;
        double x0 = (x1 + x2) / 2;
        double y0 = (y1 + y2) / 2;
        g.rotate(alpha, x0, y0);
        g.drawArc( (int) (x0 - a), (int) (y0 - b), (int) (2 * a), (int) (2 * b), 0, 360);
        g.rotate( -alpha, x0, y0);
    }


    public void paintSelection(Graphics g) {
        try {
            if (visible) {
                // Make a 2D graphics context.
                Graphics2D g2d = (Graphics2D) g;
                g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                                     RenderingHints.VALUE_ANTIALIAS_ON);

                params2points();

                if (visibleGuides) {
                    g2d.setStroke(thickStroke);
                    g.setColor(Color.black);
                    g.drawLine(xCtrlPts[1], yCtrlPts[1], xCtrlPts[0], yCtrlPts[0]);
                    g.drawLine(xCtrlPts[2], yCtrlPts[2], xCtrlPts[0], yCtrlPts[0]);
                    g.drawLine(xCtrlPts[3], yCtrlPts[3], xCtrlPts[0], yCtrlPts[0]);
                    g.drawLine(xCtrlPts[4], yCtrlPts[4], xCtrlPts[0], yCtrlPts[0]);
                    g2d.setStroke(thinStroke);
                    g.setColor(Color.white);
                    g.drawLine(xCtrlPts[1], yCtrlPts[1], xCtrlPts[0], yCtrlPts[0]);
                    g.drawLine(xCtrlPts[2], yCtrlPts[2], xCtrlPts[0], yCtrlPts[0]);
                    g.drawLine(xCtrlPts[3], yCtrlPts[3], xCtrlPts[0], yCtrlPts[0]);
                    g.drawLine(xCtrlPts[4], yCtrlPts[4], xCtrlPts[0], yCtrlPts[0]);
                }

                // draw a oval
                g2d.setStroke(thickStroke);
                g2d.setColor(Color.black);
                drawEllipse(xCtrlPts[1], yCtrlPts[1], xCtrlPts[2], yCtrlPts[2],
                            xCtrlPts[3],
                            yCtrlPts[3], xCtrlPts[4], yCtrlPts[4], g2d);
                g2d.setStroke(thinStroke);
                g.setColor(Color.white);
                drawEllipse(xCtrlPts[1], yCtrlPts[1], xCtrlPts[2], yCtrlPts[2],
                            xCtrlPts[3],
                            yCtrlPts[3], xCtrlPts[4], yCtrlPts[4], g2d);

                drawControlPoints(g2d);

                g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                                     RenderingHints.VALUE_ANTIALIAS_OFF);
            }
        } catch (Exception e) {
            System.out.println("Exception");

//            System.out.println("--plot---------------------");
//            System.out.println("x0    = " + VisionUtilities.format(_x0, 3));
//            System.out.println("y0    = " + VisionUtilities.format(_y0, 3));
//            System.out.println("a     = " + VisionUtilities.format(_a, 3));
//            System.out.println("b     = " + VisionUtilities.format(_b, 3));
//            System.out.println("alpha = " +
//                               VisionUtilities.format(_alpha * 180 / Math.PI, 3));

            System.out.println("--screen---------------------");
            System.out.println("x0    = " + StringUtil.format(x0, 3));
            System.out.println("y0    = " + StringUtil.format(y0, 3));
            System.out.println("a     = " + StringUtil.format(a, 3));
            System.out.println("b     = " + StringUtil.format(b, 3));
            System.out.println("alpha = " + StringUtil.format(alpha * 180 / Math.PI, 3));
        }
    }


    /**
     * Check that the area of the selection is non-zero.
     *
     * @return flag indicating whether the selection is valid */
    public boolean isValidSelection() {
        return visible && a != 0 && b != 0;
    }


    public double[] toArray() {
        return new double[] {SelectionType.ELLIPSE.ordinal(), _x0, _y0, _a, _b, _alpha};
    }
}
