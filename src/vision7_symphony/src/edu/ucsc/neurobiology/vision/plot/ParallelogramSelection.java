package edu.ucsc.neurobiology.vision.plot;

import java.util.*;

import java.awt.*;
import javax.swing.*;


/**
 * Selects a screen region shaped like a parallelogram.
 *
 * @author Charles Loomis
 * @author Mark Donszelmann
 */
public class ParallelogramSelection
    extends Selection implements Cloneable {

    boolean visibleGuides = false;

    /**
     * The initial starting width (in pixels) of the first and last
     * sides. */
    final private static int STARTING_WIDTH = 25;


    public ParallelogramSelection(JComponent parentComponent, Axis hAxis, Axis vAxis) {
        super(parentComponent, hAxis, vAxis);
    }


    public ParallelogramSelection(JComponent parentComponent, Axis hAxis, Axis vAxis,
                                  double[] p) {
        super(parentComponent, hAxis, vAxis);

        for (int i = 0, j = 1; i < 6; i++) {
            xPlot[i] = p[j++];
            yPlot[i] = p[j++];
        }

        for (int i = 0; i < nCtrlPts; i++) {
            xCtrlPts[i] = (int) Math.round(hAxis.getScreenCoord(xPlot[i]));
            yCtrlPts[i] = (int) Math.round(vAxis.getScreenCoord(yPlot[i]));
        }
        parentComponent.repaint();
    }


    public SelectionRegion getSelection() {
        return new Polygon2D(xPlot, yPlot, 4);
    }


    /**
     * The number of control points is 6 for the parallelogram.  The
     * four corners and two at the centerpoints of the first and last
     * sides.
     *
     * @return 6 the number of control points */
    public int getNumberOfControlPoints() {
        return 6;
    }


    /*
        public Cursor getControlPointCursor(int index) {
            int k;
            switch (index) {
                case 0:
                case 3:
                case 5:
                    k = 4;
                    break;
                case 1:
                case 2:
                case 4:
                    k = 5;
                    break;
                default:
                    return FreeHepImage.getCursor("ParallelogramCursor");
            }
            return compassCursor("Rotation", xCtrlPts[index] - xCtrlPts[k],
                                 yCtrlPts[index] - yCtrlPts[k], 8, true);
        }
     */

    /**
     * Initialize the control points for this selection given the
     * initial starting point (x,y).
     *
     * @param x the initial x-coordinate
     * @param y the initial y-coordinate */
    public void initializeControlPoints(int x, int y) {
        // Set the fifth control point to be the active one and
        // initialize all of the coordinates.
        activeCtrlPt = 4;
        Arrays.fill(yCtrlPts, y);
        xCtrlPts[0] = x - STARTING_WIDTH;
        xCtrlPts[1] = x - STARTING_WIDTH;
        xCtrlPts[2] = x + STARTING_WIDTH;
        xCtrlPts[3] = x + STARTING_WIDTH;
        xCtrlPts[4] = x;
        xCtrlPts[5] = x;
    }


    /**
     * A utility routine to get the radius from one of the control points.
     */
    private double getRadius(double x, double y, int ctrlPt) {
        double dx = x - xCtrlPts[ctrlPt];
        double dy = y - yCtrlPts[ctrlPt];
        return Math.sqrt(dx * dx + dy * dy);
    }


    /**
     * A utility routine to get the angle of the rotated rectangle.
     */
    private double getAngle() {
        // Get the angle of the rotated rectangle.
        double deltax = xCtrlPts[5] - xCtrlPts[4];
        double deltay = yCtrlPts[5] - yCtrlPts[4];
        if (deltax != 0 || deltay != 0) {
            return Math.atan2(deltay, deltax);
        } else {
            return 0.;
        }
    }


    /**
     * Move the active control point to the point (x,y).
     *
     * @param x x-coordinate of the new point
     * @param y y-coordinate of the new point */
    public void updateActiveControlPoint(int x, int y) {
        // Bring the location within bounds.
//        x = forceXCoordinateWithinBounds(x);
//        y = forceYCoordinateWithinBounds(y);

        // Change what is done depending on which control point is active.
        int width;
        double radius, angle;
        int dx, dy, dx1, dy1;

        switch (activeCtrlPt) {
            case 4:

                /*
                                 dx = x - xCtrlPts[5];
                                 dy = y - yCtrlPts[5];
                                 dx1 = xCtrlPts[4] - xCtrlPts[5];
                                 dy1 = yCtrlPts[4] - yCtrlPts[5];
                                 double cos = (dx * dx1 + dy * dy1) /
                             Math.sqrt(dx * dx + dy * dy) /
                             Math.sqrt(dx1 * dx1 + dy1 * dy1);
                                 double sin = Math.sqrt(1 - cos * cos);
                                 double x0 = xCtrlPts[5];
                                 double y0 = xCtrlPts[5];
                                 for (int i = 0; i < nCtrlPts; i++) {
                    xCtrlPts[i] += (xCtrlPts[i] - x0) * cos + (yCtrlPts[i] - y0) * sin;
                    yCtrlPts[i] += (yCtrlPts[i] - y0) * cos - (xCtrlPts[i] - x0) * sin;
                                 }
                 */

                dx = xCtrlPts[0] - xCtrlPts[4];
                dy = yCtrlPts[0] - yCtrlPts[4];

                // Update the active control point.
                xCtrlPts[activeCtrlPt] = x;
                yCtrlPts[activeCtrlPt] = y;

                // Update the control points on either side.
                xCtrlPts[0] = x + dx;
                yCtrlPts[0] = y + dy;
                xCtrlPts[3] = x - dx;
                yCtrlPts[3] = y - dy;
                break;

            case 5:
                dx = xCtrlPts[1] - xCtrlPts[5];
                dy = yCtrlPts[1] - yCtrlPts[5];

                // Update the active control point.
                xCtrlPts[activeCtrlPt] = x;
                yCtrlPts[activeCtrlPt] = y;

                // Update the control points on either side.
                xCtrlPts[1] = x + dx;
                yCtrlPts[1] = y + dy;
                xCtrlPts[2] = x - dx;
                yCtrlPts[2] = y - dy;
                break;

            case 0:
            case 3:
                xCtrlPts[activeCtrlPt] = x;
                yCtrlPts[activeCtrlPt] = y;

                // Determine the delta-x and delta-y.
                dx = xCtrlPts[activeCtrlPt] - xCtrlPts[4];
                dy = yCtrlPts[activeCtrlPt] - yCtrlPts[4];

                if (activeCtrlPt == 3) {
                    dx = -dx;
                    dy = -dy;
                }

                xCtrlPts[1] = xCtrlPts[5] + dx;
                yCtrlPts[1] = yCtrlPts[5] + dy;
                xCtrlPts[2] = xCtrlPts[5] - dx;
                yCtrlPts[2] = yCtrlPts[5] - dy;

                xCtrlPts[0] = xCtrlPts[4] + dx;
                yCtrlPts[0] = yCtrlPts[4] + dy;
                xCtrlPts[3] = xCtrlPts[4] - dx;
                yCtrlPts[3] = yCtrlPts[4] - dy;
                break;

            case 1:
            case 2:
                xCtrlPts[activeCtrlPt] = x;
                yCtrlPts[activeCtrlPt] = y;

                // Determine the delta-x and delta-y.
                dx = xCtrlPts[activeCtrlPt] - xCtrlPts[5];
                dy = yCtrlPts[activeCtrlPt] - yCtrlPts[5];

                if (activeCtrlPt == 2) {
                    dx = -dx;
                    dy = -dy;
                }

                xCtrlPts[1] = xCtrlPts[5] + dx;
                yCtrlPts[1] = yCtrlPts[5] + dy;
                xCtrlPts[2] = xCtrlPts[5] - dx;
                yCtrlPts[2] = yCtrlPts[5] - dy;

                xCtrlPts[0] = xCtrlPts[4] + dx;
                yCtrlPts[0] = yCtrlPts[4] + dy;
                xCtrlPts[3] = xCtrlPts[4] - dx;
                yCtrlPts[3] = yCtrlPts[4] - dy;

                break;

            default:
                break;
        }

        convertToPlot();
        parentComponent.repaint();
    }


    /**
     * Repaint this component.
     *
     * @param g Graphics context in which to draw */
    public void paintSelection(Graphics g) {

        // If the selection region is visible, paint it.
        if (visible) {
            // Make a 2D graphics context.
            Graphics2D g2d = (Graphics2D) g;
            g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                                 RenderingHints.VALUE_ANTIALIAS_ON);

            for (int i = 0; i < xCtrlPts.length; i++) {
                xCtrlPts[i] = (int) Math.round(hAxis.getScreenCoord(xPlot[i]));
                yCtrlPts[i] = (int) Math.round(vAxis.getScreenCoord(yPlot[i]));
            }

            // Draw a rectangle on top the the image.
            g2d.setStroke(thickStroke);
            g.setColor(Color.black);
            g.drawPolygon(xCtrlPts, yCtrlPts, 4);

            if (visibleGuides) {
                g.drawLine(xCtrlPts[4], yCtrlPts[4], xCtrlPts[5], yCtrlPts[5]);
            }

            g2d.setStroke(thinStroke);
            g.setColor(Color.white);
            g.drawPolygon(xCtrlPts, yCtrlPts, 4);

            if (visibleGuides) {
                g.drawLine(xCtrlPts[4], yCtrlPts[4], xCtrlPts[5], yCtrlPts[5]);
            }

            drawControlPoints(g2d);

            g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                                 RenderingHints.VALUE_ANTIALIAS_OFF);
        }

    }


    /**
     * Returns a boolean indicating whether or not the selected region
     * is valid.  It is valid only if the region has a non-zero area.
     *
     * @return flag indicating whether the region is valid
     */
    public boolean isValidSelection() {
        return (xCtrlPts[4] != xCtrlPts[5] || yCtrlPts[4] != yCtrlPts[5]);
    }


    public double[] toArray() {
        return new double[] {
            SelectionType.ROTATED_RECTANGLE.ordinal(),
            xPlot[0], yPlot[0],
            xPlot[1], yPlot[1],
            xPlot[2], yPlot[2],
            xPlot[3], yPlot[3],
            xPlot[4], yPlot[4],
            xPlot[5], yPlot[5]};
    }
}
