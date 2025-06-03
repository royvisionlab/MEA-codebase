package edu.ucsc.neurobiology.vision.plot;

import java.util.*;

import java.awt.*;
import javax.swing.*;


/**
 * A panel which selects a rectangular region on the screen which can
 * be arbitrarily rotated.
 *
 * @author Charles Loomis
 * @author Mark Donszelmann
 */
public class RotatedRectangularSelection
    extends Selection implements Cloneable {

    boolean visibleGuides = false;

    /**
     * The initial starting width of the first and last sides. */
    final private static int STARTING_WIDTH = 25;

    /**
     * Creates a RotatedRectangleSelectionPanel. */
    public RotatedRectangularSelection(JComponent parentComponent, Axis hAxis, Axis vAxis) {
        super(parentComponent, hAxis, vAxis);
    }


    public SelectionRegion getSelection() {
        return new Polygon2D(xPlot, yPlot, 4);
    }


    /**
     * The number of control points is 6---the four corners and the
     * centerpoints of the first and last sides.
     *
     * @return 6 the number of control points */
    public int getNumberOfControlPoints() {
        return 6;
    }


    /*
        public Cursor getControlPointCursor(int index) {
            int k;
            String type;
            switch (index) {
                case 0:
                case 3:
                    k = 4;
                    type = "Resize";
                    break;
                case 1:
                case 2:
                    k = 5;
                    type = "Resize";
                    break;
                case 4:
                    k = 5;
                    type = "Rotation";
                    break;
                case 5:
                    k = 4;
                    type = "Rotation";
                    break;
                default:
                    return FreeHepImage.getCursor("RotatedRectangleCursor");
            }
            return compassCursor(type, xCtrlPts[index] - xCtrlPts[k],
                                 yCtrlPts[index] - yCtrlPts[k], 8, true);
        }
     */

    /**
     * Initialize the control points given the starting point (x,y).
     *
     * @param x x-coordinate of the starting point
     * @param y y-coordinate of the starting point */
    public void initializeControlPoints(int x, int y) {
        // Set the fifth control point to be the active one and
        // initialize all of the coordinates.
        activeCtrlPt = 5;
        Arrays.fill(yCtrlPts, y);
        xCtrlPts[0] = x - STARTING_WIDTH;
        xCtrlPts[1] = x - STARTING_WIDTH;
        xCtrlPts[2] = x + STARTING_WIDTH;
        xCtrlPts[3] = x + STARTING_WIDTH;
        xCtrlPts[4] = x;
        xCtrlPts[5] = x;

        convertToPlot();
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

        // Change what is done depending on which control point is
        // active.
        int width;
        int dx;
        int dy;
        int deltax;
        int deltay;
        double angle;
        double radius;

        switch (activeCtrlPt) {
            case 0:

                // Determine the radius of the corners from the middle of
                // the control sides.
                radius = getRadius(x, y, 4);

                // Get the angle of the rotated rectangle.
                angle = getAngle();

                deltax = x - xCtrlPts[4];
                deltay = y - yCtrlPts[4];

                dx = (int) ( -radius * Math.sin(angle));
                dy = (int) (radius * Math.cos(angle));

                if (deltax * dx + deltay * dy < 0) {
                    dx = -dx;
                    dy = -dy;
                }

                // Update the other control points.
                xCtrlPts[0] = xCtrlPts[4] + dx;
                yCtrlPts[0] = yCtrlPts[4] + dy;
                xCtrlPts[3] = xCtrlPts[4] - dx;
                yCtrlPts[3] = yCtrlPts[4] - dy;

                xCtrlPts[1] = xCtrlPts[5] + dx;
                yCtrlPts[1] = yCtrlPts[5] + dy;
                xCtrlPts[2] = xCtrlPts[5] - dx;
                yCtrlPts[2] = yCtrlPts[5] - dy;
                break;

            case 1:

                // Determine the radius of the corners from the middle of
                // the control sides.
                radius = getRadius(x, y, 5);

                // Get the angle of the rotated rectangle.
                angle = getAngle();

                deltax = x - xCtrlPts[5];
                deltay = y - yCtrlPts[5];

                dx = (int) ( -radius * Math.sin(angle));
                dy = (int) (radius * Math.cos(angle));

                if (deltax * dx + deltay * dy < 0) {
                    dx = -dx;
                    dy = -dy;
                }

                // Update the other control points.
                xCtrlPts[0] = xCtrlPts[4] + dx;
                yCtrlPts[0] = yCtrlPts[4] + dy;
                xCtrlPts[3] = xCtrlPts[4] - dx;
                yCtrlPts[3] = yCtrlPts[4] - dy;

                xCtrlPts[1] = xCtrlPts[5] + dx;
                yCtrlPts[1] = yCtrlPts[5] + dy;
                xCtrlPts[2] = xCtrlPts[5] - dx;
                yCtrlPts[2] = yCtrlPts[5] - dy;
                break;

            case 2:

                // Determine the radius of the corners from the middle of
                // the control sides.
                radius = getRadius(x, y, 5);

                // Get the angle of the rotated rectangle.
                angle = getAngle();

                deltax = x - xCtrlPts[5];
                deltay = y - yCtrlPts[5];

                dx = (int) ( -radius * Math.sin(angle));
                dy = (int) (radius * Math.cos(angle));

                if (deltax * dx + deltay * dy < 0) {
                    dx = -dx;
                    dy = -dy;
                }

                // Update the other control points.
                xCtrlPts[0] = xCtrlPts[4] - dx;
                yCtrlPts[0] = yCtrlPts[4] - dy;
                xCtrlPts[3] = xCtrlPts[4] + dx;
                yCtrlPts[3] = yCtrlPts[4] + dy;

                xCtrlPts[1] = xCtrlPts[5] - dx;
                yCtrlPts[1] = yCtrlPts[5] - dy;
                xCtrlPts[2] = xCtrlPts[5] + dx;
                yCtrlPts[2] = yCtrlPts[5] + dy;
                break;

            case 3:

                // Determine the radius of the corners from the middle of
                // the control sides.
                radius = getRadius(x, y, 4);

                // Get the angle of the rotated rectangle.
                angle = getAngle();

                deltax = x - xCtrlPts[4];
                deltay = y - yCtrlPts[4];

                dx = (int) ( -radius * Math.sin(angle));
                dy = (int) (radius * Math.cos(angle));

                if (deltax * dx + deltay * dy < 0) {
                    dx = -dx;
                    dy = -dy;
                }

                // Update the other control points.
                xCtrlPts[0] = xCtrlPts[4] - dx;
                yCtrlPts[0] = yCtrlPts[4] - dy;
                xCtrlPts[3] = xCtrlPts[4] + dx;
                yCtrlPts[3] = yCtrlPts[4] + dy;

                xCtrlPts[1] = xCtrlPts[5] - dx;
                yCtrlPts[1] = yCtrlPts[5] - dy;
                xCtrlPts[2] = xCtrlPts[5] + dx;
                yCtrlPts[2] = yCtrlPts[5] + dy;
                break;

            case 4:
            case 5:

                // Determine the radius of the corners from the middle of
                // the control sides.
                radius = getRadius(xCtrlPts[0], yCtrlPts[0], 4);

                // Update the active control point.
                xCtrlPts[activeCtrlPt] = x;
                yCtrlPts[activeCtrlPt] = y;

                // Get the angle of the rotated rectangle.
                angle = getAngle();

                dx = (int) ( -radius * Math.sin(angle));
                dy = (int) (radius * Math.cos(angle));

                // Update the other control points.
                xCtrlPts[0] = xCtrlPts[4] + dx;
                yCtrlPts[0] = yCtrlPts[4] + dy;
                xCtrlPts[3] = xCtrlPts[4] - dx;
                yCtrlPts[3] = yCtrlPts[4] - dy;

                xCtrlPts[1] = xCtrlPts[5] + dx;
                yCtrlPts[1] = yCtrlPts[5] + dy;
                xCtrlPts[2] = xCtrlPts[5] - dx;
                yCtrlPts[2] = yCtrlPts[5] - dy;
                break;

            default:
                break;
        }

        convertToPlot();

        parentComponent.repaint();
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


    public void paintSelection(Graphics g) {
        // If the selection region is visible, paint it.
        if (visible) {
            // Make a 2D graphics context.
            Graphics2D g2d = (Graphics2D) g;
            g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                                 RenderingHints.VALUE_ANTIALIAS_ON);

            for (int i = 0; i < xCtrlPts.length; i++) {
                xCtrlPts[i] = (int) hAxis.getScreenCoord(xPlot[i]);
                yCtrlPts[i] = (int) vAxis.getScreenCoord(yPlot[i]);
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
     * Check that the area of the selection is non-zero.
     *
     * @return flag indicating whether the selection is valid */
    public boolean isValidSelection() {
        return (visible) &&
            (xCtrlPts[4] != xCtrlPts[5] || yCtrlPts[4] != yCtrlPts[5]) &&
            (xCtrlPts[0] != xCtrlPts[3] || yCtrlPts[0] != yCtrlPts[3]);
    }


    public double[] toArray() {
        return new double[] {};
    }
}
