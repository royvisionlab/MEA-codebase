package edu.ucsc.neurobiology.vision.plot;

import java.util.*;

import java.awt.*;
import javax.swing.*;

import edu.ucsc.neurobiology.vision.math.*;


/**
 * Selects a rectangular screen region.  The sides of the rectangle
 * are parallel to the x and y-axes.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class RectangularSelection
    extends Selection implements Cloneable {


    public RectangularSelection(JComponent parentComponent, Axis hAxis, Axis vAxis) {
        super(parentComponent, hAxis, vAxis);
    }


    public RectangularSelection(JComponent parentComponent, Axis hAxis, Axis vAxis,
                                double x1, double y1, double x2, double y2) {
        super(parentComponent, hAxis, vAxis);

        xPlot[0] = x1;
        yPlot[0] = y1;

        xPlot[1] = x1;
        yPlot[1] = y2;

        xPlot[2] = x2;
        yPlot[2] = y2;

        xPlot[3] = x2;
        yPlot[3] = y1;

        for (int i = 0; i < nCtrlPts; i++) {
            xCtrlPts[i] = (int) Math.round(hAxis.getScreenCoord(xPlot[i]));
            yCtrlPts[i] = (int) Math.round(vAxis.getScreenCoord(yPlot[i]));
        }

        parentComponent.repaint();
    }


    public RectangularSelection(JComponent parentComponent, Axis hAxis, Axis vAxis,
                                double[] p) {
        this(parentComponent, hAxis, vAxis, p[1], p[2], p[3], p[4]);
    }


    public void initializeControlPoints(int x, int y) {
        activeCtrlPt = 2;
        Arrays.fill(xCtrlPts, x);
        Arrays.fill(yCtrlPts, y);
    }


    /**
     * Returns the number of control points for the rectangle (4).
     *
     * @return 4, the control points at the corner of the rectangle
     */
    public int getNumberOfControlPoints() {
        return 4;
    }


    public SelectionRegion getSelection() {
        return new Polygon2D(xPlot, yPlot, xPlot.length);
    }


    /*
        public Cursor getControlPointCursor(int index) {
            if (index >= 0) {
                int k = (index + 2) % 4;
                return compassCursor(
                    "Resize", xCtrlPts[index] - xCtrlPts[k],
                    yCtrlPts[index] - yCtrlPts[k], 4, true);
            }
            return FreeHepImage.getCursor("RectangularCursor");
        }
     */


    public void updateActiveControlPoint(int x, int y) {
        if (activeCtrlPt == -1) {
            return;
        }

//        System.out.println(activeCtrlPt);

        xCtrlPts[activeCtrlPt] = x;
        yCtrlPts[activeCtrlPt] = y;

        switch (activeCtrlPt) {
            case 0:
                xCtrlPts[1] = x;
                yCtrlPts[3] = y;
                break;

            case 1:
                xCtrlPts[0] = x;
                yCtrlPts[2] = y;
                break;

            case 2:
                xCtrlPts[3] = x;
                yCtrlPts[1] = y;
                break;

            case 3:
                xCtrlPts[2] = x;
                yCtrlPts[0] = y;
                break;

            default:
                break;
        }

        convertToPlot();
        parentComponent.repaint();
    }


    public void paintSelection(Graphics _g) {
//        printControlPoints();

        if (visible) {
            Graphics2D g2d = (Graphics2D) _g;
            g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                                 RenderingHints.VALUE_ANTIALIAS_ON);

            for (int i = 0; i < nCtrlPts; i++) {
                xCtrlPts[i] = (int) Math.round(hAxis.getScreenCoord(xPlot[i]));
                yCtrlPts[i] = (int) Math.round(vAxis.getScreenCoord(yPlot[i]));
            }

            // Draw a rectangle on top the the image.
            g2d.setStroke(thickStroke);
            g2d.setColor(Color.black);
            g2d.drawLine(xCtrlPts[0], yCtrlPts[0], xCtrlPts[1], yCtrlPts[1]);
            g2d.drawLine(xCtrlPts[1], yCtrlPts[1], xCtrlPts[2], yCtrlPts[2]);
            g2d.drawLine(xCtrlPts[2], yCtrlPts[2], xCtrlPts[3], yCtrlPts[3]);
            g2d.drawLine(xCtrlPts[3], yCtrlPts[3], xCtrlPts[0], yCtrlPts[0]);

            g2d.setStroke(thinStroke);
            g2d.setColor(Color.white);
            g2d.drawLine(xCtrlPts[0], yCtrlPts[0], xCtrlPts[1], yCtrlPts[1]);
            g2d.drawLine(xCtrlPts[1], yCtrlPts[1], xCtrlPts[2], yCtrlPts[2]);
            g2d.drawLine(xCtrlPts[2], yCtrlPts[2], xCtrlPts[3], yCtrlPts[3]);
            g2d.drawLine(xCtrlPts[3], yCtrlPts[3], xCtrlPts[0], yCtrlPts[0]);

            drawControlPoints(g2d);

            g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                                 RenderingHints.VALUE_ANTIALIAS_OFF);
        }
    }


    /**
     * Determine if the selection is valid.  Return false if the area
     * is zero.
     */
    public boolean isValidSelection() {
        return (visible) && (xCtrlPts[0] != xCtrlPts[2] || yCtrlPts[0] != yCtrlPts[2]);
    }


    public double[] toArray() {
        double xMin = MathUtil.min(xPlot);
        double yMin = MathUtil.min(yPlot);
        double xMax = MathUtil.max(xPlot);
        double yMax = MathUtil.max(yPlot);
        return new double[] {SelectionType.RECTANGLE.ordinal(), xMin, yMin, xMax, yMax};
    }

}
