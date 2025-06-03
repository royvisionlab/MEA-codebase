package edu.ucsc.neurobiology.vision.plot;

import java.io.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;


/**
 * A abstract implementation of a graphical selection. Used in PlotPanel.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
abstract public class Selection
    extends MouseAdapter implements MouseMotionListener, Serializable {

    public static enum SelectionType {
        RECTANGLE, ROTATED_RECTANGLE, ELLIPSE
    }


    final protected static Stroke thinStroke =
        new BasicStroke(1.f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 3.f);
    final protected static Stroke thickStroke =
        new BasicStroke(3.f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 5.f);

    // A constant which flags that no control point was near the mouse-pressed event
    final public static int NO_CONTROL_POINT = -1;

    // Flag indicating whether or not the selection box is visible
    protected boolean visible;

    // The maximum distance from a control point the cursor can be and still be selected
    protected final static int hitThreshold = 10;

    // The size of the control point boxes
    protected final static int ctrlPtSize = 1;

    // The number of control points for this component
    protected int nCtrlPts;

    /**
     * Which control point is the active one, or which one can be
     * controlled from the arrow keys on the keyboard? */
    protected int activeCtrlPt;

    /**
     * The x-coordinates of the control points.  The first four of
     * these control points MUST define the outer boundries of the
     * selected region. */
    protected int[] xCtrlPts, yCtrlPts;
    double[] xPlot, yPlot;

    protected JComponent parentComponent;
    protected JPopupMenu menu = new JPopupMenu();
    public boolean active = true;
    Axis hAxis, vAxis;


    public Selection(JComponent parentComponent, Axis hAxis, Axis vAxis) {
        this.parentComponent = parentComponent;
        this.hAxis = hAxis;
        this.vAxis = vAxis;

        // make the selection region invisible
        visible = false;

        // Create the arrays of the x- and y-coordinates.  There must
        // be at least four control points.
        nCtrlPts = Math.max(4, getNumberOfControlPoints());
        xCtrlPts = new int[nCtrlPts];
        yCtrlPts = new int[nCtrlPts];
        xPlot = new double[nCtrlPts];
        yPlot = new double[nCtrlPts];
        activeCtrlPt = NO_CONTROL_POINT;

        // set the default cursor
        setCursor();
    }


    public static Selection makeSelection(
        JComponent parentComponent, Axis hAxis, Axis vAxis, double[] params) {

        switch (SelectionType.values()[ (int) params[0]]) {
            case RECTANGLE:
                return new RectangularSelection(parentComponent, hAxis, vAxis, params);

            case ROTATED_RECTANGLE:
                return new ParallelogramSelection(parentComponent, hAxis, vAxis, params);

            case ELLIPSE:
                return new EllipticSelection(parentComponent, hAxis, vAxis, params);

            default:
                return null;
        }
    }


    public void addPossibleAction(Action action) {
        menu.add(action);
    }


    /**
     * Sets the cursor to whatever the current active control point dictates.
     */
    private void setCursor() {
        /*
                 // active cursor
                 Cursor cursor = getControlPointCursor(activeCtrlPt);
                 // default cursor if no active cursor
                 if (cursor == null) {
            cursor = getControlPointCursor(NO_CONTROL_POINT);
                 }
                 // set only when available
                 if (cursor != null) {
            parentComponent.setCursor(cursor);
                 }*/
    }


//    public void mouseClicked(MouseEvent e) {
//    }


//    public void mouseEntered(MouseEvent e) {
//    }


//    public void mouseExited(MouseEvent e) {
//    }


    /**
     * Handle the mousePressed events.
     */
    public void mousePressed(MouseEvent e) {
        if (!SwingUtilities.isLeftMouseButton(e)) {
            resetSelection();
            return;
        }

        // If the selection box is visible AND the user has
        // clicked near one of the existing control points, make
        // the nearest one the active control point and update the
        // current selection.  Return when finished.
        if (visible) {
            int newCtrlPt = nearWhichControlPoint(e.getX(), e.getY());
            if (newCtrlPt >= 0) {
                activeCtrlPt = newCtrlPt;
                setCursor();
                parentComponent.grabFocus();
                parentComponent.repaint();
                return;
            }
        }

        // User wants to start a new selection.  So first set the
        // flag to make the selection region visible.
        visible = true;

        // The initialize method is responsible for setting all of
        // the control points to reasonable values and for setting
        // which point should be the active one.
        activeCtrlPt = NO_CONTROL_POINT;

        // MOD
//        initializeControlPoints(
//            hAxis.getPlotCoord(forceXCoordinateWithinBounds(e.getX())),
//            vAxis.getPlotCoord(forceYCoordinateWithinBounds(e.getY())));
        initializeControlPoints(e.getX(), e.getY());

        setCursor();

        // Update the display.
        parentComponent.repaint();
    }


    public void mouseReleased(MouseEvent e) {
        if (!SwingUtilities.isLeftMouseButton(e)) {
            return;
        }
        updateActiveControlPoint(e.getX(), e.getY());
        if (isValidSelection()) {
            if (active) {
                menu.show(parentComponent, e.getX(), e.getY());
            }
        } else {
            resetSelection();
            activeCtrlPt = NO_CONTROL_POINT;
        }
        setCursor();
    }


    public void mouseDragged(MouseEvent e) {
        if (!SwingUtilities.isLeftMouseButton(e)) {
            return;
        }
        updateActiveControlPoint(e.getX(), e.getY());
        setCursor();
    }


    /**
     * Changes the active control point according to mouse movements
     */
    public void mouseMoved(MouseEvent e) {
        if (!SwingUtilities.isLeftMouseButton(e)) {
            return;
        }

        if (visible) {
            int newCtrlPt = nearWhichControlPoint(e.getX(), e.getY());
            if (newCtrlPt != activeCtrlPt) {
                activeCtrlPt = newCtrlPt;
                setCursor();
                parentComponent.grabFocus();
                parentComponent.repaint();
                return;
            }
        }
    }


    /**
     * A utility method which forces the x-coordinate to be within the
     * component boundries.
     *
     * @param x x-coordinate to force within boundries
     * @return modified x-value
     */
    private double forceXCoordinateWithinBounds(double x) {
        Insets i = parentComponent.getInsets();
        int xmin = i.left;
        int xmax = parentComponent.getWidth() - i.right - 1;
        return Math.max(Math.min(x, xmax), xmin);
    }


    /**
     * A utility method which forces the y-coordinate to be within the
     * component boundries.
     *
     * @param y y-coordinate to force within boundries
     * @return modified y-value */
    private double forceYCoordinateWithinBounds(double y) {
        Insets i = parentComponent.getInsets();
        int ymin = i.top;
        int ymax = parentComponent.getHeight() - i.bottom - 1;
        return Math.max(Math.min(y, ymax), ymin);
    }


    public int nearWhichControlPoint(MouseEvent e) {
        return nearWhichControlPoint(e.getX(), e.getY());
    }


    /**
     * Check to see if the point (x,y) is near one of the control
     * points.  If it is, return the index of the nearest one,
     * otherwise return NO_CONTROL_POINT.
     *
     * @param x x-coordinate to compare to control points
     * @param y y-coordinate to compare to control points
     * still selects it
     *
     * @return the index of the nearest control point */
    public int nearWhichControlPoint(int x, int y) {
        // Initialize to no control point selected.
        int nearestCtrlPt = NO_CONTROL_POINT;
        double minDist2 = -1;

        // Loop over all control points and get the closest one.
        // (Actually calculate distance-squared here.)
        for (int i = 0; i < nCtrlPts; i++) {
            double dx = x - xCtrlPts[i];
            double dy = y - yCtrlPts[i];
            double dist = dx * dx + dy * dy;
            if (dist < minDist2 || i == 0) {
                minDist2 = dist;
                nearestCtrlPt = i;
            }
        }

        // If the closest one isn't close enough, delete the index.
        if (minDist2 > hitThreshold * hitThreshold) {
            nearestCtrlPt = NO_CONTROL_POINT;
        }
        return nearestCtrlPt;
    }


    /**
     * Make the selection box invisible. */
    public void resetSelection() {
        visible = false;
        activeCtrlPt = NO_CONTROL_POINT;
        setCursor();
        parentComponent.repaint();
    }


    /**
     * This returns whether the current selected region is valid.
     * Generally if the area has zero volume, then this method should
     * return false.
     */
    abstract public boolean isValidSelection();


    /**
     * Initialize the control points.  Subclasses must provide an
     * implementation of this method which initializes the control
     * points to reasonable values given the first mouse-pressed
     * coordinates, and must also set the activeCtrlPt to the index of
     * the control point which should be active.
     *
     * @param x x-coordinate of initial mouse-pressed event
     * @param y y-coordinate of initial mouse-pressed event */
    abstract public void initializeControlPoints(int x, int y);


    /**
     * Change the active control point to the point (x,y).  Subclasses
     * should implement this routine to get the behaviour which is
     * desired.  This is the place to impose constraints on how the
     * control points can move.  NOTE: repaintPanel() should be called
     * at the end of this method to update the display.
     *
     * @param x x-coordinate of the new point
     * @param y y-coordinate of the new point */
    abstract public void updateActiveControlPoint(int x, int y);


    /**
     * Must paint the selection.
     *
     * @param g Graphics
     */
    abstract public void paintSelection(Graphics g);


    /**
     * Useful subclasses must define the number of control points on
     * the selected region.  The first four control points define the
     * outer extent of the selected region.  This method MUST NOT
     * return a number less than four. */
    abstract public int getNumberOfControlPoints();


    /**
     * Returns the Cursor to be displayed for a certain control point
     * and the default cursor for this Selection for an index of
     * NO_CONTROL_POINT. Return of null will not change the cursor.
     * Subclasses should override this method to provide a default
     * cursor and/or to provide cursors for the different control points.
     */
    public Cursor getControlPointCursor(int index) {
        return null;
    }


    /**
     * Returns the region selected.
     *
     * @return SelectionRegion
     */
    public abstract SelectionRegion getSelection();


    public Selection getCopy() {
        Selection s = null;
        try {
            s = (Selection)this.clone();
            s.xCtrlPts = (int[])this.xCtrlPts.clone();
            s.yCtrlPts = (int[])this.yCtrlPts.clone();
            s.xPlot = (double[])this.xPlot.clone();
            s.yPlot = (double[])this.yPlot.clone();
            s.visible = true;
            s.active = false;
        } catch (CloneNotSupportedException e) {
            e.printStackTrace();
        }
        return s;
    }


    public void printControlPoints() {
        System.out.println("-------------");
        for (int i = 0; i < nCtrlPts; i++) {
            System.out.println(xCtrlPts[i] + ", " + yCtrlPts[i]);
        }
        for (int i = 0; i < nCtrlPts; i++) {
            System.out.println(xPlot[i] + ", " + yPlot[i]);
        }
    }


    private final static int d = 2;

    protected void drawControlPoints(Graphics2D g2d) {
        for (int i = 0; i < nCtrlPts; i++) {
            g2d.setStroke(thickStroke);
            g2d.setColor(Color.black);
            g2d.drawOval( (int) xCtrlPts[i] - d, (int) yCtrlPts[i] - d, 2 * d, 2 * d);
            g2d.setStroke(thinStroke);
            g2d.setColor(Color.white);
            g2d.drawOval( (int) xCtrlPts[i] - d, (int) yCtrlPts[i] - d, 2 * d, 2 * d);
        }

//        for (int i = 0; i < nCtrlPts; i++) {
//            g2d.setStroke(thickStroke);
//            g2d.setColor(PlotUtilities.getColor(i));
//            g2d.drawArc( (int) xCtrlPts[i] - d, (int) yCtrlPts[i] - d, 2 * d, 2 * d, 0,
//                        360);
//        }
    }


    public void convertToPlot() {
        for (int i = 0; i < nCtrlPts; i++) {
            xPlot[i] = hAxis.getPlotCoord(xCtrlPts[i]);
            yPlot[i] = vAxis.getPlotCoord(yCtrlPts[i]);
        }
    }


    public void componentResized() {
        for (int i = 0; i < nCtrlPts; i++) {
            xCtrlPts[i] = (int) hAxis.getScreenCoord(xPlot[i]);
            yCtrlPts[i] = (int) vAxis.getScreenCoord(yPlot[i]);
        }
    }


    /**
     * Must returns the selection as a simple list of points. The first element in the
     * array must be a unique integer serving as a selection ID.
     *
     * @return double[]
     */
    public abstract double[] toArray();
}
