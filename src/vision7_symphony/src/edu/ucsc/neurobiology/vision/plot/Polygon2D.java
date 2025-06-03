package edu.ucsc.neurobiology.vision.plot;

import java.util.*;

import java.awt.*;
import java.awt.geom.*;

import edu.ucsc.neurobiology.vision.math.*;
import sun.awt.geom.*;


/**
 * The <code>Polygon2D</code> class encapsulates a description of a
 * closed, two-dimensional region within a coordinate space. This
 * region is bounded by an arbitrary number of line segments, each of
 * which is one side of the polygon. Internally, a polygon
 * comprises of a list of (<i>x</i>,&nbsp;<i>y</i>)
 * coordinate pairs, where each pair defines a <i>vertex</i> of the
 * polygon, and two successive pairs are the endpoints of a
 * line that is a side of the polygon. The first and final
 * pairs of (<i>x</i>,&nbsp;<i>y</i>) points are joined by a line segment
 * that closes the polygon.  This <code>Polygon2D</code> is defined with
 * an even-odd winding rule.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 * code copied from java.awt and modified
 * @author Matthew Grivich, University of California, Santa Cruz
 * @owner Matthew Grivich
 */
public class Polygon2D
    implements Shape, SelectionRegion {

    /**
     * The total number of points.
     * This value can be NULL.
     */
    public int nPoints = 0;

    /**
     * The array of <i>x</i> coordinates.
     */
    public double xPoints[] = new double[4];

    /**
     * The array of <i>y</i> coordinates.
     */
    public double yPoints[] = new double[4];

    /**
     * Bounds of the polygon.
     * This value can be NULL.
     * Please see the javadoc comments getBounds().
     */
    private Rectangle2D.Double bounds = null;

    private boolean close = false;


    /**
     * Creates an empty polygon.
     */
    public Polygon2D(boolean close) {
        this.close = close;
    }


    public Polygon2D() { }


    /**
     * Constructs and initializes a <code>Polygon2D</code> from the specified
     * parameters.
     * @param xPoints an array of <i>x</i> coordinates
     * @param yPoints an array of <i>y</i> coordinates
     * @param nPoints the total number of points
     */
    public Polygon2D(double xPoints[], double yPoints[], int nPoints) {
        this.nPoints = nPoints;
        this.xPoints = new double[nPoints];
        this.yPoints = new double[nPoints];
        System.arraycopy(xPoints, 0, this.xPoints, 0, nPoints);
        System.arraycopy(yPoints, 0, this.yPoints, 0, nPoints);
    }


    /**
     * Translates the vertices of the <code>Polygon2D</code> by
     * <code>deltaX</code> along the x axis and by
     * <code>deltaY</code> along the y axis.
     * @param deltaX the amount to translate along the <i>x</i> axis
     * @param deltaY the amount to translate along the <i>y</i> axis
     */
    public void translate(double deltaX, double deltaY) {
        for (int i = 0; i < nPoints; i++) {
            xPoints[i] += deltaX;
            yPoints[i] += deltaY;
        }
        if (bounds != null) {
            bounds.x += deltaX;
            bounds.y += deltaY;
        }
    }


    /*
     * Calculates the bounding box of the points passed to the constructor.
     * Sets <code>bounds</code> to the result.
     * @param xPoints[] array of <i>x</i> coordinates
     * @param yPoints[] array of <i>y</i> coordinates
     * @param nPoints the total number of points
     */
    private void calculateBounds(double xPoints[], double yPoints[], int nPoints) {
        double boundsMinX = Double.POSITIVE_INFINITY;
        double boundsMinY = Double.POSITIVE_INFINITY;
        double boundsMaxX = Double.NEGATIVE_INFINITY;
        double boundsMaxY = Double.NEGATIVE_INFINITY;

        for (int i = 0; i < nPoints; i++) {
            if (xPoints[i] < boundsMinX) {
                boundsMinX = xPoints[i];
            }
            if (xPoints[i] > boundsMaxX) {
                boundsMaxX = xPoints[i];
            }
            if (yPoints[i] < boundsMinY) {
                boundsMinY = yPoints[i];
            }
            if (yPoints[i] > boundsMaxY) {
                boundsMaxY = yPoints[i];
            }
//            double x = xPoints[i];
//            boundsMinX = Math.min(boundsMinX, x);
//            boundsMaxX = Math.max(boundsMaxX, x);
//            double y = yPoints[i];
//            boundsMinY = Math.min(boundsMinY, y);
//            boundsMaxY = Math.max(boundsMaxY, y);
        }

        bounds = new Rectangle2D.Double(
            boundsMinX, boundsMinY,
            boundsMaxX - boundsMinX,
            boundsMaxY - boundsMinY);
    }


    /*
     * Resizes the bounding box to accomodate the specified coordinates.
     * @param x,&nbsp;y the specified coordinates
     */
    private void updateBounds(double x, double y) {
        if (x < bounds.x) {
            bounds.width = bounds.width + (bounds.x - x);
            bounds.x = x;
        } else {
            bounds.width = Math.max(bounds.width, x - bounds.x);
            // bounds.x = bounds.x;
        }

        if (y < bounds.y) {
            bounds.height = bounds.height + (bounds.y - y);
            bounds.y = y;
        } else {
            bounds.height = Math.max(bounds.height, y - bounds.y);
            // bounds.y = bounds.y;
        }
    }


    /**
     * Appends the specified coordinates to this <code>Polygon2D</code>.
     * If an operation that calculates the bounding box of this
     * <code>Polygon2D</code> has already been performed, such as
     * <code>getBounds</code> or <code>contains</code>, then this
     * method updates the bounding box.
     */
    public void addPoint(double x, double y) {
        if (nPoints == xPoints.length) {
            double tmp[];

            tmp = new double[nPoints * 2];
            System.arraycopy(xPoints, 0, tmp, 0, nPoints);
            xPoints = tmp;

            tmp = new double[nPoints * 2];
            System.arraycopy(yPoints, 0, tmp, 0, nPoints);
            yPoints = tmp;
        }
        xPoints[nPoints] = x;
        yPoints[nPoints] = y;
        nPoints++;
        if (bounds != null) {
            updateBounds(x, y);
        }
    }


    public void clear() {
        nPoints = 0;
    }


    /**
     * Gets the bounding box of this <code>Polygon2D</code>.
     * The bounding box is the smallest {@link Rectangle} whose
     * sides are parallel to the x and y axes of the
     * coordinate space, and can completely contain the <code>Polygon2D</code>.
     * @return a <code>Rectangle</code> that defines the bounds of this
     * <code>Polygon2D</code>.
     */
    public Rectangle getBounds() {
        if (bounds == null) {
            calculateBounds(xPoints, yPoints, nPoints);
        }
        return new Rectangle( (int) bounds.x, (int) bounds.y, (int) bounds.width,
                             (int) bounds.height);
    }


    /**
     * Determines whether the specified {@link Point} is inside this
     * <code>Polygon2D</code>.
     * @param p the specified <code>Point</code> to be tested
     * @return <code>true</code> if the <code>Polygon2D</code> contains the
     * 			<code>Point</code>; <code>false</code> otherwise.
     * @see #contains(double, double)
     */
    public boolean contains(Point p) {
        return contains(p.x, p.y);
    }


    /**
     * Determines whether the specified coordinates are inside this
     * <code>Polygon2D</code>.
     * <p>
     * @param x  the specified coordinates to be tested
     * @param y  the specified coordinates to be tested
     * @return  <code>true</code> if this <code>Polygon2D</code> contains
     * 			the specified coordinates, (<i>x</i>,&nbsp;<i>y</i>);
     * 			<code>false</code> otherwise.
     */
    public boolean contains(int x, int y) {
        return contains( (double) x, (double) y);
    }


    /**
     * Returns the high precision bounding box of the {@link Shape}.
     * @return a {@link Rectangle2D} that precisely
     *		bounds the <code>Shape</code>.
     */
    public Rectangle2D getBounds2D() {
//        if (bounds == null) {
        calculateBounds(xPoints, yPoints, nPoints);
//        }
        return new Rectangle2D.Double(bounds.x, bounds.y, bounds.width, bounds.height);
    }


    public boolean contains(double x, double y) {
        boolean contained = false;
        int i, j;

        for (i = 0, j = nPoints - 1; i < nPoints; j = i++) {
            if ( ( ( (yPoints[i] <= y) && (y < yPoints[j])) ||
                  ( (yPoints[j] <= y) && (y < yPoints[i]))) &&
                (x <
                 (xPoints[j] - xPoints[i]) * (y - yPoints[i]) /
                 (yPoints[j] - yPoints[i]) + xPoints[i])) {
                contained = !contained;
            }
        }

        return contained;
    }


    private Crossings getCrossings(double xlo, double ylo, double xhi, double yhi) {
        Crossings cross = new Crossings.EvenOdd(xlo, ylo, xhi, yhi);
        double lastx = xPoints[nPoints - 1];
        double lasty = yPoints[nPoints - 1];
        double curx, cury;

        // Walk the edges of the polygon
        for (int i = 0; i < nPoints; i++) {
            curx = xPoints[i];
            cury = yPoints[i];
            if (cross.accumulateLine(lastx, lasty, curx, cury)) {
                return null;
            }
            lastx = curx;
            lasty = cury;
        }

        return cross;
    }


    /**
     * Tests if a specified {@link Point2D} is inside the boundary of this
     * <code>Polygon2D</code>.
     * @param p a specified <code>Point2D</code>
     * @return <code>true</code> if this <code>Polygon2D</code> contains the
     * 		specified <code>Point2D</code>; <code>false</code> otherwise.
     */
    public boolean contains(Point2D p) {
        return contains(p.getX(), p.getY());
    }


    /**
     * Tests if the interior of this <code>Polygon2D</code> intersects the
     * interior of a specified set of rectangular coordinates.
     * @param x the x coordinate of the specified rectangular shape's top-left corner
     * @param y the y coordinate of the specified rectangular shape's top-left corner
     * @param w the width of the specified rectangular shape
     * @param h the height of the specified rectangular shape
     * @return <code>true</code> if the interior of this
     *		<code>Polygon2D</code> and the interior of the specified set of
     * rectangular coordinates intersect each other; <code>false</code> otherwise.
     */
    public boolean intersects(double x, double y, double w, double h) {
        if (nPoints <= 0 || !getBounds().intersects(x, y, w, h)) {
            return false;
        }

        Crossings cross = getCrossings(x, y, x + w, y + h);
        return (cross == null || !cross.isEmpty());
    }


    /**
     * Tests if the interior of this <code>Polygon2D</code> intersects the
     * interior of a specified <code>Rectangle2D</code>.
     * @param r a specified <code>Rectangle2D</code>
     * @return <code>true</code> if this <code>Polygon2D</code> and the
     * 			interior of the specified <code>Rectangle2D</code>
     * 			intersect each other; <code>false</code>
     * 			otherwise.
     */
    public boolean intersects(Rectangle2D r) {
        return intersects(r.getX(), r.getY(), r.getWidth(), r.getHeight());
    }


    /**
     * Tests if the interior of this <code>Polygon2D</code> entirely
     * contains the specified set of rectangular coordinates.
     * @param x the x coordinat of the top-left corner of the
     * 			specified set of rectangular coordinates
     * @param y the y coordinat of the top-left corner of the
     * 			specified set of rectangular coordinates
     * @param w the width of the set of rectangular coordinates
     * @param h the height of the set of rectangular coordinates
     * @return <code>true</code> if this <code>Polygon2D</code> entirely
     * 			contains the specified set of rectangular
     * 			coordinates; <code>false</code> otherwise.
     */
    public boolean contains(double x, double y, double w, double h) {
        if (nPoints <= 0 || !getBounds().intersects(x, y, w, h)) {
            return false;
        }

        Crossings cross = getCrossings(x, y, x + w, y + h);
        return (cross != null && cross.covers(y, y + h));
    }


    /**
     * Tests if the interior of this <code>Polygon2D</code> entirely
     * contains the specified <code>Rectangle2D</code>.
     * @param r the specified <code>Rectangle2D</code>
     * @return <code>true</code> if this <code>Polygon2D</code> entirely
     * 			contains the specified <code>Rectangle2D</code>;
     *			<code>false</code> otherwise.
     */
    public boolean contains(Rectangle2D r) {
        return contains(r.getX(), r.getY(), r.getWidth(), r.getHeight());
    }


    /**
     * Returns an iterator object that iterates along the boundary of this
     * <code>Polygon2D</code> and provides access to the geometry
     * of the outline of this <code>Polygon2D</code>.  An optional
     * {@link AffineTransform} can be specified so that the coordinates
     * returned in the iteration are transformed accordingly.
     * @param at an optional <code>AffineTransform</code> to be applied to the
     * 		coordinates as they are returned in the iteration, or
     *		<code>null</code> if untransformed coordinates are desired
     * @return a {@link PathIterator} object that provides access to the
     *		geometry of this <code>Polygon2D</code>.
     */
    public PathIterator getPathIterator(AffineTransform at) {
        return new Polygon2DPathIterator(this, at);
    }


    /**
     * Returns an iterator object that iterates along the boundary of
     * the <code>Shape</code> and provides access to the geometry of the
     * outline of the <code>Shape</code>.  Only SEG_MOVETO, SEG_LINETO, and
     * SEG_CLOSE point types are returned by the iterator.
     * Since polygons are already flat, the <code>flatness</code> parameter
     * is ignored.  An optional <code>AffineTransform</code> can be specified
     * in which case the coordinates returned in the iteration are transformed
     * accordingly.
     * @param at an optional <code>AffineTransform</code> to be applied to the
     * 		coordinates as they are returned in the iteration, or
     *		<code>null</code> if untransformed coordinates are desired
     * @param flatness the maximum amount that the control points
     * 		for a given curve can vary from colinear before a subdivided
     *		curve is replaced by a straight line connecting the
     * 		endpoints.  Since polygons are already flat the
     * 		<code>flatness</code> parameter is ignored.
     * @return a <code>PathIterator</code> object that provides access to the
     * 		<code>Shape</code> object's geometry.
     */
    public PathIterator getPathIterator(AffineTransform at, double flatness) {
        return getPathIterator(at);
    }


    class Polygon2DPathIterator
        implements PathIterator {
        Polygon2D poly;
        AffineTransform transform;
        int index = 0;

        public Polygon2DPathIterator(Polygon2D pg, AffineTransform at) {
            poly = pg;
            transform = at;
        }


        /**
         * Returns the winding rule for determining the interior of the
         * path.
         * @return an integer representing the current winding rule.
         * @see PathIterator#WIND_NON_ZERO
         */
        public int getWindingRule() {
            return WIND_EVEN_ODD;
        }


        /**
         * Tests if there are more points to read.
         * @return <code>true</code> if there are more points to read;
         *          <code>false</code> otherwise.
         */
        public boolean isDone() {
            return index > poly.nPoints;
        }


        /**
         * Moves the iterator forwards, along the primary direction of
         * traversal, to the next segment of the path when there are
         * more points in that direction.
         */
        public void next() {
            index++;
        }


        /**
         * Returns the coordinates and type of the current path segment in
         * the iteration.
         * The return value is the path segment type:
         * SEG_MOVETO, SEG_LINETO, or SEG_CLOSE.
         * A <code>float</code> array of length 2 must be passed in and
         * can be used to store the coordinates of the point(s).
         * Each point is stored as a pair of <code>float</code> x,&nbsp;y
         * coordinates.  SEG_MOVETO and SEG_LINETO types return one
         * point, and SEG_CLOSE does not return any points.
         * @param coords a <code>float</code> array that specifies the
         * coordinates of the point(s)
         * @return an integer representing the type and coordinates of the
         * 		current path segment.
         * @see PathIterator#SEG_MOVETO
         * @see PathIterator#SEG_LINETO
         * @see PathIterator#SEG_CLOSE
         */
        public int currentSegment(float[] coords) {
            if (index >= poly.nPoints) {
                if (close) {
                    return SEG_CLOSE;
                } else {
                    return SEG_MOVETO;
                }
            }
            coords[0] = (float) poly.xPoints[index];
            coords[1] = (float) poly.yPoints[index];
            if (transform != null) {
                transform.transform(coords, 0, coords, 0, 1);
            }
            return (index == 0 ? SEG_MOVETO : SEG_LINETO);
        }


        /**
         * Returns the coordinates and type of the current path segment in
         * the iteration.
         * The return value is the path segment type:
         * SEG_MOVETO, SEG_LINETO, or SEG_CLOSE.
         * A <code>double</code> array of length 2 must be passed in and
         * can be used to store the coordinates of the point(s).
         * Each point is stored as a pair of <code>double</code> x,&nbsp;y
         * coordinates.
         * SEG_MOVETO and SEG_LINETO types return one point,
         * and SEG_CLOSE does not return any points.
         * @param coords a <code>double</code> array that specifies the
         * coordinates of the point(s)
         * @return an integer representing the type and coordinates of the
         * 		current path segment.
         * @see PathIterator#SEG_MOVETO
         * @see PathIterator#SEG_LINETO
         * @see PathIterator#SEG_CLOSE
         */
        public int currentSegment(double[] coords) {
            if (index >= poly.nPoints) {
                return SEG_CLOSE;
            }
            coords[0] = poly.xPoints[index];
            coords[1] = poly.yPoints[index];
            if (transform != null) {
                transform.transform(coords, 0, coords, 0, 1);
            }
            return (index == 0 ? SEG_MOVETO : SEG_LINETO);
        }
    }


    //////////////////////////////////////////////////////////////////////////////////////

    //Only works for a single simply-connected polygon, returns a negative area for counter-clockwise,
    //positive for clockwise. Multiple polygons can be added together if they do not intersect,
    //external polygons have CW chirality and internal polygons have CCW chirality.
    public double getArea() {
        double a = 0;
        for (int i = 0; i < nPoints - 1; i++) {
            a += -xPoints[i] * yPoints[i + 1] + xPoints[i + 1] * yPoints[i];
        }
        return a / 2;
    }


    public double[] getCentroid() {
        double x = 0, y = 0, a = getArea();
        for (int i = 0; i < nPoints - 1; i++) {
            x += (xPoints[i] + xPoints[i + 1]) * (xPoints[i + 1]*yPoints[i] - xPoints[i]*yPoints[i + 1]);
            y += (yPoints[i] + yPoints[i + 1]) * (xPoints[i + 1]*yPoints[i] - xPoints[i]*yPoints[i + 1]);
        }
        return new double[] {x / (6.0 * a), y / (6.0 * a)};
    }


    //Changes CCW to CW and vice versa.
    //When doing area calculations, changes sign of area.
    public void switchChirality() {

        Polygon2D newPoly = new Polygon2D();
        for (int i = 0; i < nPoints; i++) {
            newPoly.addPoint(xPoints[nPoints - 1 - i], yPoints[nPoints - 1 - i]);
        }
        clear();
        for (int i = 0; i < newPoly.nPoints; i++) {
            addPoint(newPoly.xPoints[i], newPoly.yPoints[i]);
        }

    }


    /**
     * Takes multiple contours that have been put into a single polygon, terminated
     * by their initial points, and breaks it into multiple polygons.  Last point in
     * the incoming list is the same as the last point in the incoming list.
     * External polygons are given CW chirality.  Internal polygons are given CCW
     * chirality.
     * @return Polygon2D
     */
    public Polygon2D[] divideIntoSeparatePolygons() {
        ArrayList<Polygon2D> polygons = new ArrayList<Polygon2D> ();

        Polygon2D polygon_current = new Polygon2D();

        double xFirst = 0;
        double yFirst = 0;
        boolean first = true;
//        System.out.println("big break");
        //Break polygon into sub-polygons.
        for (int ipoint = 0; ipoint < nPoints - 1; ipoint++) {
            if (first) {
                xFirst = xPoints[ipoint];
                yFirst = yPoints[ipoint];
            }
            polygon_current.addPoint(xPoints[ipoint],
                                     yPoints[ipoint]);
//            System.out.println("x: " + xPoints[ipoint]);
//            System.out.println("y: " + yPoints[ipoint]);
            if (xPoints[ipoint] == xFirst &&
                yPoints[ipoint] == yFirst && !first) {
//                System.out.println("small break");
//                System.out.println("Area: " + polygon_current.getArea());
                polygons.add(polygon_current);
                polygon_current = new Polygon2D();
                first = true;

            } else {
                first = false;

            }
        }

        //Switch all polygons to CW chirality
        for (int i = 0; i < polygons.size(); i++) {
            polygon_current = polygons.get(i);
            if (polygon_current.getArea() < 0) {
                polygon_current.switchChirality();
            }
        }
        //Change the chirality of every contained polygon once for every
        //time that it is contained.
        //This makes it so that donut holes have negative area
        // Note that this works for interlocking donuts and donuts with
        // multiple holes as well.
        for (int i = 0; i < polygons.size(); i++) {
            for (int j = i + 1; j < polygons.size(); j++) {
                if (polygons.get(i).contains(
                    polygons.get(j).xPoints[0], polygons.get(j).yPoints[0])) {
                    polygons.get(j).switchChirality();

                }

                else if (polygons.get(j).contains(
                    polygons.get(i).xPoints[0], polygons.get(i).yPoints[0])) {
                    polygons.get(i).switchChirality();

                }

            }

        }
        Polygon2D[] polygonsArray = new Polygon2D[polygons.size()];
        for (int i = 0; i < polygons.size(); i++) {
            polygonsArray[i] = polygons.get(i);
        }

        return polygonsArray;

    }


    /**
     * Creates polygon envelope that is guaranteed to be convex.
     * @return Polygon2D
     */
    public Polygon2D simplifyPolygon() {
        if(nPoints>2) {
        Polygon2D simple = new Polygon2D();

        //Find highest point. Guaranteed to be convex.
        double max = Double.NEGATIVE_INFINITY;
        int maxIndex = -1;
        for (int i = 0; i < nPoints; i++) {
            if (yPoints[i] > max) {
                max = yPoints[i];
                maxIndex = i;
            }
        }
        simple.addPoint(xPoints[maxIndex], yPoints[maxIndex]);

        boolean[] pointDone = new boolean[nPoints];

        int currentIndex = maxIndex;
        double[] currentVector = new double[] {1.0, 0};
        double[] newVector = new double[] {0.0, 0.0};
        double[] nextVector = {0.0, 0.0};
        int nextPoint = -1;

        do {
            double maxDot = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < nPoints; i++) {
                if (!pointDone[i] && i != currentIndex) {
                    newVector = new double[] {xPoints[i] - xPoints[currentIndex]
                                , yPoints[i] - yPoints[currentIndex]};
                    double newVectorMagnitude = Math.sqrt(newVector[0] * newVector[0] +
                        newVector[1] * newVector[1]);
                    newVector[0] = newVector[0] / newVectorMagnitude;
                    newVector[1] = newVector[1] / newVectorMagnitude;
                    double dot = currentVector[0] * newVector[0] +
                                 currentVector[1] * newVector[1];
                    if (dot > maxDot) {
                        nextPoint = i;
                        maxDot = dot;
                        nextVector = newVector;
                    }
                }
            }
            pointDone[nextPoint] = true;
            simple.addPoint(xPoints[nextPoint], yPoints[nextPoint]);
            currentIndex = nextPoint;
            currentVector = nextVector;

            //while you have not made it back to the starting point
        } while (! (simple.xPoints[0] == simple.xPoints[simple.nPoints - 1] &&
                    simple.yPoints[0] == simple.yPoints[simple.nPoints - 1]));
//        MathUtil.resize(simple.xPoints, simple.nPoints);
//        MathUtil.resize(simple.yPoints, simple.nPoints);
        return new Polygon2D(simple.xPoints, simple.yPoints, simple.nPoints);
        } else {
            return new Polygon2D(xPoints, yPoints, nPoints);
        }
    }


    private boolean intersect(double x1, double y1, double x2, double y2,
                              double x3, double y3, double x4, double y4) {
        double f = (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1);
        if (f == 0) {
            return false;
        }
        double ua = ( (x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / f;
        double ub = ( (x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / f;
        if (ua > 0 && ua < 1 && ub > 0 && ub < 1) {
            return true;
        }

        return false;
    }


    public boolean isSelfIntersecting() {
        for (int i = 0; i < nPoints - 1; i++) {
            for (int j = 0; j < nPoints - 1; j++) {
                if (intersect(xPoints[i], yPoints[i], xPoints[i + 1], yPoints[i + 1],
                              xPoints[j], yPoints[j], xPoints[j + 1], yPoints[j + 1])) {
                    return true;
                }
            }
        }

        return false;
    }


    public double getMinimumDistance(double x, double y) {
        double minD = Double.POSITIVE_INFINITY;
        for (int i = 0; i < nPoints - 1; i++) {
            double d = MathUtil.getPointToLineDistance(
                x, y, xPoints[i], yPoints[i], xPoints[i + 1], yPoints[i + 1], true);
            if (d < minD) {
                minD = d;
            }
        }
        double d = MathUtil.getPointToLineDistance(
            x, y, xPoints[nPoints - 1], yPoints[nPoints - 1],
            xPoints[0], yPoints[0], true);
        if (d < minD) {
            minD = d;
        }
        return minD;
    }
    
    public void transform(double xScaleFactor, double yScaleFactor, 
            double xPositionIncrement, double yPositionIncrement) { 
        MathUtil.multiply(xPoints, xScaleFactor);
        MathUtil.multiply(yPoints, yScaleFactor);
        MathUtil.add(xPoints, xPositionIncrement);
        MathUtil.add(yPoints, yPositionIncrement);

        // FIXME: Need to check call hierarchy; people may have gotten used to this not happening.
        // bounds = null;
    }
    
    public void scale(double xFactor, double yFactor) {
        double[] centroid = getCentroid();
        translate(-centroid[0], -centroid[1]);
        transform(xFactor, yFactor, centroid[0], centroid[1]);

        // FIXME: Really this should be done in transform, but need to check call hierarchy...
        bounds = null;
    }
    
    public void grow(double xIncrement, double yIncrement) {
        Rectangle2D bounds = getBounds2D();
        double xFactor = (bounds.getWidth()  + 2*xIncrement) / bounds.getWidth();
        double yFactor = (bounds.getHeight() + 2*yIncrement) / bounds.getHeight();
        scale(xFactor, yFactor);
    }
    
    public void grow(double increment) { grow(increment, increment); }
    
    /**
     * Designed to match ElectrodeMap#flipAndRotate
     * 
     * @param flipX
     * @param flipY
     * @param angle
     * @param centerX
     * @param centerY
     * 
     * NOTE: If you make any change here, be sure to make matching change to ElectrodeMap#flipAndRotate, or 519
     * streaming will break.
     * 
     * FIXME: Should actually use THE SAME code to do both, so make an abstracted call or encapsulate the x and y points
     * in both objects into a Coordinates object that has this method...
     */
    public void flipAndRotate(boolean flipX, boolean flipY, double angle, double centerX, double centerY) {
        for (int i = 0; i < nPoints; i++) {
            xPoints[i] = flipX ? - xPoints[i] : xPoints[i];
            yPoints[i] = flipY ? - yPoints[i] : yPoints[i];
        }

        double newX, newY;
        for (int i = 0; i < nPoints; i++) {
            //modifies variables in super.  Local xCoord, yCoord are final "pure" definition.
            newX = xPoints[i] * Math.cos(angle) - yPoints[i] * Math.sin(angle) + centerX;
            newY = xPoints[i] * Math.sin(angle) + yPoints[i] * Math.cos(angle) + centerY;
            xPoints[i] = (float) newX;
            yPoints[i] = (float) newY;
        }		 

    }
}
