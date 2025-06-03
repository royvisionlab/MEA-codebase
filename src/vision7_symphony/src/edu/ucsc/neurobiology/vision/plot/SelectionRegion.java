package edu.ucsc.neurobiology.vision.plot;

import java.awt.geom.*;


/**
 * Defines a selection region in the PlotPanel.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public interface SelectionRegion {
    /**
     * Does this selection contain the given point?
     *
     * @param x double
     * @param y double
     * @return boolean
     */
    public abstract boolean contains(double x, double y);


    /**
     * What is the smallest rectangle that bounds the region?
     *
     * @return Rectangle2D
     */
    public abstract Rectangle2D getBounds2D();


    /**
     * What is the minimum distance from this point to any point on the region?
     *
     * @param x double
     * @param y double
     * @return double
     */
    public abstract double getMinimumDistance(double x, double y);

}
