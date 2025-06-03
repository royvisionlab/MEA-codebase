package edu.ucsc.neurobiology.vision.plot;

import java.io.*;

import java.awt.*;
import java.awt.image.*;


/**
 * This interface is implemented by a concrete class to allow it to act as a plotter
 * for plots of particular kinds. The plotters have to be registered in the
 * <tt>PlotPanel</tt> before the data gets added.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public interface DataPlotter {
    /**
     * Returns true if the given data and style can be handled by this
     * <tt>DataPlotter</tt> and false otherwise.
     */
    public boolean accept(Object data, Object style);


    /**
     * Draws the given data using the given style on the graphics context g.
     * Translations from real space to display space can be done by using the given
     * axes.
     */
    public void draw(Axis horAxis, Axis verAxis, Graphics g, Raster raster, Object data, Object style);
}
