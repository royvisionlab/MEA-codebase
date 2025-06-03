package edu.ucsc.neurobiology.vision.plot;

import java.awt.geom.*;


/**
 * Used to define and draw poligons.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public interface PolygonData
    extends PlotData {

    public void reset();


    public boolean hasNextPoint();


    public Point2D nextPoint();
}
