package edu.ucsc.neurobiology.vision.plot;

import java.awt.geom.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class Polygon2DAdapter
    implements PolygonData {
    int i;
    Point2D.Double p = new Point2D.Double(0, 0);
    Polygon2D polygon;


    public void setPolygon(Polygon2D polygon) {
        this.polygon = polygon;
    }


    public void reset() {
        i = 0;
    }


    public boolean hasNextPoint() {
        if (polygon == null) {
            return false;
        }
        return (i < polygon.nPoints);
    }


    public Point2D nextPoint() {
        p.setLocation(polygon.xPoints[i], polygon.yPoints[i]);
        i++;
        return p;
    }


    public String getDescription() {
        return "pd";
    }
}
