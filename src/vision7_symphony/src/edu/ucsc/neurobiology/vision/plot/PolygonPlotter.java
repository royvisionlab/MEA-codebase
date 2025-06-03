package edu.ucsc.neurobiology.vision.plot;

import java.awt.*;
import java.awt.geom.*;
import java.awt.image.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class PolygonPlotter
    implements DataPlotter {

    public synchronized boolean accept(Object data, Object style) {
        if ( (data instanceof PolygonData)) {
            return true;
        } else {
            return false;
        }
    }


    public synchronized void draw(
        Axis horAxis, Axis verAxis, Graphics g, Raster raster, Object _data,
        Object _style) {

        FunctionStyle style = (FunctionStyle) _style;
        PolygonData data = (PolygonData) _data;

        g.setColor(style.getLineColor());
        ( (Graphics2D) g).setStroke(new BasicStroke(style.getLineThickness(),
            BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
        data.reset();
        int xOld = 0, yOld = 0;
        boolean first = true;

        while (data.hasNextPoint()) {
            Point2D p = data.nextPoint();
            int x = (int) Math.round(horAxis.getScreenCoord(p.getX()));
            int y = (int) Math.round(verAxis.getScreenCoord(p.getY()));

            if (first) {
                xOld = x;
                yOld = y;
                first = false;
            }

            g.drawLine(xOld, yOld, x, y);

            xOld = x;
            yOld = y;
        }
    }

}
