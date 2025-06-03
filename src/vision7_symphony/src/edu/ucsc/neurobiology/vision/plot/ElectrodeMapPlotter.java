package edu.ucsc.neurobiology.vision.plot;

import java.awt.*;
import java.awt.geom.*;
import java.awt.image.*;

import edu.ucsc.neurobiology.vision.electrodemap.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ElectrodeMapPlotter
    implements DataPlotter {
    
    public synchronized boolean accept(Object data, Object style) {
        if ( (data instanceof ElectrodeMap) && (style instanceof ElectrodeMapStyle)) {
            return true;
        } else {
            return false;
        }
    }
    
    public synchronized void draw(Axis horAxis, Axis verAxis, Graphics _g, Raster raster,
                                  Object data, Object style) {
        
        Graphics2D g = (Graphics2D) _g;
        ElectrodeMap electrodeMap = (ElectrodeMap) data;
        double w = 4;
        g.setColor(Color.black);
        
        Point2D.Double p = null;
        int n = electrodeMap.getNumberOfElectrodes();
        Rectangle2D.Double r = new Rectangle2D.Double();
        for (int i = 0; i < n; i++) {
            p = (Point2D.Double) electrodeMap.getPosition(i, p);
            p.x = horAxis.getScreenCoord(p.x);
            p.y = verAxis.getScreenCoord(p.y);

            r.setFrame(p.x - w / 2, p.y - w / 2, w, w);
            g.fill(r);
            
            String label = Integer.toString(i);
            g.drawString(label, (int) p.x + 2, (int) p.y - 2);
        }
    }
}
