package edu.ucsc.neurobiology.vision.plot;

import java.awt.*;
import java.awt.image.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class Histogram2DPlotter
    implements DataPlotter {

    public synchronized boolean accept(Object data, Object style) {
        if ( (data instanceof Histogram2DData) && (style instanceof HistogramStyle)) {
            return true;
        } else {
            return false;
        }
    }


    public synchronized void draw(
        Axis horAxis, Axis verAxis, Graphics graphics, Raster raster, Object _data,
        Object _style) {

        Histogram2DData data = (Histogram2DData) _data;
        HistogramStyle style = (HistogramStyle) _style;

//        double dx = (data.getMaxX() - data.getMinX())/data.getBinCountX();
//        double dy = (data.getMaxY() - data.getMinY())/data.getBinCountY();
        double dx = data.getBinIntervalX();
        double dy = data.getBinIntervalY();

        double max = Double.NEGATIVE_INFINITY;
        double min = Double.POSITIVE_INFINITY;
        for (int i = 0; i < data.getBinCountX(); i++) {
            for (int j = 0; j < data.getBinCountY(); j++) {
                double v = data.getBin(i, j);
                if (v > max) {
                    max = v;
                }
                if (v < min) {
                    min = v;
                }
            }
        }

        for (int i = 0; i < data.getBinCountX(); i++) {
            for (int j = 0; j < data.getBinCountY(); j++) {
                int x1 = (int) Math.round(horAxis.getScreenCoord(data.getMinX() +
                    i * dx));
                int y1 = (int) Math.round(verAxis.getScreenCoord(data.getMinY() +
                    j * dy));
                int x2 = (int) Math.round(horAxis.getScreenCoord(data.getMinX() +
                    (i + 1) * dx));
                int y2 = (int) Math.round(verAxis.getScreenCoord(data.getMinY() +
                    (j + 1) * dy));

                float c = (float) ( (data.getBin(i, j) - min) / (max - min));
                graphics.setColor(new Color(c, c, c));
                graphics.fillRect(Math.min(x1, x2), Math.min(y1, y2),
                                  Math.abs(x2 - x1), Math.abs(y2 - y1));
            }
        }
    }

}
