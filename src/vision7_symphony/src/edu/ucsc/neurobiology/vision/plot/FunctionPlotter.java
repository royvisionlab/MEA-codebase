package edu.ucsc.neurobiology.vision.plot;

import java.awt.*;
import java.awt.image.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class FunctionPlotter implements DataPlotter {

    Polygon2D polygon1 = new Polygon2D(false);
    Polygon2D polygon2 = new Polygon2D(true);
    
    
    public synchronized boolean accept(Object data, Object style) {
        if ( (style instanceof FunctionStyle) &&
            (data instanceof FunctionData || data instanceof ParametricFunctionData)) {
            return true;
        } else {
            return false;
        }
    }


    public synchronized void draw(Axis horAxis, Axis verAxis, Graphics g, Raster raster,
                                  Object data,
                                  Object style) {
        try {
            if (data instanceof FunctionData) {
                drawSPI(horAxis, verAxis, (Graphics2D) g, (FunctionData) data,
                        (FunctionStyle) style);
            } else if (data instanceof ParametricFunctionData) {
                drawSPI(horAxis, verAxis, (Graphics2D) g, (ParametricFunctionData) data,
                        (FunctionStyle) style);
            }
        } catch (CannotEvaluateException e) {
            Vision.reportException("", e);
        }
    }


    private synchronized void drawSPI(Axis horAxis, Axis verAxis,
                                      Graphics2D g, FunctionData data,
                                      FunctionStyle style) throws CannotEvaluateException {
        polygon1.clear();
        int x1 = Math.max(horAxis.s1, (int) horAxis.getScreenCoord(style.getMinimum()));
        int x2 = Math.min(horAxis.s2, (int) horAxis.getScreenCoord(style.getMaximum()));
        for (double x = x1; x < x2; x += .1) {
            double y = verAxis.getScreenCoord(data.getValueAt(horAxis.getPlotCoord(x)));
            polygon1.addPoint(x, y);
        }

        g.setColor(style.getLineColor());
        String dash = style.getDashPattern();
        if (dash == null || dash.length() == 0) {
            g.setStroke(new BasicStroke(style.getLineThickness(), BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 1));
        } else {
            String[] s = StringUtil.decomposeString(dash, ",");
            float[] d = new float[s.length];
            for (int i = 0; i < s.length; i++) {
                d[i] = Float.parseFloat(s[i]);
            }
            g.setStroke(new BasicStroke(
                style.getLineThickness(), BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 1,
                d, 0));
        }
        g.draw(polygon1);
    }


    private synchronized void drawSPI(Axis horAxis, Axis verAxis,
                                      Graphics2D g, ParametricFunctionData data,
                                      FunctionStyle style) {
        polygon2.clear();
        g.setStroke(new BasicStroke(style.getLineThickness()));

        int n = 50;
        double dp = (data.getMaxParamValue() - data.getMinParamValue()) / n;
        double p = data.getMinParamValue();

        for (int i = 0; i < n + n / 20; i++, p += dp) {
            double[] point = data.getPointFor(p);
            polygon2.addPoint(horAxis.getScreenCoord(point[0]),
                              verAxis.getScreenCoord(point[1]));
        }
        g.setColor(style.getLineColor());
        g.draw(polygon2);
        
        String label = data.getLabel();
        if (label != null) {
            g.setFont(style.getLabelFont());
            double[] point = data.getLabelPosition();
            g.drawString(label, (int) horAxis.getScreenCoord(point[0]), (int) verAxis.getScreenCoord(point[1]));
        }
    }

}
