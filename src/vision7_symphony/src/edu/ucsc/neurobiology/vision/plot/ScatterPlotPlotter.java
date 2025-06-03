package edu.ucsc.neurobiology.vision.plot;

import java.awt.*;
import java.awt.geom.*;
import java.awt.image.*;

import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ScatterPlotPlotter
    implements DataPlotter {

    Polygon2D p = new Polygon2D(false);


    public synchronized boolean accept(Object data, Object style) {
        if ( (data instanceof ScatterPlotData) && (style instanceof ScatterPlotStyle)) {
            return true;
        } else {
            return false;
        }
    }


    public synchronized void draw(
        Axis horAxis, Axis verAxis, Graphics _g, Raster raster, Object _data,
        Object _style) {

        Graphics2D g = (Graphics2D) _g;
        Ellipse2D.Double ellipse = new Ellipse2D.Double();
        Line2D.Double line = new Line2D.Double();
        Rectangle2D.Double rect = new Rectangle2D.Double();

        ScatterPlotData data = (ScatterPlotData) _data;
        ScatterPlotStyle style = (ScatterPlotStyle) _style;

        int nPoints = data.getPointCount();
        double[] point = new double[4];

        final int s = style.getSymbolSize();
        final double ss = s / 2.0;

        final double x0 = horAxis.getScreenCoord(0);
        final double y0 = verAxis.getScreenCoord(0);

        // draw connecting lines
        if (style.getConnectingPoints()) {
            p.clear();
            int connectionPeriod = style.getConnectionPeriod();

            g.setColor(style.getConnectionLineColor());

            // set stroke
//            g.setStroke(new BasicStroke(style.getConnectionLineThickness()));
            String dash = style.getDashPattern();
            if (dash == null || dash.length() == 0) {
                g.setStroke(new BasicStroke(
                    style.getConnectionLineThickness(), BasicStroke.CAP_BUTT,
                    BasicStroke.JOIN_MITER, 1));
            } else {
                String[] str = StringUtil.decomposeString(dash, ",");
                float[] d = new float[str.length];
                for (int i = 0; i < str.length; i++) {
                    d[i] = Float.parseFloat(str[i]);
                }
                g.setStroke(new BasicStroke(
                    style.getConnectionLineThickness(), BasicStroke.CAP_BUTT,
                    BasicStroke.JOIN_MITER, 1, d, 0));
            }

            for (int i = 0; i < nPoints; i++) {
                data.getDataPoint(i, point);
                double x = horAxis.getScreenCoord(point[0]);
                double y = verAxis.getScreenCoord(point[1]);

//                if (i % connectionPeriod != 0) {
                p.addPoint(x, y);
//                }

//                if (i != 0 && i % connectionPeriod != 0) {
//                    g.drawLine(xOld, yOld, x, y);
//                }
            }
            g.draw(p);
        }

        int errorTick = 1;

        // draw the error bars
        if (data.hasXErrors()) {
            g.setColor(style.getSymbolColor());
            g.setStroke(new BasicStroke(style.getConnectionLineThickness()));

            switch (style.getErrorSymbolType()) {
                case VERTICAL_LINE:
                    for (int i = 0; i < nPoints; i++) {
                        data.getDataPoint(i, point);
                        if (point[3] != 0) {
                            double xLeft = horAxis.getScreenCoord(point[0] - point[3]);
                            double xRight = horAxis.getScreenCoord(point[0] + point[3]);

                            double y = verAxis.getScreenCoord(point[1]);
                            line.setLine(xLeft, y, xRight, y);
                            g.draw(line);

                            line.setLine(xLeft, y - errorTick, xLeft, y + errorTick);
                            g.draw(line);

                            line.setLine(xRight, y - errorTick, xRight, y + errorTick);
                            g.draw(line);
                        }
                    }
                    break;
                    
                case NONE:
                    for (int i = 0; i < nPoints; i++) {
                        data.getDataPoint(i, point);
                        if (point[3] != 0) {
                            double xLeft = horAxis.getScreenCoord(point[0] - point[3]);
                            double xRight = horAxis.getScreenCoord(point[0] + point[3]);

                            double y = verAxis.getScreenCoord(point[1]);
                            line.setLine(xLeft, y, xRight, y);
                            g.draw(line);

                        }
                    }
                    break;
            }
        }

        // draw the error bars
        if (data.hasYErrors()) {
            g.setColor(style.getSymbolColor());
            g.setStroke(new BasicStroke(style.getConnectionLineThickness()));

            switch (style.getErrorSymbolType()) {
                case VERTICAL_LINE:
                    for (int i = 0; i < nPoints; i++) {
                        data.getDataPoint(i, point);
                        if (point[2] != 0) {
                            double yUp = verAxis.getScreenCoord(point[1] + point[2]);
                            double yDown = verAxis.getScreenCoord(point[1] - point[2]);
                            double x = horAxis.getScreenCoord(point[0]);

                            line.setLine(x, yUp, x, yDown);
                            g.draw(line);

                            line.setLine(x - errorTick, yUp, x + errorTick, yUp);
                            g.draw(line);

                            line.setLine(x - errorTick, yDown, x + errorTick, yDown);
                            g.draw(line);
                        }
                    }
                    break;

                case DISK:
                    for (int i = 0; i < nPoints; i++) {
                        data.getDataPoint(i, point);
                        double x = horAxis.getScreenCoord(point[0]);
                        double y = verAxis.getScreenCoord(point[1]);
                        double rx = 3 * Math.abs(horAxis.getScreenCoord(point[2]) - x0);
                        double ry = 3 * Math.abs(verAxis.getScreenCoord(point[2]) - y0);
                        ellipse.setFrame(x - rx / 2, y - ry / 2, rx, ry);
                        g.fill(ellipse);
                    }
                    break;

                case CIRCLE:
                    for (int i = 0; i < nPoints; i++) {
                        data.getDataPoint(i, point);
                        double x = horAxis.getScreenCoord(point[0]);
                        double y = verAxis.getScreenCoord(point[1]);
                        double rx = 3 * Math.abs(horAxis.getScreenCoord(point[2]) - x0);
                        double ry = 3 * Math.abs(verAxis.getScreenCoord(point[2]) - y0);
                        ellipse.setFrame(x - rx / 2, y - ry / 2, rx, ry);
                        g.draw(ellipse);
                    }
                    break;
                    
                case NONE:
                    for (int i = 0; i < nPoints; i++) {
                        data.getDataPoint(i, point);
                        if (point[2] != 0) {
                            double yUp = verAxis.getScreenCoord(point[1] + point[2]);
                            double yDown = verAxis.getScreenCoord(point[1] - point[2]);
                            double x = horAxis.getScreenCoord(point[0]);

                            line.setLine(x, yUp, x, yDown);
                            g.draw(line);

                        }
                    }
                    break;
            }
        }

        int color = style.getSymbolColor().getRGB();
        final SymbolType sType = style.getSymbolType();
        g.setColor(style.getSymbolColor());
        g.setStroke(new BasicStroke(1f / 2f));

        switch (sType) {
            case FILLED_SQUARE:
                for (int i = 0; i < nPoints; i++) {
                    data.getDataPoint(i, point);
                    double x = horAxis.getScreenCoord(point[0]);
                    double y = verAxis.getScreenCoord(point[1]);
                    rect.setRect(x - ss, y - ss, s, s);
                    g.fill(rect);
                }
                break;

            case SQUARE:
                if (raster == null) {
                    for (int i = 0; i < nPoints; i++) {
                        data.getDataPoint(i, point);
                        double x = horAxis.getScreenCoord(point[0]);
                        double y = verAxis.getScreenCoord(point[1]);
                        rect.setRect(x - ss, y - ss, s, s);
                        g.draw(rect);
                    }
                } else {
                    DataBufferInt dataBufferByte = (DataBufferInt) raster.getDataBuffer();
                    int[] bank = dataBufferByte.getData(0);
                    int w = raster.getWidth();
                    int h = raster.getHeight();
                    int size = w * raster.getHeight();
                    Rectangle r = g.getClipBounds();
                    int x2 = r.x + r.width;
                    int y2 = r.y + r.height;

                    if (s == 1) {
                        for (int i = 0; i < nPoints; i++) {
                            data.getDataPoint(i, point);
                            int x = (int) horAxis.getScreenCoord(point[0]);
                            int y = (int) verAxis.getScreenCoord(point[1]);
                            if (x >= r.x && x < x2 && y >= r.y && y < y2) {
                                bank[y * w + x] = color;
                            }
                        }
                    } else {
                        for (int i = 0; i < nPoints; i++) {
                            data.getDataPoint(i, point);
                            int x = (int) horAxis.getScreenCoord(point[0]);
                            int y = (int) verAxis.getScreenCoord(point[1]);
                            for (int xi = x - s + 1; xi <= x + s - 1; xi++) {
                                for (int yi = y - s + 1; yi <= y + s - 1; yi++) {
                                    if (xi >= r.x && xi < x2 && yi >= r.y && yi < y2) {
                                        bank[yi * w + xi] = color;
                                    }
                                }
                            }
                        }
                    }
                }
                break;

            case CIRCLE:
                for (int i = 0; i < nPoints; i++) {
                    data.getDataPoint(i, point);
                    double x = horAxis.getScreenCoord(point[0]);
                    double y = verAxis.getScreenCoord(point[1]);
                    ellipse.setFrame(x - ss, y - ss, s, s);
                    g.draw(ellipse);
                }
                break;

            case DISK:
                for (int i = 0; i < nPoints; i++) {
                    data.getDataPoint(i, point);
                    double x = horAxis.getScreenCoord(point[0]);
                    double y = verAxis.getScreenCoord(point[1]);
                    ellipse.setFrame(x - ss, y - ss, s, s);
                    g.fill(ellipse);
                }
                break;

            case VERTICAL_LINE:
                for (int i = 0; i < nPoints; i++) {
                    data.getDataPoint(i, point);
                    double x = horAxis.getScreenCoord(point[0]);
                    double y = verAxis.getScreenCoord(point[1]);
                    line.setLine(x, y - ss, x, y + ss);
                    g.draw(line);
                }
                break;

            case HORIZONTAL_LINE:
                for (int i = 0; i < nPoints; i++) {
                    data.getDataPoint(i, point);
                    double x = horAxis.getScreenCoord(point[0]);
                    double y = verAxis.getScreenCoord(point[1]);
                    line.setLine(x - ss, y, x + ss, y);
                    g.draw(line);
                }
                break;

            case CROSS:
                for (int i = 0; i < nPoints; i++) {
                    data.getDataPoint(i, point);
                    double x = horAxis.getScreenCoord(point[0]);
                    double y = verAxis.getScreenCoord(point[1]);
                    line.setLine(x - ss, y - ss, x + ss, y + ss);
                    g.draw(line);
                    line.setLine(x + ss, y - ss, x - ss, y + ss);
                    g.draw(line);
                }
                break;
        }
    }

}
