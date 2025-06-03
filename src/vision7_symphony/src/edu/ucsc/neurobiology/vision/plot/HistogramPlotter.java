package edu.ucsc.neurobiology.vision.plot;

import java.util.*;

import java.awt.*;
import java.awt.geom.*;
import java.awt.image.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class HistogramPlotter
    implements DataPlotter {

    Rectangle2D.Double rect = new Rectangle2D.Double();
    Line2D.Double line = new Line2D.Double();


    public synchronized boolean accept(Object data, Object style) {
        if ( (data instanceof HistogramData) && (style instanceof HistogramStyle)) {
            return true;
        } else {
            return false;
        }
    }


    public synchronized void draw(Axis horAxis, Axis verAxis, Graphics g, Raster raster,
                                  Object _data, Object _style) {

        HistogramData data = (HistogramData) _data;
        HistogramStyle style = (HistogramStyle) _style;

        switch (style.getCoordinatesType()) {
            case CARTEZIAN:
                if (style.isFastAndSimple && raster != null) {
                    drawCartezianDirect(
                        horAxis, verAxis, (Graphics2D) g, raster, data, style);
                } else {
                    drawCartezian(horAxis, verAxis, (Graphics2D) g, data, style);
                }
                break;

            case POLAR:
                drawPolar(horAxis, verAxis, (Graphics2D) g, data, style);
                break;

            default:
                throw new IllegalStateException("Wrong coordinates");
        }
    }


    private synchronized void drawPolar(Axis horAxis, Axis verAxis, Graphics2D g,
                                        HistogramData data, HistogramStyle style) {

        g.setStroke(new BasicStroke(style.getOutlineThickness()));

        int[] x = horAxis.getScreenRange();
        int[] y = verAxis.getScreenRange();
        final int x0 = (x[0] + x[1]) / 2;
        final int y0 = (y[0] + y[1]) / 2;
        final int rMax = Math.min(x[1] - x0, y[1] - y0);
        final int rMin = rMax / 10;
        final int nBins = data.getBinCount();
        double vMax = Double.NEGATIVE_INFINITY;
//        double vMin = Double.POSITIVE_INFINITY;
        for (int i = 0; i < nBins; i++) {
            if (data.getBin(i) > vMax) {
                vMax = data.getBin(i);
//            if (data.getBin(i) < vMin) vMin = data.getBin(i);
            }
        }

        final double dPhi = data.getBinInterval();
        for (int i = 0; i < nBins; i++) {
            double r = rMin + (rMax - rMin) * data.getBin(i) / vMax;

            double phi1 = (data.getMin() + (i - 0.5 + 0.1) * dPhi) * Math.PI / 180;
            double phi2 = (data.getMin() + (i + 0.5 - 0.1) * dPhi) * Math.PI / 180;

            Polygon2D p = new Polygon2D();
            p.addPoint(x0 + rMin * Math.cos(phi1), y0 - rMin * Math.sin(phi1));
            p.addPoint(x0 + rMin * Math.cos(phi2), y0 - rMin * Math.sin(phi2));
            p.addPoint(x0 + r * Math.cos(phi2), y0 - r * Math.sin(phi2));
            p.addPoint(x0 + r * Math.cos(phi1), y0 - r * Math.sin(phi1));

            g.setPaint(Color.yellow);
            g.fill(p);
            g.setPaint(Color.black);
            g.draw(p);
        }

        g.drawOval(x0 - rMin, y0 - rMin, 2 * rMin, 2 * rMin);
    }


    private synchronized void drawCartezian(Axis horAxis, Axis verAxis, Graphics2D g,
                                            HistogramData data, HistogramStyle style) {

        int nBins = data.getBinCount();
        double y0 = verAxis.getScreenCoord(0);
        final double binInterval = data.getBinInterval();
        final double min = data.getMin();

        double xOld = horAxis.getScreenCoord(min);
        double x = horAxis.getScreenCoord(min + binInterval);
        double y, yOld = -1;

        for (int i = 0; i < nBins; i++) {
            y = verAxis.getScreenCoord(data.getBin(i));

            // draw the outline
            g.setStroke(new BasicStroke(
                style.getOutlineThickness(), BasicStroke.CAP_SQUARE,
                BasicStroke.JOIN_MITER));
            switch (style.getOutlineType()) {
                case RECTANGULAR:
                    if (style.getFillingTowers()) {
                        g.setColor(style.getFillingColor());
                        rect.setRect(xOld, y, x - xOld, y0 - y);
                        g.fill(rect);
                    }

                    g.setColor(style.getOutlineColor());
                    rect.setRect(xOld, y, x - xOld, y0 - y);
                    g.draw(rect);
                    break;

                case LINEAR:
                    if (style.getFillingTowers()) {
                        g.setColor(style.getFillingColor());
                        rect.setRect(xOld, y, x - xOld, y0 - y);
                        g.fill(rect);
                    }

                    g.setColor(style.getOutlineColor());
                    line.setLine(xOld, y, x, y);
                    g.draw(line);
                    if (yOld != -1) {
                        line.setLine(xOld, yOld, xOld, y);
                        g.draw(line);
                    }
                    break;
            }

            // draw the error
            if (data instanceof HistogramErrors) {
                double error = ( (HistogramErrors) data).getBinError(i);
                double xMiddle = (xOld + x) / 2;
                double yUp = verAxis.getScreenCoord(data.getBin(i) + error / 2);
                double yDown = verAxis.getScreenCoord(data.getBin(i) - error / 2);
                line.setLine(xMiddle, yUp, xMiddle, yDown);
                g.draw(line);
                line.setLine(xMiddle - 2, yUp, xMiddle + 2, yUp);
                g.draw(line);
                line.setLine(xMiddle - 2, yDown, xMiddle + 2, yDown);
                g.draw(line);
            }

            if (style.getConnectingTowers() && (yOld != -1)) {
                g.setStroke(new BasicStroke(
                    style.getConnectionLineThickness(), BasicStroke.CAP_SQUARE,
                    BasicStroke.JOIN_MITER));
                g.setColor(style.getConnectionLineColor());
                double w = (x - xOld) / 2;
                line.setLine(xOld - w, yOld, x - w, y);
                g.draw(line);
            }

            xOld = x;
            x = horAxis.getScreenCoord(min + binInterval * (i + 2));
            yOld = y;
        }
    }


    public void drawHLine(int[] bank, int w, int y, int x1, int x2, int color,
                          int X1, int X2, int Y1, int Y2) {

        if (y < Y1 || y > Y2 || (x1 < X1 && x2 < X1) || (x1 > X2 && x2 > X2)) {
            // the line is not visible
            return;
        }

        if (x1 < X1) {
            x1 = X1;
        }

        if (x2 > X2) {
            x2 = X2;
        }

        Arrays.fill(bank, y * w + x1, y * w + x2, color);
    }


    public void drawVLine(int[] bank, int w, int x, int y1, int y2, int color,
                          int X1, int X2, int Y1, int Y2) {

        if (y2 < y1) {
            int temp = y1;
            y1 = y2;
            y2 = temp;
        }

        if (x < X1 || x > X2 || (y1 < Y1 && y2 < Y1) || (y1 > Y2 && y2 > Y2)) {
            // the line is not visible
            return;
        }

        if (y1 < Y1) {
            y1 = Y1;
        }

        if (y2 > Y2) {
            y2 = Y2;
        }

        if (y2 < y1) {
            throw new Error();
        }

        int index = y1 * w + x;
        for (int y = y1; y <= y2; y++) {
            bank[index] = color;
            index += w;
        }
    }


    private synchronized void drawCartezianDirect(Axis horAxis, Axis verAxis,
                                                  Graphics2D g, Raster raster,
                                                  HistogramData data,
                                                  HistogramStyle style) {

        // prepare the image buffer
        DataBufferInt dataBufferByte = (DataBufferInt) raster.getDataBuffer();
        int[] bank = dataBufferByte.getData(0);
        int w = raster.getWidth();
        int h = raster.getHeight();
        int size = w * raster.getHeight();
        Rectangle r = g.getClipBounds();
        int X1 = r.x;
        int Y1 = r.y;
        int X2 = r.x + r.width;
        int Y2 = r.y + r.height;

//        System.err.println(X1 + ", " + X2);
//        System.err.println(Y1 + ", " + Y2);

        // prepare data variables
        final int nBins = data.getBinCount();
        final double y0 = verAxis.getScreenCoord(0);
        final double binInterval = data.getBinInterval();
        final double min = data.getMin();
        double xOld = horAxis.getScreenCoord(min);
        double x = horAxis.getScreenCoord(min + binInterval);
        double y, yOld = -1;

        // get style properties
        int outlineColor = style.getOutlineColor().getRGB();

        for (int i = 0; i < nBins; i++) {
            y = verAxis.getScreenCoord(data.getBin(i));

            // draw the outline
            g.setStroke(new BasicStroke(
                style.getOutlineThickness(), BasicStroke.CAP_SQUARE,
                BasicStroke.JOIN_MITER));
            switch (style.getOutlineType()) {
                case RECTANGULAR:

//                    if (style.getFillingTowers()) {
//                        g.setColor(style.getFillingColor());
//                        rect.setRect(xOld, y, x - xOld, y0 - y);
//                        g.fill(rect);
//                    }

                    g.setColor(style.getOutlineColor());
                    rect.setRect(xOld, y, x - xOld, y0 - y);
                    g.draw(rect);
                    break;

                case LINEAR:

//                    if (style.getFillingTowers()) {
//                        g.setColor(style.getFillingColor());
//                        rect.setRect(xOld, y, x - xOld, y0 - y);
//                        g.fill(rect);
//                    }

//                    g.setColor(style.getOutlineColor());

//                line.setLine(xOld, y, x, y);
//                g.draw(line);
//                if (yOld != -1) {
//                    line.setLine(xOld, yOld, xOld, y);
//                    g.draw(line);
//                }


//                    g.drawLine( (int) xOld, (int) y, (int) x, (int) y);
                    drawHLine(bank, w, (int) y, (int) xOld, (int) x, outlineColor,
                              X1, X2, Y1, Y2);
                    if (yOld != -1) {
//                        g.drawLine( (int) xOld, (int) yOld, (int) xOld, (int) y);
                        drawVLine(bank, w, (int) xOld, (int) yOld, (int) y, outlineColor,
                                  X1, X2, Y1, Y2);
                    }

                    break;
            }

            // draw the error
            if (data instanceof HistogramErrors) {
                double error = ( (HistogramErrors) data).getBinError(i);
                double xMiddle = (xOld + x) / 2;
                double yUp = verAxis.getScreenCoord(data.getBin(i) + error / 2);
                double yDown = verAxis.getScreenCoord(data.getBin(i) - error / 2);
                line.setLine(xMiddle, yUp, xMiddle, yDown);
                g.draw(line);
                line.setLine(xMiddle - 2, yUp, xMiddle + 2, yUp);
                g.draw(line);
                line.setLine(xMiddle - 2, yDown, xMiddle + 2, yDown);
                g.draw(line);
            }

            if (style.getConnectingTowers() && (yOld != -1)) {
                g.setStroke(new BasicStroke(
                    style.getConnectionLineThickness(), BasicStroke.CAP_SQUARE,
                    BasicStroke.JOIN_MITER));
                g.setColor(style.getConnectionLineColor());
                double ww = (x - xOld) / 2;
                line.setLine(xOld - w, yOld, x - ww, y);
                g.draw(line);
            }

            xOld = x;
            x = horAxis.getScreenCoord(min + binInterval * (i + 2));
            yOld = y;
        }
    }

}
