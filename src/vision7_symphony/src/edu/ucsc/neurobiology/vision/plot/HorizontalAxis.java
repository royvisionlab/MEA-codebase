package edu.ucsc.neurobiology.vision.plot;

import static java.lang.Math.*;

import java.awt.*;
import java.awt.geom.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class HorizontalAxis
    extends Axis {

    public static final int CENTER_TICK_ALLIGNMENT = -1;
    public static final int RIGHT_TICK_ALLIGNMENT = -2;


    public HorizontalAxis() {
    }


    public double getK() {
        return (convert(x2) - convert(x1)) / (s2 - s1);
    }


    public double getScreenCoord(double x) {
        return s1 + (1 / getK()) * (convert(x) - convert(x1));
    }


    public double getPlotCoord(double s) {
        if (type == AxisType.LINEAR) {
            return x1 + getK() * (s - s1);
        } else {
            return x1 * Math.pow(10, getK() * (s - s1));
        }
    }


    synchronized public int updateTicks() {
        double[] x = AxisUtil.performAutoScale(x1, x2);
        maxTickMarkSize = max(metrics.stringWidth("" + x[0]),
                              metrics.stringWidth("" + x[1]));

        if (x1 != 0 || x2 != 0) { // create the ticks
            if (type == AxisType.LINEAR) {
                if(manualTicks != null) {
                    createManualTicks();
                } else if (!Double.isNaN(fixedTick)) {
                    createLinearFixedTicks();
                } else {
                    int factor = 1;
                    again:while (true) {
                        createLinearTicks(factor);

                        int lastTick = Integer.MIN_VALUE;
                        int lastW = 0;
                        for (Double _tick : ticksAndLabels.keySet()) {
                            String label = ticksAndLabels.get(_tick);
                            if (label == null) {
                                continue;
                            }

                            int tick = (int) Math.round(getScreenCoord(_tick));
                            int w = metrics.stringWidth(label);
                            if (abs( (tick - w / 2) - (lastTick + lastW / 2)) <
                                minTickSeparation) {
                                factor *= 2;
                                continue again;
                            }
                            lastTick = tick;
                            lastW = w;
                        }

                        break;
                    }
                }
            } else {
                createLogTicks(0);
            }
        }

        int thickness = 0;
        if (showTicks) {
            thickness += tickSize;
        }

        if (showTicks) {
            thickness += tickToMarkSpacing + metrics.getAscent() - metrics.getDescent();
        }

        if (showLabel) {
            if (labelComponent.getText() != null &&
                labelComponent.getText().trim().length() > 0) {

                thickness += markToLabelSpacing + labelComponent.getHeight();
            }
        }

        return thickness;
    }


    synchronized public void paintAxis(Graphics2D g, int x0, int y0, int width,
                                       int height) {
        g.setFont(currentFont);
        g.setColor(Color.black);

        float yLocation = y0;
        if (showTicks) {
            yLocation += metrics.getAscent() + tickToMarkSpacing;
            // draw the tick lines and the tick labels
            for (Double tick : ticksAndLabels.keySet()) {
                float xTick = (float) getScreenCoord(tick);

                String label = ticksAndLabels.get(tick);
                if (label != null) {
                    int sign = (tick.doubleValue() < 0) ? metrics.stringWidth("-") : 0;
                    int w = metrics.stringWidth(label);
                    float xTickScreen;
//                    if (shiftRightLabel && xTickScreen + w > width) {
//                        xTickScreen = width - w - 1;
//                    }
                    switch (tickMarkAllignment) {
                        case HorizontalAxis.RIGHT_TICK_ALLIGNMENT:
                            xTickScreen = xTick - w - 1;
                            break;

                        default:
                            xTickScreen = xTick - (w - sign) / 2.f - sign;
                            break;
                    }

                    g.drawString(label, xTickScreen, yLocation);
                    g.draw(new Line2D.Float(xTick, y0, xTick, y0 + tickSize));
                } else {
                    g.draw(new Line2D.Float(xTick, y0, xTick, y0 + tickSize / 2));
                }
            } // for
        }

        // draw the axis label
//        String _label = (secondaryLabel == null) ? label : label + " " + secondaryLabel;

        // draw the label
        if (showLabel) {
            if (fixedLabelPosition < 0) {
                yLocation += tickSize - metrics.getDescent() + markToLabelSpacing;
            } else {
                yLocation = fixedLabelPosition;
            }
            double x = (s2 + s1) / 2.f - labelComponent.getWidth() / 2.f;
            g.translate(x, yLocation);
            g.setColor(Color.black);
            labelComponent.paint(g);
            g.translate( -x, -yLocation);
        }

        // draw the grid lines
        if (showGridLines) {
            g.setColor(Color.lightGray);
            g.setStroke(new BasicStroke(
                1, BasicStroke.CAP_SQUARE,
                BasicStroke.JOIN_MITER, 2, new float[] {5, 5, 1}, 0));
            for (Double tick : ticksAndLabels.keySet()) {
                int xTick = (int) Math.round(getScreenCoord(tick));
                if (Math.abs(xTick - s1) > gridPadding &&
                    Math.abs(xTick - s2) > gridPadding) {
                    g.drawLine(xTick, 0, xTick, y0);
                }
            }
        }
    }

}
