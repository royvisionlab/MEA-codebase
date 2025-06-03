package edu.ucsc.neurobiology.vision.plot;

import static java.lang.Math.*;

import java.awt.*;
import java.awt.geom.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class VerticalAxis
    extends Axis {


    public VerticalAxis() {
    }


    public double getK() {
        return (convert(x1) - convert(x2)) / (s2 - s1);
    }


    public double getScreenCoord(double x) {
        return s1 + (1 / getK()) * (convert(x) - convert(x2));
    }


    public double getPlotCoord(double s) {
        if (type == AxisType.LINEAR) {
            return x2 + getK() * (s - s1);
        } else {
            return x2 * Math.pow(10, getK() * (s - s1));
        }
    }


    synchronized public int updateTicks() {
        if (x1 != 0 || x2 != 0) { // create the ticks
            
            if (type == AxisType.LINEAR) {
                if (manualTicks != null) {
                    createManualTicks();
                } else if (!Double.isNaN(fixedTick)) {
                    createLinearFixedTicks();
                } else {
                    int factor = 1;
                    again:while (true) {
                        createLinearTicks(factor);

                        int h = metrics.getAscent();
                        int lastTick = Integer.MIN_VALUE;
                        for (Double _tick : ticksAndLabels.keySet()) {
                            String label = ticksAndLabels.get(_tick);
                            if (label == null) {
                                continue;
                            }

                            int tick = (int) Math.round(getScreenCoord(_tick));
                            if (abs( (lastTick - h / 2) - (tick + h / 2)) <
                                minTickSeparation) {
                                factor *= 2;
                                continue again;
                            }
                            lastTick = tick;
                        }

                        break;
                    }
                }
            } else if (type == AxisType.LOG10) {
                createLogTicks(0);
            }
        }

        maxTickMarkSize = 0;
        for (String label : ticksAndLabels.values()) {
            if (label != null) {
                int w = metrics.stringWidth(label);
                if (w > maxTickMarkSize) {
                    maxTickMarkSize = w;
                }
            }
        }

        int thickness = 0;
        if (showTicks) {
            thickness += maxTickMarkSize + tickSize + tickToMarkSpacing;
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
        float xLocation = x0;

        // draw the tick lines and the tick labels
        if (showTicks) {
            xLocation -= tickSize + tickToMarkSpacing;
            maxTickMarkSize = -1;
            for (Double tick : ticksAndLabels.keySet()) {
                float yTick = (float) getScreenCoord(tick);

                // draw label
                String label = ticksAndLabels.get(tick);
                if (label != null) {
                    int w = metrics.stringWidth(label);
                    if (w > maxTickMarkSize) {
                        maxTickMarkSize = w;
                    }
                    g.drawString(label, xLocation - w,
                                 yTick + metrics.getAscent() / 2.8f);
                    g.draw(new Line2D.Float(x0 - tickSize, yTick, x0, yTick));
                } else {
                    g.draw(new Line2D.Float(x0 - tickSize / 2, yTick, x0, yTick));
                }
            }

            xLocation -= maxTickMarkSize;
        }

        // draw the label
        if (showLabel) {
            if (fixedLabelPosition < 0) {
                xLocation -= markToLabelSpacing + labelComponent.getHeight();
            } else {
                xLocation = fixedLabelPosition;
            }

            double y;
            if (labelComponent.getWidth() <= Math.abs(s2 - s1)) {
                y = (s2 + s1) / 2. + labelComponent.getWidth() / 2.;
            } else {
                y = labelComponent.getWidth();
//                y = height / 2. + labelComponent.getWidth() / 2.;
            }

            g.translate(xLocation, y);
            g.rotate( -Math.PI / 2);
            labelComponent.paint(g);
            g.rotate( +Math.PI / 2);
            g.translate( -xLocation, y);
        }

//        g.setFont(currentLabelFont);
//        int x = labelFontMetrics.getAscent();
//        int y = (s2 + s1) / 2 + labelFontMetrics.stringWidth(label) / 2;
//        g.transform(AffineTransform.getRotateInstance( -Math.PI / 2, x, y));
//        g.drawString(label, x, y);
//        g.transform(AffineTransform.getRotateInstance(Math.PI / 2, x, y));
        /*
                // draw the grid lines
                if (showGridLines) {
                    g.setColor(Color.lightGray);
                    g.setStroke(new BasicStroke(
                        1, BasicStroke.CAP_SQUARE,
                        BasicStroke.JOIN_MITER, 2, new float[] {5, 5, 1}, 0));
                    for (Double tick : tickaAndLabels.keySet()) {
                        int yTick = (int) Math.round(getScreenCoord(tick));
                        if (Math.abs(yTick - s1) > gridPadding &&
                            Math.abs(yTick - s2) > gridPadding) {
                            g.drawLine(x0, yTick, thickness + width - 1, yTick);
                        }
                    }
                }
         */
    }

}
