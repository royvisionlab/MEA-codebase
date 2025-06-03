package edu.ucsc.neurobiology.vision.plot;

import java.awt.BasicStroke;
import java.awt.Component;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Insets;
import java.awt.RenderingHints;

import javax.swing.border.Border;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class AxesBorder implements Border {

    private Axis xAxis;
    private Axis yAxis;

    int topPadding = 0;
    int bottomPadding = 0;
    int rightPadding = 0;
    int leftPadding = 0;


    public AxesBorder() {
        this.xAxis = new HorizontalAxis();
        this.yAxis = new VerticalAxis();
    }


    /**
     * Needs a repaint() to take effect !
     */
    public void setXAxisVisible(boolean xAxisVisible) {
        xAxis.showTicks = xAxisVisible;
        xAxis.showLabel = xAxisVisible;
    }


    public void setYAxisVisible(boolean yAxisVisible) {
        yAxis.showTicks = yAxisVisible;
        yAxis.showLabel = yAxisVisible;
    }


    public Axis getHorizontalAxis() {
        return xAxis;
    }


    public Axis getVerticalAxis() {
        return yAxis;
    }


    public boolean isBorderOpaque() {
        return false;
    }


    public void setPadding(int topPadding, int bottomPadding, int leftPadding,
                           int rightPadding) {
        this.topPadding = topPadding;
        this.bottomPadding = bottomPadding;
        this.rightPadding = rightPadding;
        this.leftPadding = leftPadding;
    }


    public void setTopPadding(int topPadding) {
        this.topPadding = topPadding;
    }


    public int getTopPadding() {
        return topPadding;
    }


    public void setBottomPadding(int bottomPadding) {
        this.bottomPadding = bottomPadding;
    }


    public int getBottomPadding() {
        return bottomPadding;
    }


    public void setRightPadding(int rightPadding) {
        this.rightPadding = rightPadding;
    }


    public int getRightPadding() {
        return rightPadding;
    }


    public void setLeftPadding(int leftPadding) {
        this.leftPadding = leftPadding;
    }


    public int getLeftPadding() {
        return leftPadding;
    }


    public Insets getBorderInsets(Component c) {
        int topInset = topPadding;
        int rightInset = rightPadding;
        int leftInset = yAxis.getThickness() + leftPadding;
        int bottomInset = xAxis.getThickness() + bottomPadding;

        Insets insets = new Insets(topInset, leftInset, bottomInset, rightInset);
        xAxis.setScreenRange(insets.left, c.getWidth() - insets.right);
        yAxis.setScreenRange(insets.top, c.getHeight() - insets.bottom);

        return insets;
    }


    public void paintBorder(Component c, Graphics g, int x, int y, int width, int height) {
        Graphics2D g2 = (Graphics2D) g;
        g2.setRenderingHint(
            RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        g2.setStroke(new BasicStroke(0.5f));

        Insets insets = getBorderInsets(c);
        xAxis.setScreenRange(insets.left, c.getWidth() - insets.right);
        yAxis.setScreenRange(insets.top, c.getHeight() - insets.bottom);

        xAxis.paintAxis(g2, insets.left, height - insets.bottom, width, height);
        yAxis.paintAxis(g2, insets.left, insets.top, width, height);

        g2.setRenderingHint(
            RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_OFF);
    }


}
