package edu.ucsc.neurobiology.vision.plot;

import java.awt.*;


/**
 * This interface specifies the drawing style of a Function.
 * One instance of this interface is passed to the <tt>PlotPanel</tt> class when a
 * Function or a ParametricFunction gets added to the <tt>PlotPanel</tt>.
 * The style includes the colors and the line thickness.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class FunctionStyle implements PlotStyle {

    private Color lineColor;
    private float lineThickness;
    private double min = Double.NEGATIVE_INFINITY, max = Double.POSITIVE_INFINITY;
    private String dashPattern;
    private final String description;
    private final Font labelFont;


    public FunctionStyle(String description) {
        this(description, Color.black, 1);
    }

    public FunctionStyle(String description, Color lineColor, float lineThickness) {
        this(description, lineColor, lineThickness, "Arial 10");
    }
    
    public FunctionStyle(String description, Color lineColor, float lineThickness, String fontString) {
        this(description, lineColor, lineThickness, Font.decode(fontString));
    }

    public FunctionStyle(String description, Color lineColor, float lineThickness, Font labelFont) {
        this.lineColor = lineColor;
        this.lineThickness = lineThickness;
        this.description = description;
        this.labelFont = labelFont;
    }


    public String getDescription() {
        return description;
    }


    /**
     * Returns the line color as a <tt>java.awt.Color</tt> class (could be any).
     */
    public Color getLineColor() {
        return lineColor;
    }


    /**
     * Returns the thickness in pixels of the lines used to draw the function
     * (usualy 1-4 pixels).
     */
    public float getLineThickness() {
        return lineThickness;
    }


    public void setLineColor(Color lineColor) {
        this.lineColor = lineColor;
    }


    public void setLineThickness(float lineThickness) {
        this.lineThickness = lineThickness;
    }


    public double getMinimum() {
        return min;
    }


    public double getMaximum() {
        return max;
    }


    public void setMinimum(double min) {
        this.min = min;
    }


    public void setMaximum(double max) {
        this.max = max;
    }


    public String getDashPattern() {
        return dashPattern;
    }


    public void setDashPattern(String dashPattern) {
        this.dashPattern = dashPattern;
    }
    
    public Font getLabelFont() { return labelFont; }
}
