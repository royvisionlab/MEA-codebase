package edu.ucsc.neurobiology.vision.plot;

import java.lang.reflect.*;

import java.awt.*;

import edu.ucsc.neurobiology.vision.util.*;


/**
 * This interface specifies the drawing style of a ScatterPlot.
 * One instance of this interface is passed to the <tt>PlotPanel</tt> class when a
 * ScatterPlot gets added to the <tt>PlotPanel</tt>.
 * The style include colors, symbol types and sizes, etc.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ScatterPlotStyle
    implements PlotStyle {

    public enum LineStyle {
        SOLID, DOTTED, DASHED, DOTTED_DASHED
    }


    private String dashPattern;
    private SymbolType symbolType;
    private SymbolType errorSymbolType = SymbolType.VERTICAL_LINE;
    private int symbolSize;
    private Color symbolColor;
    private boolean connectingPoints;
    private Color connectionLineColor;
    private float connectionLineThickness;
    private int connectionPeriod = Integer.MAX_VALUE;
    private LineStyle lineStyle = LineStyle.SOLID;
    private final String description;
    private static int styleCount = 1;


    public void setSymbolStyle(String style) {
        if (style.equals("")) { // no symbol info
            this.setSymbolType(SymbolType.NONE);
            this.setSymbolColor(Color.black);
        } else {
            String[] s = StringUtil.decomposeString(style, " ");
            if (s.length != 3) {
                throw new Error("3 params should be given for symbols");
            }

            this.setSymbolType(SymbolType.valueOf(s[0]));
            this.setSymbolSize(Integer.parseInt(s[1]));
            this.setSymbolColor(getColor(s[2]));
        }
    }


    public void setLineStyle(String style) {
        if (style.equals("")) {
            this.setConnectingPoints(false);
            this.setConnectionLineColor(Color.black);
        } else {
            String[] s = StringUtil.decomposeString(style, " ");
            if (s.length != 3) {
                throw new Error("3 params should be given for symbols");
            }

            this.setConnectingPoints(true);
            this.setConnectionLineStyle(LineStyle.valueOf(s[0]));
            this.setConnectionLineThickness( (float) Double.parseDouble(s[1]));
            this.setConnectionLineColor(getColor(s[2]));
        }
    }


    public ScatterPlotStyle(String styleString) {
        String[] subStyles = StringUtil.decomposeString(styleString, ",");

        if (subStyles.length == 0) {
//            System.err.println("0p");
            setSymbolStyle("");
            setLineStyle("");
        }

        if (subStyles.length == 1) {
//            System.err.println("1p");
            setSymbolStyle(subStyles[0]);
        }

        if (subStyles.length == 2) {
//            System.err.println("2p");
            setSymbolStyle(subStyles[0]);
            setLineStyle(subStyles[1]);
        }

        this.description = "SCatterPlot " + styleCount++;
    }


    public ScatterPlotStyle() {
        this("SCatterPlot " + styleCount++, SymbolType.FILLED_SQUARE, 1, Color.black, false,
             Color.black, 1);
    }


    /**
     *  For unconnected plots
     */
    public ScatterPlotStyle(String name, SymbolType symbolType, int symbolSize,
                            Color symbolColor) {
        this(name, symbolType, symbolSize, symbolColor, false, Color.black, 1);
    }


    /**
     *  For unconnected plots
     */
    public ScatterPlotStyle(SymbolType symbolType, int symbolSize, Color symbolColor) {
        this("SCatterPlot " + styleCount++, symbolType, symbolSize, symbolColor, false,
             Color.black, 1);
    }


    /**
     * For connected plots with no symbols
     */
    public ScatterPlotStyle(
        boolean connectingPoints, Color connectionLineColor,
        float connectionLineThickness) {

        this("SCatterPlot " + styleCount++, SymbolType.NONE, 1, Color.black,
             connectingPoints, connectionLineColor, connectionLineThickness);
    }


    public ScatterPlotStyle(
        SymbolType symbolType, int symbolSize, Color symbolColor,
        boolean connectingPoints, Color connectionLineColor,
        float connectionLineThickness) {

        this("SCatterPlot " + styleCount++, symbolType, symbolSize, symbolColor,
             connectingPoints,
             connectionLineColor, connectionLineThickness);
    }


    public ScatterPlotStyle(
        String description, SymbolType symbolType, int symbolSize, Color symbolColor,
        boolean connectingPoints, Color connectionLineColor,
        float connectionLineThickness) {

        this.symbolType = symbolType;
        this.symbolSize = symbolSize;
        this.symbolColor = symbolColor;
        this.connectingPoints = connectingPoints;
        this.connectionLineColor = connectionLineColor;
        this.connectionLineThickness = connectionLineThickness;
        this.description = description;
    }


    public String getDashPattern() {
        return dashPattern;
    }


    public void setDashPattern(String dashPattern) {
        this.dashPattern = dashPattern;
    }


    public String getDescription() {
        return description;
    }


    /**
     * Returns the symbol type used to draw the scatter plot. One of SYMBOLTYPE_NONE,
     * SQUARE or CIRCLE.
     */
    public SymbolType getSymbolType() {
        return symbolType;
    }


    public SymbolType getErrorSymbolType() {
        return errorSymbolType;
    }


    /**
     * Returns the size in pixels of the symbols (usually from 0 to 10).
     */
    public int getSymbolSize() {
        return symbolSize;
    }


    /**
     * Returns the <tt>Color</tt> of symbols (could be any).
     */
    public Color getSymbolColor() {
        return symbolColor;
    }


    /**
     * Returns <tt>true</tt> if the points of the scatter plot should be connected with
     * lines, <tt>false</tt> otherwise.
     */
    public boolean getConnectingPoints() {
        return connectingPoints;
    }


    /**
     * Returns the <tt>Color</tt> of connection lines (could be any).
     */
    public Color getConnectionLineColor() {
        return connectionLineColor;
    }


    /**
     * Returns the thickness in pixels of connection lines (usualy 1-4 pixels).
     */
    public float getConnectionLineThickness() {
        return connectionLineThickness;
    }


    public void setSymbolType(SymbolType symbolType) {
        this.symbolType = symbolType;
    }


    public void setErrorSymbolType(SymbolType errorSymbolType) {
        this.errorSymbolType = errorSymbolType;
    }


    public void setSymbolSize(int symbolSize) {
        this.symbolSize = symbolSize;
    }


    public void setSymbolColor(Color symbolColor) {
        this.symbolColor = symbolColor;
    }


    public void setConnectingPoints(boolean connectingPoints) {
        this.connectingPoints = connectingPoints;
    }


    public void setConnectionLineColor(Color connectionLineColor) {
        this.connectionLineColor = connectionLineColor;
    }


    public void setConnectionLineThickness(float connectionLineThickness) {
        this.connectionLineThickness = connectionLineThickness;
    }


    public int getConnectionPeriod() {
        return connectionPeriod;
    }


    public void setConnectionPeriod(int period) {
        this.connectionPeriod = period;
    }


    public float[] getConnectionLineDash() {
        switch (lineStyle) {
            case DASHED:
                return new float[] {
                    6, 6};

            case DOTTED:
                return new float[] {
                    1, 6};

            case DOTTED_DASHED:
                return new float[] {
                    6, 5, 1, 5};

            default:
                return new float[] {
                    1};
        }
    }


    public void setConnectionLineStyle(LineStyle lineStyle) {
        this.lineStyle = lineStyle;
    }


    public LineStyle getConnectionLineStyle() {
        return lineStyle;
    }


    public static Color getColor(String c) {
        Field[] f = Color.class.getFields();
        for (int i = 0; i < f.length; i++) {
            if (f[i].getDeclaringClass().equals(Color.class)) {
                if (f[i].getName().equals(c)) {
                    try {
                        return (Color) f[i].get(Color.white);
                    } catch (IllegalAccessException ex) {
                        ex.printStackTrace();
                    }
                }
            }
        }

        return null;
    }
}
