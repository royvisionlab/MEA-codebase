

package edu.ucsc.neurobiology.vision.plot;

import java.awt.*;


/**
 * This interface specifies the drawing style of a Histogram.
 * One instance of this interface is passed to the <tt>PlotPanel</tt> class when a
 * Histogram gets added to the <tt>PlotPanel</tt>.
 * The style includes colors, outline type, line sizes, etc.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 * @owner Matthew Grivich, University of California, Santa Cruz
 */
public class HistogramStyle
    implements PlotStyle {

    public enum CoordinatesType {
        CARTEZIAN, POLAR
    }


    /**
     * This is a possible return value for <tt>getOutlineType()</tt> method.
     * It specifies that the hisogram will be drawn as towers.
     */
    public enum OutlineType {
        RECTANGULAR, LINEAR
    }


    public boolean isFastAndSimple = false;

    private final String description;
    private CoordinatesType coordinatesType;
    private OutlineType outlineType;
    private Color outlineColor;
    private float outlineThickness;
    private boolean fillingTowers;
    private Color fillingColor;
    private boolean connectingTowers;
    private Color connectionLineColor;
    private float connectionLineThickness;


    public HistogramStyle(String description) {
        this(description, OutlineType.LINEAR, Color.black, 1, false, Color.yellow, false,
             Color.blue, 1);
    }


    public HistogramStyle() {
        this("", OutlineType.LINEAR, Color.black, 1, false, Color.yellow, false,
             Color.blue, 1);
    }


    public HistogramStyle(CoordinatesType coordinatesType) {
        this("", OutlineType.LINEAR, Color.black, 1, false, Color.yellow, false,
             Color.blue,
             1);
        this.coordinatesType = coordinatesType;
    }


    public HistogramStyle(Color outlineColor, int outlineThickness) {
        this("", OutlineType.LINEAR, outlineColor, outlineThickness,
             false, Color.yellow, false, Color.blue, 1);
    }

    public HistogramStyle(Color outlineColor, int outlineThickness, boolean fill) {
          this("", OutlineType.LINEAR, outlineColor, outlineThickness,
               fill, outlineColor, false, Color.blue, 1);
      }



    public HistogramStyle(
        String description, OutlineType outlineType, Color outlineColor,
        float outlineThickness, boolean fillingTowers, Color fillingColor,
        boolean connectingTowers, Color connectionLineColor,
        float connectionLineThickness) {

        this.coordinatesType = CoordinatesType.CARTEZIAN;
        this.outlineType = outlineType;
        this.outlineColor = outlineColor;
        this.outlineThickness = outlineThickness;
        this.fillingTowers = fillingTowers;
        this.fillingColor = fillingColor;
        this.connectingTowers = connectingTowers;
        this.connectionLineColor = connectionLineColor;
        this.connectionLineThickness = connectionLineThickness;
        this.description = description;
    }


    public String getDescription() {
        return description;
    }


    public CoordinatesType getCoordinatesType() {
        return coordinatesType;
    }


    /**
     * Returns the type of the outline used to draw the histogram. The return value can
     * be one of RECTANGULAR_OUTLINE or LINEAR_OUTLINE;
     */
    public OutlineType getOutlineType() {
        return outlineType;
    }


    /**
     * Returns the color of the histogram outline as a <tt>java.awt.Color</tt>.
     */
    public Color getOutlineColor() {
        return outlineColor;
    }


    /**
     * Returns the thickness of the outline in pixels. Usualy 1, 2, 3 or 4 pixels.
     */
    public float getOutlineThickness() {
        return outlineThickness;
    }


    /**
     * Controls whether or not to fill the towers with color.
     *
     * @return <b>true</b> if filling is wanted and <b>false</b> otherwise
     */
    public boolean getFillingTowers() {
        return fillingTowers;
    }


    /**
     * Returns the fill color of the histogram bins as a <tt>java.awt.Color</tt>. The
     * return value of this method is used only if <tt>isFillingTowers()</tt> returns
     * true.
     */
    public Color getFillingColor() {
        return fillingColor;
    }


    /**
     * Controls whether or not to connect the towers with a line.
     *
     * @return <b>true</b> if the connection line is needed and <b>false</b> otherwise
     */
    public boolean getConnectingTowers() {
        return connectingTowers;
    }


    /**
     * Returns the color used to draw the connecting line as a <tt>java.awt.Color</tt>.
     * The return value of this method is used only if <tt>isConnectingTowersTowers()</tt>
     * returns true.
     */
    public Color getConnectionLineColor() {
        return connectionLineColor;
    }


    /**
     * Returns the thickness of the connecting line in pixels. Usualy 1, 2, 3 or 4 pixels.
     */
    public float getConnectionLineThickness() {
        return connectionLineThickness;
    }


    public void setCoordinatesType(CoordinatesType coordinatesType) {
        this.coordinatesType = coordinatesType;
    }


    public void setOutlineType(OutlineType outlineType) {
        this.outlineType = outlineType;
    }


    public void setOutlineColor(Color outlineColor) {
        this.outlineColor = outlineColor;
    }


    public void setOutlineThickness(float outlineThickness) {
        this.outlineThickness = outlineThickness;
    }


    public void setFillingTowers(boolean fillingTowers) {
        this.fillingTowers = fillingTowers;
    }


    public void setFillingColor(Color fillingColor) {
        this.fillingColor = fillingColor;
    }


    public void setConnectingTowers(boolean connectingTowers) {
        this.connectingTowers = connectingTowers;
    }


    public void setConnectionLineColor(Color connectionLineColor) {
        this.connectionLineColor = connectionLineColor;
    }


    public void setConnectionLineThickness(float connectionLineThickness) {
        this.connectionLineThickness = connectionLineThickness;
    }
}

