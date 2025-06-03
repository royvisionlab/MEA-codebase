package edu.ucsc.neurobiology.vision.plot;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ColorPlotStyle
    implements PlotStyle {

    private ColorStyle colorStyle;

    public enum ColorStyle {
        RGB, R, G, B, R_AND_G, R_MINUS_G
    }


    private final String description;


    public ColorPlotStyle(String description) {
        colorStyle = ColorStyle.RGB;
        this.description = description;
    }


    public String getDescription() {
        return description;
    }


//    public boolean getGridPainted() {
//        return false;
//    }
//

//    public Color getGridColor() {
//        return Color.white;
//    }


    public int getLineThickness() {
        return 1;
    }


    public ColorStyle getColorStyle() {
        return colorStyle;
    }


    public void setColorStyle(ColorStyle colorStyle) {
        this.colorStyle = colorStyle;
    }
}
