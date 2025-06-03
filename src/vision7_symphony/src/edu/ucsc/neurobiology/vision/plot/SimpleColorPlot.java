package edu.ucsc.neurobiology.vision.plot;

import java.io.*;

import edu.ucsc.neurobiology.vision.stimulus.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class SimpleColorPlot
implements ColorPlotData, ChangeableData {

    public static final int ADD_GRAY = 0;
    public static final int NORMALIZE = 1;
    public static final int UNCHANGED = 2;
    public static final int NORMALIZE_ALTERNATE = 3;

    private final int style;
    private float scale;
    private ImageFrame frame;
    private DrawingControl drawingControl;
    private double minColor, maxColor;
    public float overscale = 1;
    private double minX = Double.NaN, maxX = Double.NaN, minY = Double.NaN, maxY = Double.NaN;


    public SimpleColorPlot(int style) {
        this.style = style;
    }


    /**
     * Calculates the factor needed for scaling the STA's.
     */
    private double getScaleFor(Movie movie) throws IOException {
        double maxColor = Double.NEGATIVE_INFINITY, minColor = Double.POSITIVE_INFINITY;

        for (int f = 0; f < movie.size(); f++) {
            ImageFrame frame = movie.getFrame(f);
            for (int i = 0; i < movie.getWidth(); i++) {
                for (int j = 0; j < movie.getHeight(); j++) {
                    for (int c = 0; c < 3; c++) {
                        double color = frame.getPixel(i, j, c);
//						System.out.println(color);
                        if (color > maxColor) {
                            maxColor = color;
                        }
                        if (color < minColor) {
                            minColor = color;
                        }
                    }
                }
            }

        }
        this.minColor = minColor;
        this.maxColor = maxColor;
        // divided by 2 because the max value of the factor should be 0.5
        return Math.min(1 / Math.abs(2 * minColor), 1 / Math.abs(2 * maxColor));
    }


    public void setFrame(Movie movie, int f) throws IOException {
        if (style == NORMALIZE || style == NORMALIZE_ALTERNATE) {
            scale = (float) getScaleFor(movie);
        }
        this.frame = movie.getFrame(f);

        if (drawingControl != null) {
            drawingControl.updateNeeded(this);
        }
    }


    public void setDrawingControl(DrawingControl drawingControl) {
        this.drawingControl = drawingControl;
    }


    public int getRowsCount() {
        return frame.getHeight();
    }


    public int getColumnsCount() {
        return frame.getWidth();
    }


    /**
     * Must behave reasonably when read from multiple threads concurrently.  
     * FIXME: Should get further thread-safing.
     */
    public float[] getCell(int x, int y, float[] pixel) {
        frame.getPixel(x, y, pixel);
        
        switch (style) {
        case NORMALIZE:
            for (int i = 0; i < pixel.length; i++) {
                float v = (0.5f + pixel[i] * scale * overscale);
                if (v < 0) {
                    v = 0;
                } else if (v > 1) {
                    v = 1;
                }
                pixel[i] = v;
            }
            return pixel;

        case NORMALIZE_ALTERNATE:
            for (int i = 0; i < pixel.length; i++) {
                float v = (pixel[i] - (float)this.minColor) /
                ( (float) this.maxColor - (float)this.minColor);
                if (v < 0) {
                    v = 0;
                } else if (v > 1) {
                    v = 1;
                }
                pixel[i] = v;
            }
            return pixel;

        case ADD_GRAY:
            for (int i = 0; i < pixel.length; i++) {
                pixel[i] += 0.5f;
            }
            return pixel;

        case UNCHANGED:
            return pixel;

        default:
            throw new IllegalArgumentException();
        }
    }

    public float[] getCell(int x, int y) {
        return getCell(x, y, new float[3]);
    }

    public String getDescription() {
        return "STAFramePlot";
    }


    public double getMinX() {
        if (Double.isNaN(minX)) 
            return 0;
        else 
            return minX;
    }


    public double getMaxX() {
        if (Double.isNaN(maxX)) 
            return frame.getWidth() * frame.getStixelWidth();
        else 
            return maxX;
    }


    public double getMinY() {
        if(Double.isNaN(minY))
            return 0;
        else
            return minY;
    }


    public double getMaxY() {
        if(Double.isNaN(maxY))
            return frame.getHeight() * frame.getStixelHeight();
        else
            return maxY;
    }
    
    /*
     * Sets the data bounds of the plot (not the view range).
     * 
     */
    public void setBounds(double minX, double maxX, double minY, double maxY) {
        this.minX = minX;
        this.maxX = maxX;
        this.minY = minY;
        this.maxY = maxY;
    }

}