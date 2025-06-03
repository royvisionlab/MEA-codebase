package edu.ucsc.neurobiology.vision.plot;


/**
 * This interface is used to draw color plots and is used for STA displays.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public interface ColorPlotData
    extends PlotData {

    public int getRowsCount();


    public int getColumnsCount();


    public float[] getCell(int x, int y, float[] pixel);
    public float[] getCell(int x, int y);

    
    public double getMinX();


    public double getMaxX();


    public double getMinY();


    public double getMaxY();
}
