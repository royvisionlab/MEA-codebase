package edu.ucsc.neurobiology.vision.plot;


/**
 * An interface for defining 2D histograms. Implemented by DoubleHistogram2D.
 *
 * @see HistogramData
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public interface Histogram2DData
    extends PlotData {

    public int getBinCountX();


    public int getBinCountY();


//    public int getBinCount();

    public double getBinIntervalX();


    public double getBinIntervalY();


    public double getBin(int x, int y);


    // public void clear();

    public String getDescription();


    public double getMinX();


    public double getMaxX();


    public double getMinY();


    public double getMaxY();
}
