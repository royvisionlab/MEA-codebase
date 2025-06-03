package edu.ucsc.neurobiology.vision.plot;

import java.io.*;


/**
 * This interface specifies the functionality any scatter plot MUST implement to
 * allow it to be drawn using the <tt>PlotPanel<tt> class. A
 * scatter plot represents any number of freely distributed points.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public interface ScatterPlotData
    extends PlotData, Serializable {

    public int getPointCount();


    /**
     * This method returns a instance of <tt>ScatterPlotIterator</tt> class which
     * will iterate over the available data points.
     */
    public void getDataPoint(int i, double[] point);


    public boolean hasXErrors();


    public boolean hasYErrors();

}
