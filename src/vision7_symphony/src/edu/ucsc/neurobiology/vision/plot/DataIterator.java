package edu.ucsc.neurobiology.vision.plot;

import java.io.*;


/**
 * The <tt>DataIterator</tt> interface can be implemented by a data source to allow
 * for drawing of series of plots of the same type.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public interface DataIterator
    extends PlotData, Serializable {

    /**
     * Resets the state of the iterator. Called every time before drawing.
     */
    public void reset();


    /**
     * Returns true if there still are plots to draw, false otherwise.
     */
    public boolean hasNextData();


    /**
     * Returns the next plot to be drawn.
     */
    public Object nextData();

}
