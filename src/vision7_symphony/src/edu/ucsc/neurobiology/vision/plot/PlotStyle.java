package edu.ucsc.neurobiology.vision.plot;

import java.io.*;


/**
 * This interface contains general methods common for all kinds of data plots. Up
 * to now there is only one common method: <tt>getDescription()</tt> with should return
 * a human understandable desciption of the plot.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public interface PlotStyle {

    /**
     * Returns a human understandable desciption of the data plot.
     *
     * @return the desciption as a <tt>String</tt>.
     */
    public String getDescription();

}
