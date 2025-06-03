package edu.ucsc.neurobiology.vision.plot;


/**
 * This interface contains general methods common for all kinds of data plots. Up
 * to now there is only one common method: <tt>getDescription()</tt> with should return
 * a human understandable description of the plot.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public interface PlotData {

    /**
     * Returns a human understandable desciption of the data plot.
     *
     * @return the description as a <tt>String</tt>.
     */
    public String getDescription();
}
