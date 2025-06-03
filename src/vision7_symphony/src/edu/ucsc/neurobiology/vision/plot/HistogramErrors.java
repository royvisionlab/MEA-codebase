package edu.ucsc.neurobiology.vision.plot;

/**
 * This interface specifies the functionality any histogram MUST implement to
 * allow it to be drawn using the <tt>PlotPanel<tt> class.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public interface HistogramErrors
    extends PlotData {

    /**
     * Returns the error of the bin indexed by <tt>bin</tt>.
     */
    public double getBinError(int bin);

}
