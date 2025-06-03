package edu.ucsc.neurobiology.vision.plot;

/**
 * This interface specifies the functionality any histogram MUST implement to
 * allow it to be drawn using the <tt>PlotPanel<tt> class.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public interface HistogramData
    extends PlotData {

    /**
     * Returns the number of bins in this histogram.
     */
    public int getBinCount();


    /**
     * Returns the bin indexed by <tt>bin</tt>.
     */
    public double getBin(int bin);


    /**
     * Returns the lower real space coordinate. At this position the histogram
     * will start.
     */
    public double getMin();


    /**
     * Returns the upper real space coordinate. At this position the histogram
     * will stop.
     */
    public double getMax();


    public double getBinInterval();

}
