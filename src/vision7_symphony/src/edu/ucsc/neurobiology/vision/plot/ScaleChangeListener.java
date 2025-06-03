package edu.ucsc.neurobiology.vision.plot;


/**
 * Used by PlotPanel to inform interested objects that the scales changed.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public interface ScaleChangeListener {

    public void xScaleChanged(double xMin, double xMax);


    public void yScaleChanged(double yMin, double yMax);

}
