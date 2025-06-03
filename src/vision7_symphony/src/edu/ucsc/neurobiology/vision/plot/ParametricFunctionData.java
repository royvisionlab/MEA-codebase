package edu.ucsc.neurobiology.vision.plot;


/**
 * This interface specifies the functionality any parametric function MUST implement
 * to allow it to be drawn using the <tt>PlotPanel<tt> class.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public interface ParametricFunctionData extends PlotData {

    /**
     * Returns the function point (x, y) fot the given parameter <tt>p</tt>.
     */
    public double[] getPointFor(double p);


    /**
     * Returns the lower limit of the parameter.
     */
    public double getMinParamValue();


    /**
     * Returns the upper limit of the parameter.
     */
    public double getMaxParamValue();
    
    /**
     * Set label
     */
    public void setLabel(String l);
    
    /**
     * Returns the label, or null if no label.
     */
    public String getLabel();
    
    /**
     * Point for label position
     */
    public double[] getLabelPosition();
}
