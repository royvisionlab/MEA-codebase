package edu.ucsc.neurobiology.vision.plot;

import edu.ucsc.neurobiology.vision.math.*;


/**
 * This interface specifies the functionality any analitycal function MUST implement
 * to allow it to be drawn using the <tt>PlotPanel<tt> class.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public interface FunctionError {

    /**
     * Returns the value of the function at the coordinate x (real space).
     */
    public Num getValueAndErrorAt(double x) throws CannotEvaluateException;
}
