package edu.ucsc.neurobiology.vision.math.fitting;

import edu.ucsc.neurobiology.vision.math.*;


/**
 * An abstract multidimensional function.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public abstract class Function {

    /**
     * This method MUST be properly implemented by subclasses to provide the value
     * of the model function.
     *
     * @param x the vector of coordinates
     * @return the value of the function
     */
    public abstract double getValue(final double[] x) throws CannotEvaluateException;

}
