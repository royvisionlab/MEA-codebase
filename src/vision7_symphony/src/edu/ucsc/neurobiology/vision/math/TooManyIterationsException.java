package edu.ucsc.neurobiology.vision.math;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class TooManyIterationsException
    extends CannotEvaluateException {

    public TooManyIterationsException(String msg) {
        super(msg);
    }

}
