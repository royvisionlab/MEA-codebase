package edu.ucsc.neurobiology.vision.math.fitting;

/**
 * This exception gets thrown by the <code>fit</code> routine of the
 * <code>Fitter</code> class when the data passed to it is invalid
 * (null data arrays, wrong array length and so on).
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class InvalidFitDataException
    extends RuntimeException {

    public InvalidFitDataException() {
        super();
    }


    public InvalidFitDataException(String message) {
        super(message);
    }

}
