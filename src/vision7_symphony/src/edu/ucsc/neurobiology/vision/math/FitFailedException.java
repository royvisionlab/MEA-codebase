package edu.ucsc.neurobiology.vision.math;


/**
 * This exception gets thrown by the <code>fit</code> routine of the
 * <code>Fitter</code> class when the fitting process cannot be continued
 * due to some error (singular matrix, fit fails to converge and so on).
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class FitFailedException
    extends Exception {

    Exception cause;


    public FitFailedException() {
        super();
    }


    public FitFailedException(String message) {
        super(message);
    }


    public FitFailedException(String message, Exception cause) {
        super(message);
        this.cause = cause;
    }


    public void printStackTrace() {
        StackTraceElement[] st = this.getStackTrace();
        System.err.println("Exception occured: " + this.getClass().getName() + ": " +
                           this.getMessage());
        for (int i = 0; i < st.length; i++) {
            System.err.println(st[i]);
        }
        System.err.println("because: ");
        cause.printStackTrace();
    }

}
