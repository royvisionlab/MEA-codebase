package edu.ucsc.neurobiology.vision.math;

/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
@SuppressWarnings("serial")
public class CannotEvaluateException
    extends Exception {

    Exception cause;

    public CannotEvaluateException() {
    }


    public CannotEvaluateException(Exception cause) {
        this.cause = cause;
    }


    public CannotEvaluateException(String msg) {
        super(msg);
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
