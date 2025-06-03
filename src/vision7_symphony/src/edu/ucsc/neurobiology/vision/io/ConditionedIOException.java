package edu.ucsc.neurobiology.vision.io;

import java.io.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
@SuppressWarnings("serial")
public class ConditionedIOException
    extends IOException {

    public ConditionedIOException(Exception e) {
    }


    public ConditionedIOException(String msg, Exception e) {
    }
}
