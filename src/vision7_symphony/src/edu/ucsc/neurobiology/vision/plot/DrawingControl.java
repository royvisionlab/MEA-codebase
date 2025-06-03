package edu.ucsc.neurobiology.vision.plot;

import java.io.*;


/**
 * Any graphical component that displays plots should have an instance of this class
 * to allow for plots that change in time.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public interface DrawingControl
    extends Serializable {

    /**
     * Called by the various plots to announce the change of data and the need
     * of redrawing.
     *
     * @param data the plot that heeds to be redrawn
     */
    public void updateNeeded(Object data);
}
