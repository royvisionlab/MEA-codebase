package edu.ucsc.neurobiology.vision.gui;

/**
 * This interface is to be implemented by a graphics compenent displayable in a
 * Vision internal window. Vision will only close when all internal windows return true
 * when the canClose() method is called.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public interface Closable {

    /**
     * Should return true if Vision can exit and false if it cannot exit because
     * there is unsaved data (for example).
     */
    public boolean canClose();
}
