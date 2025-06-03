package edu.ucsc.neurobiology.vision.plot;


/**
 * An adapter with a default implementation of getDescription().
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public abstract class FunctionDataAdapter
    implements FunctionData {

    /**
     * @return String ""
     */
    public String getDescription() {
        return "";
    }
}
