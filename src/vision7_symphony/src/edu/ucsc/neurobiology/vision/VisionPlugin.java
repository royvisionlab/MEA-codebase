package edu.ucsc.neurobiology.vision;

/**
 * Any class that implements this interface can serve as the entry point of a Vision plugin.
 * This class also must have a no-argument contuctor, used to initialize it.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public interface VisionPlugin {
    public String getPluginName();
}
