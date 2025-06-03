package edu.ucsc.neurobiology.vision.plot;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ElectrodeMapStyle
    implements PlotStyle {

    private final String description;

    public ElectrodeMapStyle(String description) {
        this.description = description;
    }


    public String getDescription() {
        return description;
    }
}
