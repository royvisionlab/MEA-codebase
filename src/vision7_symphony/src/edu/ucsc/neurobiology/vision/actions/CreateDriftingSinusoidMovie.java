package edu.ucsc.neurobiology.vision.actions;

import java.awt.event.*;
import javax.swing.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.parameters.*;
import edu.ucsc.neurobiology.vision.stimulus.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class CreateDriftingSinusoidMovie
    extends AbstractAction {


    public CreateDriftingSinusoidMovie() {
        super("Create Drifting Sinusoid Movie");
//        LinkedList a;
    }


    public void actionPerformed(ActionEvent event) {
        Vision app = Vision.getInstance();
        ParametersTable table = app.getConfig().showDialog(
            "CreateDriftingSinusoidMovie", "Create Drifting Sinusoid Movie",
            app.getMainFrame());
        if (table == null) {
            return;
        }

        int refreshInterval = table.getIntParameter("refreshInterval");
        int nFrames = table.getIntParameter("nFrames");
        int width = table.getIntParameter("width");
        int height = table.getIntParameter("height");
        double pixelSize = table.getDoubleParameter("pixelSize");
        double spatialPeriod = table.getDoubleParameter("spatialPeriod");
        double temporalPeriod = table.getDoubleParameter("temporalPeriod");

        new DriftingSinusoidMovie(
            nFrames, width, height, pixelSize, refreshInterval * 8.34,
            spatialPeriod, temporalPeriod, 0.48,
            100000);
    }

}
