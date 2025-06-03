package edu.ucsc.neurobiology.vision.neuronviewer;

import java.awt.*;
import javax.swing.tree.*;

import edu.ucsc.neurobiology.vision.analysis.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class FlashesPlotMaker
    extends PlotMaker {


    public FlashesPlotMaker() {
        super("Flashes", CLASS_AND_NEURON_PLOT);
    }


    public Component makePlot(IntegerList list, int plotType, TreePath classPath) {
        return FlashesClassification.makeFlashesPanel(paramsFile, list.toArray());
    }

}
