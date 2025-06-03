package edu.ucsc.neurobiology.vision.neuronviewer;

import java.awt.*;
import javax.swing.tree.*;

import edu.ucsc.neurobiology.vision.analysis.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class DriftingGratingsPlotMaker
    extends PlotMaker {

    public DriftingGratingsPlotMaker() {
        super("Drifting Gratings", CLASS_AND_NEURON_PLOT);
    }


    public void initialize(NeuronViewer viewer) {
        super.initialize(viewer);
    }


    public Component[] makePlot(IntegerList list, int plotType, TreePath classPath) {

        if (plotType == CLASS_PLOT) {

            return new Component[] {
                DriftingGratings.makeSinusoidsPanel(list.toArray(), paramsFile),
                DriftingGratings.makeClassAverageGratingPanel2D(
                    list.toArray(), paramsFile, true)
            };

        } else {

            return new Component[] {
                DriftingGratings.makeGratingPanel(list.toArray()[0], paramsFile, true),
                DriftingGratings.makeAverageGratingPanel(list.toArray()[0], paramsFile, true),
                DriftingGratings.makeAverageGratingPanel2D(list.toArray()[0], paramsFile, true)
            };

        }
    }


}
