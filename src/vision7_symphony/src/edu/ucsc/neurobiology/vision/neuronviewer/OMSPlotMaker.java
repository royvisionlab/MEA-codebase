package edu.ucsc.neurobiology.vision.neuronviewer;

import java.awt.*;
import javax.swing.tree.*;

import edu.ucsc.neurobiology.vision.analysis.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class OMSPlotMaker
    extends PlotMaker {
    OMSClassification omsClassification;
    public OMSPlotMaker() {
        super("OMS", CLASS_AND_NEURON_PLOT);
    }


    public void initialize(NeuronViewer viewer) {
        super.initialize(viewer);
        omsClassification = new OMSClassification(paramsFile);
    }


    public Component[] makePlot(IntegerList list, int plotType, TreePath classPath) {
        if (plotType == CLASS_PLOT) {
            return new Component[] {omsClassification.makeOMSHistgrams(list.toArray())};
        } else {
            return new Component[] {
                omsClassification.makeOMSScatter(list.toArray()[0]),
                omsClassification.makeOMSHistogram(list.toArray()[0])
            };

        }
    }


}
