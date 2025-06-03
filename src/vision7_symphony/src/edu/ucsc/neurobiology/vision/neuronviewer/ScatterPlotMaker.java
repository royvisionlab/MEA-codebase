package edu.ucsc.neurobiology.vision.neuronviewer;

import java.awt.*;
import javax.swing.*;
import javax.swing.tree.*;

import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ScatterPlotMaker
    extends PlotMaker {

    public String expression1 = "";
    public String expression2 = "";

    public ScatterPlotMaker() {
        super("Scatter Plot", CLASS_PLOT);
    }


    public Component makePlot(IntegerList list, int plotType, TreePath classPath) {
        if (list.size() == 0) {
            return new JLabel("No neurons selected");
        }
        try {
            return viewer.makeScatterPlot(expression1, expression2, classPath, list.toArray());
        } catch (CannotEvaluateException ex) {
            return new JLabel("Could not make histogram");
        }
    }

}
