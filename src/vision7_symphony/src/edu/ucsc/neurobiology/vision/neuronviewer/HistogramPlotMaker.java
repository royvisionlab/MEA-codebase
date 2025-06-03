package edu.ucsc.neurobiology.vision.neuronviewer;

import java.awt.Component;

import javax.swing.JLabel;
import javax.swing.tree.TreePath;

import edu.ucsc.neurobiology.vision.math.CannotEvaluateException;
import edu.ucsc.neurobiology.vision.util.IntegerList;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class HistogramPlotMaker
    extends PlotMaker {

    public double xAxisType = 0;
    public double yAxisType = 0;

    public String expression = "";
    public double binWidth = 0.1;
    public double x1 = 0;
    public double x2 = 1;

    public HistogramPlotMaker() {
        super("Histogram Plot", CLASS_PLOT);
    }


    public Component makePlot(IntegerList list, int plotType, TreePath classPath) {
        if (list.size() == 0) {
            return new JLabel("No neurons selected");
        }
        try {
            return NeuronViewer.makeHistogram(viewer, expression,
                                              binWidth, x1, x2,
                                              classPath, list.toArray());
        } catch (CannotEvaluateException ex) {
            return new JLabel("Could not make histogram");
        }
    }

}
