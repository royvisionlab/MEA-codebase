package edu.ucsc.neurobiology.vision.neuronviewer;

import java.util.*;

import static java.awt.Color.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.tree.*;

import edu.ucsc.neurobiology.vision.math.*;
import static edu.ucsc.neurobiology.vision.plot.SymbolType.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class AveragePlotMaker
    extends PlotMaker {


    public String expression = "GreenTimeCourse";

    private ScatterPlotStyle averageStyle = new ScatterPlotStyle(
        "Average", NONE, 0, black, true, black, 1);

    public AveragePlotMaker() {
        super("Average Plot", CLASS_PLOT);
    }


    public Component makePlot(IntegerList list, int plotType,
                              final TreePath classTreePath) {
        if (list.size() == 0) {
            return new JLabel("No neurons selected");
        }

        HashMap<Integer, double[]> v;
        try {
            v = paramsFile.evaluate(expression, list.toArray());
        } catch (CannotEvaluateException ex) {
            return new JLabel("Cannot evaluate this expression");
        }

        // average
        int size = v.get(list.get(0)).length;
        final double[] average = new double[size];
        int n = 0;
        for (double[] x : v.values()) {
            if (x != null) {
                MathUtil.add(average, x);
                n++;
            }
        }
        if (n == 0) {
            return new JLabel("No neurons selected");
        }
        MathUtil.divide(average, n);

        PlotPanel p = new PlotPanel();

        p.addPopupAction(new AbstractAction("Print Points") {
            public void actionPerformed(ActionEvent e) {
                for (int i = 0; i < average.length; i++) {
                    System.err.println(average[i]);
                }
            }
        });

        p.addData(new ScatterPlot(average), averageStyle);
        p.setLabels(expression, "");
        p.autoscale();

        return p;
    }


}
