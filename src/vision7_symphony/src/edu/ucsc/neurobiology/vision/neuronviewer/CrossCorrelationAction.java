package edu.ucsc.neurobiology.vision.neuronviewer;

import java.io.IOException;

import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreePath;

import edu.ucsc.neurobiology.vision.dataview.CrossCorrelator;
import edu.ucsc.neurobiology.vision.plot.PlotUtil;
import edu.ucsc.neurobiology.vision.plot.ScatterPlot;
import edu.ucsc.neurobiology.vision.plot.ScatterPlotStyle;
import edu.ucsc.neurobiology.vision.util.IntegerList;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class CrossCorrelationAction
    extends CalculationAction {

    public CrossCorrelationAction() {
        super("Correlaton Analysis", CalculationAction.CLASS_ACTION);
    }


    public void doAction(final IntegerList idList, final TreePath classTreePath) {
        if (viewer.leftSelection instanceof String &&
            viewer.rightSelection instanceof String) {
            TreePath p1 = viewer.leftTree.getSelectionPath();
            TreePath p2 = viewer.rightTree.getSelectionPath();
            IntegerList list1 = new IntegerList();
            IntegerList list2 = new IntegerList();
            InteractiveTree.getNeuronsInClass(
                (DefaultMutableTreeNode) p1.getLastPathComponent(), list1, false);
            InteractiveTree.getNeuronsInClass(
                (DefaultMutableTreeNode) p2.getLastPathComponent(), list2, false);
            boolean sameClass = p1.equals(p2);

            ScatterPlot sp = new ScatterPlot();

            for (int i = 0; i < list1.size(); i++) {
                int id1 = list1.get(i);
                for (int j = sameClass ? i + 1 : 0; j < list2.size(); j++) {
                    int id2 = list2.get(j);
                    double x1 = paramsFile.getDoubleCell(id1, "x0");
                    double y1 = paramsFile.getDoubleCell(id1, "y0");
                    double x2 = paramsFile.getDoubleCell(id2, "x0");
                    double y2 = paramsFile.getDoubleCell(id2, "y0");
                    double d = Math.sqrt( (x2 - x1) * (x2 - x1) +
                                         (y2 - y1) * (y2 - y1));
                    try {
                        sp.add(d, CrossCorrelator.getCoupling(
                            neuronFile.getSpikeTimes(id1),
                            neuronFile.getSpikeTimes(id2), 5, 1000));
                    } catch (IOException ex) {
                        ex.printStackTrace();
                    }
                }
            }

            ScatterPlotStyle s = new ScatterPlotStyle();
            s.setSymbolSize(2);
            PlotUtil.showData(classTreePath.toString(), sp, s).pad();
        }

        /*
                 ParametersTable params =
            viewer.configuration.showDialog("Axon Speed", "Axon Speed", viewer);
                 if (params == null) {
            return;
                 }

                 try {
            AxonSpeed t = new AxonSpeed(imgFile);
            t.setAxonTemplate(params.getIntParameter("Template Neuron ID"),
                              params.getIntParameter("Template Electrode"));

            double minSpikeAmplitude = params.getDoubleParameter("Min Spike Amplitude");
            double maxChiSquared = params.getDoubleParameter("Min Chi Squared");
            int minElectrodes = params.getIntParameter("Min Electrodes");

            final HashMap<Integer, Double> valueMap = new HashMap();
            for (int i = 0; i < idList.size(); i++) {
                double speed = t.getSpeed(idList.get(i), idList.size() == 1,
         minSpikeAmplitude, maxChiSquared, minElectrodes);
                if (speed > 0) {
                    valueMap.put(idList.get(i), speed);
                }
            }

            if (idList.size() > 1) {
                PlotPanel p = NeuronViewer.makeHistogram(viewer, valueMap, 0.05, 0, 5,
                    classTreePath, "speed");
                PlotUtil.showData(InteractiveTree.pathToString(classTreePath) +
                                  " axon propagation speed", p);
                p.setLabels("Axon Propagation Speed (mm/ms)", "Number of Cells");
            }
                 } catch (Exception ex) {
            ex.printStackTrace();
                 }
         */
    }


}
