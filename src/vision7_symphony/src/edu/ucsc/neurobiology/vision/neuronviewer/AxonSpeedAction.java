package edu.ucsc.neurobiology.vision.neuronviewer;

import java.util.HashMap;

import javax.swing.tree.TreePath;

import edu.ucsc.neurobiology.vision.analysis.AxonSpeed;
import edu.ucsc.neurobiology.vision.parameters.ParametersTable;
import edu.ucsc.neurobiology.vision.plot.PlotPanel;
import edu.ucsc.neurobiology.vision.plot.PlotUtil;
import edu.ucsc.neurobiology.vision.util.IntegerList;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class AxonSpeedAction
    extends CalculationAction {

    public AxonSpeedAction() {
        super("Calculate Axon Speed", CalculationAction.CLASS_AND_NEURON_ACTION);
    }


    public void doAction(final IntegerList idList, final TreePath classTreePath) {
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

            final HashMap<Integer, Double> valueMap = new HashMap<Integer, Double>();
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
    }

}
