package edu.ucsc.neurobiology.vision.neuronviewer;

import java.io.IOException;

import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreePath;

import edu.ucsc.neurobiology.vision.Vision;
import edu.ucsc.neurobiology.vision.anf.NeuronCleaning;
import edu.ucsc.neurobiology.vision.parameters.ParametersTable;
import edu.ucsc.neurobiology.vision.util.IntegerList;
import edu.ucsc.neurobiology.vision.util.SwingWorker;



/**
 * @author Matthew Grivich, The Salk Institute
 */
public class NeuronCleaningAction
extends CalculationAction {


    public NeuronCleaningAction() {
        super("Clean Neurons", CalculationAction.CLASS_ACTION);
    }


    public void doAction(final IntegerList list, final TreePath classTreePath) {


        SwingWorker worker = new SwingWorker() {
            public Object construct() {

                ParametersTable table = viewer.configuration.showDialog("Neuron Cleaning: Viewer", "Neuron Cleaning", viewer.mainFrame);
                if(table != null) {
                    int minSpikesCount = table.getIntParameter("Minimun Number of Spikes");
                    double maxContamination = table.getDoubleParameter("Maximum Contamination");
                    int coincidenceTime = table.getIntParameter("Coincidence Time");
                    double maxCorrelation = table.getDoubleParameter("Maximum Correlation");
                    try {

                        IntegerList[] lists = NeuronCleaning.removeLowSpikeAndContam(neuronFile, list, minSpikesCount, maxContamination);
                        int[] neurons = lists[0].toArray();
                        if (neurons != null) {
                            viewer.classificationHappened(neurons,
                                    classTreePath.pathByAddingChild(
                                            new DefaultMutableTreeNode("low spike and contam", true)));
                        }

                        neurons = ((IntegerList) (NeuronCleaning.removeCorrelated(neuronFile, lists[1], coincidenceTime, maxCorrelation)[0])).toArray();
                        if (neurons != null) {
                            viewer.classificationHappened(neurons,
                                    classTreePath.pathByAddingChild(
                                            new DefaultMutableTreeNode("coincidence", true)));
                        }
                    } catch(IOException e) {
                        e.printStackTrace();
                    }
                }
                Vision.getInstance().sendMessage("Done.");
                return null;
            }
        };
    worker.start();
    
}






}




