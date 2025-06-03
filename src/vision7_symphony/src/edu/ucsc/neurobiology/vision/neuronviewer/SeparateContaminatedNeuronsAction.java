package edu.ucsc.neurobiology.vision.neuronviewer;

import java.io.File;
import java.io.IOException;

import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreePath;

import edu.ucsc.neurobiology.vision.Vision;
import edu.ucsc.neurobiology.vision.analysis.AutocorrelationCalculator;
import edu.ucsc.neurobiology.vision.io.NeuronFile;
import edu.ucsc.neurobiology.vision.parameters.DoubleParameter;
import edu.ucsc.neurobiology.vision.parameters.ParametersTable;
import edu.ucsc.neurobiology.vision.parameters.StringParameter;
import edu.ucsc.neurobiology.vision.util.IntegerList;
import edu.ucsc.neurobiology.vision.util.SwingWorker;


/**
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class SeparateContaminatedNeuronsAction
    extends CalculationAction {

    public SeparateContaminatedNeuronsAction() {
        super("Separate Contaminated Neurons", CalculationAction.CLASS_ACTION);
    }


    public void doAction(final IntegerList list, final TreePath classTreePath) {
        SwingWorker worker = new SwingWorker() {
            double maxContamination;
            String dataSet;

            public Object construct() {
                ParametersTable table = configuration.showDialog(
                    "Separate Contaminated Neurons", "Separate Contaminated Neurons",
                    viewer.mainFrame);
                if (table == null) {
                    return null;
                }
                maxContamination = ( (DoubleParameter) table.getParameter(
                    "Max Contamination")).getValue();
                dataSet = ( (StringParameter) table.getParameter("Data Set")).
                          getValue();

                Vision.getInstance().sendMessage("Removing contaminated neurons...");
                try {
                    return separateContaminatedNeurons(list, maxContamination,
                        dataSet);
                } catch (IOException e) {
                    Vision.reportException(e);
                    return null;
                } finally {
                    Vision.getInstance().sendMessage("Done.");
                    Vision.getInstance().setProgress(0);
                }
            }


            public void finished() {
                int[] neurons = (int[]) getValue();
                if (neurons != null) {
                    viewer.classificationHappened(
                        neurons, classTreePath.pathByAddingChild(
                            new DefaultMutableTreeNode("contaminated", true)));
                }
            }
        };
        worker.start();
    }


    private int[] separateContaminatedNeurons(IntegerList list, double maxContamination,
                                              String dataSet) throws
        IOException {

        File ff = new File(viewer.filePathRoot);
        NeuronFile contamNeuronFile = new NeuronFile(ff.getParentFile().getParentFile().
            getAbsolutePath() + File.separator +
            dataSet + File.separator + dataSet +
            ".neurons");

        IntegerList idList = new IntegerList();

        for (int i = 0; i < list.size(); i++) {
            Vision.getInstance().setProgress(100 * i / list.size());

            int id = list.get(i);
            if (
                AutocorrelationCalculator.getContamination(id, contamNeuronFile) >
                maxContamination) {
                idList.add(id);
            }

        }

        return idList.toArray();
    }
}
