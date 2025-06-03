package edu.ucsc.neurobiology.vision.neuronviewer;

import java.io.*;

import java.awt.event.*;
import javax.swing.*;
import javax.swing.tree.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.math.fitting.*;
import edu.ucsc.neurobiology.vision.parameters.*;
import edu.ucsc.neurobiology.vision.stimulus.*;
import edu.ucsc.neurobiology.vision.util.*;
import edu.ucsc.neurobiology.vision.util.SwingWorker;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class SeparateWeakTimecoursesAction
    extends CalculationAction {

    public SeparateWeakTimecoursesAction() {
        super("Separate Weak Timecourses", CalculationAction.CLASS_ACTION);
    }


    public void doAction(final IntegerList list, final TreePath classTreePath) {
        SwingWorker worker = new SwingWorker() {
            double chiSquared;

            public Object construct() {
                ParametersTable table = viewer.configuration.showDialog(
                    "Separate Bad Neurons", "Separate Bad Neurons", viewer.mainFrame);
                if (table == null) {
                    return null;
                }
                chiSquared = ( (DoubleParameter) table.getParameter(
                    "Chi Squared")).getValue();

                Vision.getInstance().sendMessage("Removing bad neurons...");
                try {
                    return separateBadNeurons(list, chiSquared);
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
                    viewer.classificationHappened(neurons,
                                                  classTreePath.pathByAddingChild(
                        new DefaultMutableTreeNode("weak " + chiSquared, true)));
                }
            }
        };
        worker.start();
    }


    private int[] separateBadNeurons(IntegerList list, double chiSquared) throws
        IOException {

        // separate weak timecourses
        IntegerList noTCList = new IntegerList();
        double[][] x = new double[1][viewer.staCollection.getSTADepth()];
        for (int k = 0; k < viewer.staCollection.getSTADepth(); k++) {
            x[0][k] = k;
        }

        for (int i = 0; i < list.size(); i++) {
            Vision.getInstance().setProgress(100 * i / list.size());

            int id = list.get(i);
            STA sta = viewer.staCollection.getSTA(id);
            double[][] y = sta.getTimeFilters(3);
            double[][] err = sta.getTimeFiltersError();
            if (y == null) {
                noTCList.add(id);
            } else {
                boolean bad = true;
                for (int c = 0; c < 3; c++) {
                    Fitter fit = new Fitter();
                    Linear1DFunction f = new Linear1DFunction(0, 0);
                    try {
                        fit.fit(f, x, y[c], err[c]);
                    } catch (FitFailedException e) {
                        f.setChiSquared(Double.POSITIVE_INFINITY);
                    }
                    if (f.getChiSquared() > chiSquared) {
                        bad = false;
                    }
                }
                if (bad) {
                    noTCList.add(id);
                }
            }
        }

        return noTCList.toArray();
    }


    public KeyStroke getKeystroke() {
        return KeyStroke.getKeyStroke(KeyEvent.VK_C,
                                      KeyEvent.CTRL_MASK + KeyEvent.SHIFT_MASK);
    }
}
