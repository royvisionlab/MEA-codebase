package edu.ucsc.neurobiology.vision.neuronviewer;

import java.awt.Component;
import java.awt.GridLayout;
import java.awt.event.KeyEvent;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.KeyStroke;
import javax.swing.tree.TreePath;

import edu.ucsc.neurobiology.vision.math.CannotEvaluateException;
import edu.ucsc.neurobiology.vision.util.IntegerList;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ShowClassificationPlotsAction extends CalculationAction {

    public static String sizeExp = "2*((SigmaX*SigmaY)^0.5) *116";
    public static String tc1Exp =
        "pca(norm(RedTimeCourse#GreenTimeCourse#BlueTimeCourse), 0)";
    public static String tc2Exp =
        "pca(norm(RedTimeCourse#GreenTimeCourse#BlueTimeCourse), 1)";
    public static String auto1Exp = "pca(norm(Auto), 0)";
    public static String auto2Exp = "pca(norm(Auto), 1)";
    public static String f1f20 = "pca(norm(T1reversingF1#T1reversingF2), 0)";
    public static String f1f21 = "pca(norm(T1reversingF1#T1reversingF2), 1)";


    public ShowClassificationPlotsAction() {
        super("Show Classification Plots", CalculationAction.CLASS_ACTION);
    }


    public void doAction(final IntegerList list, final TreePath classTreePath) {
        int nColumns = 3;

        if (!NeuronViewer.isValidFolder(classTreePath)) {
            System.err.println("No class selected");
            return;
        }

        JPanel panel = new JPanel(new GridLayout(0, nColumns));
        try {
            panel.add(viewer.makeScatterPlot(tc1Exp, tc2Exp, "TF1", "TF2",
                                             classTreePath));
        } catch (CannotEvaluateException ex) {
        }
        try {
            panel.add(viewer.makeScatterPlot(f1f20, f1f21,
                                             "Reversing Gratings Response 1",
                                             "Reversing Gratings Response 2",
                                             classTreePath));
        } catch (CannotEvaluateException ex1) {
        }
        try {
            panel.add(viewer.makeScatterPlot(sizeExp, tc1Exp, "RF Diameter (\u03bcm)",
                                             "TF1",
                                             classTreePath));
        } catch (CannotEvaluateException ex2) {
        }

        try {
            panel.add(viewer.makeScatterPlot(sizeExp, f1f20, "RF Diameter (\u03bcm)",
                                             "Reversing Gratings Response 1",
                                             classTreePath));
        } catch (CannotEvaluateException ex3) {
        }
        try {
            panel.add(viewer.makeScatterPlot(tc1Exp, f1f20, "TF1",
                                             "Reversing Gratings Response 1",
                                             classTreePath));
        } catch (CannotEvaluateException ex4) {
        }
        try {
            Component p = viewer.makeHistogram(viewer, sizeExp, 0.01, classTreePath);
            if (p != null) {
                panel.add(p);
            }
        } catch (CannotEvaluateException ex5) {
        }

        try {
            panel.add(viewer.makeScatterPlot(auto1Exp, auto2Exp, "Autocorrelation 1",
                                             "Autocorrelation 2", classTreePath));
        } catch (CannotEvaluateException ex6) {
        }
//                panel.add(makeScatterPlot(auto0, f1f20, "Autocorrelation 1", "Reversing Gratings Response 1", classTreePath));
        try {
            panel.add(viewer.makeScatterPlot(sizeExp, auto1Exp,
                                             "RF Diameter (\u03bcm)",
                                             "Autocorrelation 1", classTreePath));
        } catch (CannotEvaluateException ex7) {
        }
        try {
            panel.add(viewer.makeScatterPlot(tc1Exp, auto1Exp, "TF1",
                                             "Autocorrelation 1", classTreePath));
        } catch (CannotEvaluateException ex8) {
        }

        JFrame frame = new JFrame();
        frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        frame.getContentPane().add(panel);
        frame.setBounds(0, 0, 1200, 1000);
        frame.setVisible(true);
    }


    public KeyStroke getKeystroke() {
        return KeyStroke.getKeyStroke(KeyEvent.VK_F9, KeyEvent.CTRL_MASK + KeyEvent.SHIFT_MASK);
    }
}
