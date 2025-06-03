package edu.ucsc.neurobiology.vision.neuronviewer;

import javax.swing.tree.*;

import edu.ucsc.neurobiology.vision.util.*;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class PrintIDsAction
    extends CalculationAction {

    public PrintIDsAction() {
        super("Print IDs", CalculationAction.CLASS_ACTION);
    }


    public void doAction(final IntegerList list, final TreePath classTreePath) {
        if (list.size() == 0) {
            System.err.println(InteractiveTree.pathToString(classTreePath) +
                               " - no neurons");
        } else {
            System.err.println(InteractiveTree.pathToString(classTreePath) + " - ID list");
            for (int i = 0; i < list.size(); i++) {
                System.err.println(list.get(i));
            }
        }
    }

}
