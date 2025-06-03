package edu.ucsc.neurobiology.vision.neuronviewer;

import javax.swing.tree.*;

import edu.ucsc.neurobiology.vision.util.*;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class PrintCommentAction extends CalculationAction {
    
    public PrintCommentAction() {
        super("Print Comment", CalculationAction.NEURON_ACTION);
    }
    
    public void doAction(final IntegerList list, final TreePath classTreePath) {
        String comment = paramsFile.getStringCell(list.get(0), "comment");
        System.err.println(comment == null ? "no comment" : comment);
    }
    
}