package edu.ucsc.neurobiology.vision.analysis;

import javax.swing.tree.*;

import edu.ucsc.neurobiology.vision.util.*;
import edu.ucsc.neurobiology.vision.neuronviewer.*;

/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public abstract class CalculationAction {
    public static final int CLASS_ACTION = 0;
    public static final int NEURON_ACTION = 1;
    private String name;
    private int action;

    public NeuronViewer viewer;


    public CalculationAction(String name, int actionType) {
        this.name = name;
        this.action = actionType;
    }

    public abstract void doAction(IntegerList list, TreePath classPath);


    public String getName() {
        return name;
    }


    public int getAction() {
        return action;
    }
}
