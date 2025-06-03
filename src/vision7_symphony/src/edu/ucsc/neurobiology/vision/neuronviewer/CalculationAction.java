package edu.ucsc.neurobiology.vision.neuronviewer;

import javax.swing.*;
import javax.swing.tree.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.electrodemap.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * Abstract action used in the NeuronViewer. Actions can be of 3 types, CLASS_ACTION,
 * NEURON_ACTION and CLASS_AND_NEURON_ACTION. Non-abstract subclasses have to be declared
 * in the "Neuron Viewer Actions" group in config.xml to be accessible in NeuronViewer.
 * The class provides a lot of information needed to build usefull actions (pointers to
 * the NeuronViewer, NeuronFile, Congig, etc).
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public abstract class CalculationAction {
    public static final int UNSUPPORTED_TYPE = -1;
    public static final int CLASS_ACTION = 0;
    public static final int NEURON_ACTION = 1;
    public static final int CLASS_AND_NEURON_ACTION = 2;

    protected NeuronViewer viewer;
    protected ParametersFile paramsFile;
    protected NeuronFile neuronFile;
    protected STACollection staFile;
    protected PhysiologicalImagingFile imgFile;
    protected ElectrodeMap electrodeMap;
    protected Config configuration;

    private String name;
    private int action;


    public CalculationAction(String name, int actionType) {
        this.name = name;
        this.action = actionType;
    }


    /**
     * Called when NeuronViewer is opened.
     *
     * @param viewer NeuronViewer
     */
    public void initialize(NeuronViewer viewer) {
        this.viewer = viewer;
        this.paramsFile = viewer.paramsFile;
        this.neuronFile = viewer.neuronFile;
        this.staFile = viewer.staCollection;
        this.imgFile = viewer.imgFile;
        this.electrodeMap = viewer.electrodeMap;
        this.configuration = viewer.configuration;
    }


    /**
     * The real action must be implemented here.
     *
     * @param list IntegerList
     * @param classPath TreePath
     */
    public abstract void doAction(IntegerList list, TreePath classPath);


    /**
     * Returns the name of the action.
     * @return String
     */
    public String getName() {
        return name;
    }


    /**
     * Returns the action type: CLASS_ACTION, NEURON_ACTION or CLASS_AND_NEURON_ACTION
     * @return int
     */
    public int getAction() {
        return action;
    }


    /**
     * A global keystroke can be mapped to the action by return one from after overwriting
     * this method.
     *
     * @return KeyStroke
     */
    public KeyStroke getKeystroke() {
        return null;
    }
}
