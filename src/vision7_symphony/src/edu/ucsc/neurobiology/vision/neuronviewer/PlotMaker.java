package edu.ucsc.neurobiology.vision.neuronviewer;

import java.util.List;

import javax.swing.AbstractAction;
import javax.swing.KeyStroke;
import javax.swing.tree.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.electrodemap.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * Abstract PlotMaker used in the NeuronViewer. PlotMakers can be of 3 types, CLASS_PLOT,
 * NEURON_PLOT and CLASS_AND_NEURON_PLOT. Non-abstract subclasses have to be declared
 * in the "Neuron Viewer Plots" group in config.xml to be displayable in NeuronViewer.
 * The class provides a lot of information needed to build useful plots (pointers to
 * the NeuronViewer, NeuronFile, Config, etc).
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public abstract class PlotMaker {
    public static final int CLASS_PLOT = 0;
    public static final int NEURON_PLOT = 1;
    public static final int CLASS_AND_NEURON_PLOT = 2;
    protected String name;
    protected int plotType;
    protected NeuronViewer viewer;
    protected ParametersFile paramsFile;
    protected NeuronFile neuronFile;
    protected STACollection staFile;
    protected PhysiologicalImagingFile imgFile;
    protected ElectrodeMap electrodeMap;
    public Config configuration;
    public boolean isSelected;


    public PlotMaker(String name, int plotType) {
        this.name = name;
        this.plotType = plotType;
    }


    /**
     * Called when NeuronViewer is opened by the NeuronViewer. Do not call otherwise.
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
     * Returns the name of the plot.
     * @return String
     */
    public String getName() {
        return name;
    }


    /**
     * Returns the plot type: CLASS_ACTION, NEURON_ACTION or CLASS_AND_NEURON_ACTION
     * @return int
     */
    public int getPlotType() {
        return plotType;
    }


    /**
     * Here the plot must be made and returned.
     *
     * @param list IntegerList : neuron ids
     * @param classPath TreePath : the selected folder
     * @return Component
     */
    public abstract Object makePlot(IntegerList list, int plotType, TreePath classPath);

    
    /**
     * Implementers should provide a list of keybindings they want added to a calling ancestor.
     * The implementer is a "builder" in the sense that the keybinding should already be built and
     * ready to be returned by the time this is called, thus getKeybindings does not take arguments.
     * 
     * Since we assume that the implementer is some PlotMaker, and the keybinding is a hotkey relevant
     * to the particular PlotMaker, the expectation is that the keybinding will be built in the course 
     * of a call to makePlot.
     * 
     * @author peterli
     *
     */
    interface KeybindingBuilder {
        List<Keybinding> getKeybindings();
        void clearKeybindings();
    }
    
    static class Keybinding {
        public Object handle;
        public KeyStroke key;
        public AbstractAction action;
    }
}