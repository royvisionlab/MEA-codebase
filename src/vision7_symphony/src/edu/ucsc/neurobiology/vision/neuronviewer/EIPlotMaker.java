package edu.ucsc.neurobiology.vision.neuronviewer;

import static java.awt.Color.white;

import java.awt.Component;
import java.io.File;
import java.io.IOException;

import javax.swing.JLabel;
import javax.swing.tree.TreePath;

import edu.ucsc.neurobiology.vision.util.IntegerList;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class EIPlotMaker
    extends PlotMaker {
    
    //GUI/config.xml editable parameters
    public boolean showAmplitudeSpinner = true;       // hide amplitude spinner?
    public boolean showAnimationDelaySpinner = false; // hide animation framerate spinner?
    public boolean showFrameSlider = true;            // hide animation frame slider?
    public int defaultAnimationDelay = 40;            // Determines how fast movie plays
    public boolean negativeAmplitudes = true;	      // Show negative amplitudes as different color in EI movie playback
    public double minDisplayAmplitude = 0;            // min display amplitude
    public double defaultAmplitude = 4;               // value to start frame spinner at?

    public EIPlotMaker() {
        super("EI", NEURON_PLOT);
    }


    public void initialize(NeuronViewer viewer) {
        super.initialize(viewer);
    }


    public Component makePlot(IntegerList list, int plotType, TreePath classPath) {
        int neuronID = list.get(0);

        float[][][] img = null;
        if (imgFile != null) {
            try {
                img = imgFile.getImage(neuronID);
            } catch (IOException ex) {
            }
        }
        if (img != null) {
            return new PhysiologicalImagePanel(
                img, null, 2, electrodeMap, neuronFile.getElectrode(neuronID),
                new File(viewer.filePathRoot).getParent(), neuronID + ".gif", false, true, 1, showAmplitudeSpinner, showAnimationDelaySpinner, showFrameSlider, defaultAnimationDelay, negativeAmplitudes, defaultAmplitude, minDisplayAmplitude);
        } else {
            JLabel label = new JLabel("No EI for this neuron", JLabel.CENTER);
            label.setBackground(white);
            return label;
        }
    }
}