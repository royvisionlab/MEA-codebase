package edu.ucsc.neurobiology.vision.actions;

import java.io.*;

import java.awt.event.*;
import javax.swing.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.dataview.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.parameters.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class OpenNeuronWaveformViewer
    extends AbstractAction {

    public OpenNeuronWaveformViewer() {
        super("Neuron Waveform Viewer");
    }


    public void actionPerformed(ActionEvent event) {
        Vision app = Vision.getInstance();

        ParametersTable table = app.getConfig().showDialog(
            "NeuronWaveformViewer", "Open Neuorn Waveform Viewer", app.getMainFrame());
        if (table == null) {
            return;
        }

        try {
            String rawDataFileName = ( (FileParameter) table.getParameter(
                "Raw_Data_File")).getValue();
            String neuronFileName = ( (FileParameter) table.getParameter(
                "Neuron_File")).getValue();
            RawDataFile rawDataFile = new RawDataFile(new File(rawDataFileName));
            NeuronFile spikes = new NeuronFile(neuronFileName);

            new NeuronWaveformView(rawDataFile, spikes);
        } catch (IOException e) {
            Vision.reportException(e);
        }
    }

}
