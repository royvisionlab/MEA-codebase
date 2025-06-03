package edu.ucsc.neurobiology.vision.dataview;

import java.io.*;
import java.util.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.electrodemap.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.plot.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class NeuronWaveformView
    extends JPanel {

    private PlotPanel panel;
    private ScatterPlotStyle style;
    private double[][] data;
    private int nElectrodes;
    int nSpikes;
    int neuronID;
    int nLeftPoints, nRightPoints, nPoints;
    double minAmplitude, maxAmplitude;
    NeuronFile neurons;
    RawDataFile rawData;
    JSpinner nLeftPointsBox, nRightPointsBox, nSpikesBox, neuronIDBox;
    JCheckBox subtractMeanBox, goodSpikesBox, microvoltBox;
    JTextField electrodesBox;
    boolean subtractMean, goodSpikes;
    int[] electrodes;


    public NeuronWaveformView(RawDataFile rawData, NeuronFile neurons) throws IOException {
        super(new BorderLayout());
        this.neurons = neurons;
        this.rawData = rawData;

        style = new ScatterPlotStyle();
        style.setSymbolType(SymbolType.NONE);
        style.setConnectingPoints(true);

        int arrayID = rawData.getHeader().getArrayID();
        ElectrodeMap electrodeMap = ElectrodeMapFactory.getElectrodeMap(arrayID);
        nElectrodes = electrodeMap.getNumberOfElectrodes();

        JPanel controlPanel = new JPanel(new BorderLayout());
        controlPanel.add(makeUpperPanel(), BorderLayout.NORTH);
        controlPanel.add(makeLowerPanel(), BorderLayout.SOUTH);
        this.add(controlPanel, BorderLayout.NORTH);

        panel = new PlotPanel();
        this.add(panel, BorderLayout.CENTER);

        Vision.getInstance().createFrame(this, null, null, "Waveforms");
    }


    private JPanel makeLowerPanel() {
        JPanel controlPanel = new JPanel(new GridLayout(1, 0));

        controlPanel.add(new JLabel("Neuron ID", JLabel.RIGHT));
        neuronIDBox = new JSpinner(new SpinnerNumberModel(0, 0, 100000, 1));
        controlPanel.add(neuronIDBox);

        controlPanel.add(new JLabel("Spikes", JLabel.RIGHT));
        nSpikesBox = new JSpinner(new SpinnerNumberModel(200, 0, 100000, 1));
        controlPanel.add(nSpikesBox);

        controlPanel.add(new JLabel("Points", JLabel.RIGHT));
        nLeftPointsBox = new JSpinner(new SpinnerNumberModel(10, 0, 500, 1));
        nRightPointsBox = new JSpinner(new SpinnerNumberModel(15, 0, 500, 1));
        controlPanel.add(nLeftPointsBox);
        controlPanel.add(nRightPointsBox);

        subtractMeanBox = new JCheckBox("Zero Mean", false);
        controlPanel.add(subtractMeanBox);

        goodSpikesBox = new JCheckBox("Good Spikes", false);
        controlPanel.add(goodSpikesBox);

        microvoltBox = new JCheckBox("\u03BCV", false);
        controlPanel.add(microvoltBox);

        return controlPanel;
    }


    private JPanel makeUpperPanel() {
        JPanel lowerPanel = new JPanel(new BorderLayout());

        lowerPanel.add(new JLabel("Electrodes"), BorderLayout.WEST);
        electrodesBox = new JTextField("1");
        lowerPanel.add(electrodesBox, BorderLayout.CENTER);

        JButton updateButton = new JButton("Do");
        updateButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                nLeftPoints = ( (Integer) nLeftPointsBox.getValue()).intValue();
                nRightPoints = ( (Integer) nRightPointsBox.getValue()).intValue();
                nPoints = nLeftPoints + nRightPoints + 1;
                style.setConnectionPeriod(nPoints);
                nSpikes = ( (Integer) nSpikesBox.getValue()).intValue();
                neuronID = ( (Integer) neuronIDBox.getValue()).intValue();
                subtractMean = subtractMeanBox.isSelected();
                goodSpikes = goodSpikesBox.isSelected();
                electrodes = getelectrodes();
                IOUtil.printArray(electrodes);
                update();
            }
        });
        lowerPanel.add(updateButton, BorderLayout.EAST);

        return lowerPanel;
    }


    private int[] getelectrodes() {
        ArrayList<String> list = new ArrayList<String>();
        StringTokenizer t = new StringTokenizer(
            electrodesBox.getText(), ",", false);
        while (t.hasMoreTokens()) {
            list.add(t.nextToken().trim());
        }
        int[] electrodes = new int[list.size()];
        for (int i = 0; i < electrodes.length; i++) {
            try {
                electrodes[i] = Integer.parseInt( (String) list.get(i));
            } catch (NumberFormatException e) {
                JOptionPane.showMessageDialog(
                    Vision.getInstance().getMainFrame(),
                    "Please correct: " + list.get(i) + " !",
                    "Error", JOptionPane.ERROR_MESSAGE);
                return null;
            }
        }
        return electrodes;
    }


    /*
        private int[] getelectrodes() {
            int el;
            try {
                el = Integer.parseInt(electrodesBox.getText());
            } catch (NumberFormatException e) {
                JOptionPane.showMessageDialog(
                    Vision.getInstance().getMainFrame(),
                    "Please correct: " + " !",
                    "Error", JOptionPane.ERROR_MESSAGE);
                return null;
            }
         ElectrodeMap map = ElectrodeMapFactory.getElectrodeMap(rawData.getHeader());
            IntegerList list = new IntegerList();
            list.add(el);
            for (int i = 0; i < map.getNumberOfElectrodes(); i++) {
                if (i != el && map.areAdjacent(i, el)) {
                    list.add(i);
                }
            }
            return list.toArray();
        }
     */

    private void update() {
        if (electrodes == null) {
            return;
        }

        boolean microvolt = microvoltBox.isSelected();
        double factor;
        if (microvolt) {
            panel.setLabels("Time (samples)", "Amplitude (\u03BCV)");
            factor = 0.7;
        } else {
            panel.setLabels("Time (samples)", "Amplitude (ADC Count)");
            factor = 1;
        }

        try {
            int[] time = neurons.getSpikeTimes(neuronID);

            // load in all the required data
            data = new double[nSpikes][nPoints * electrodes.length];
            short[] samples = new short[nPoints];

            for (int spikeIndex = 0; spikeIndex < nSpikes; spikeIndex++) {
                final int sampleNumber = time[spikeIndex] - nLeftPoints;
                if (sampleNumber < 0) {
                    continue;
                }

                // load the data for this particular electrode
                for (int el = 0; el < electrodes.length; el++) {
                    rawData.getData(electrodes[el], sampleNumber, samples);
                    for (int i = 0; i < nPoints; i++) {
                        data[spikeIndex][el * nPoints + i] = samples[i] * factor;
                    }
                }
            }
        } catch (IOException e) {
            Vision.reportException(e);
        }

        // create the supperposition plots
        ScatterPlot superposition = new ScatterPlot();
        for (int spikeIndex = 0; spikeIndex < nSpikes; spikeIndex++) {
            // add the waveform to the scatterplot
            for (int sample = 0; sample < nPoints * electrodes.length; sample++) {
                superposition.add(sample, data[spikeIndex][sample]);
            }
        }

        double[] bounds = PlotUtil.getPlotBounds(superposition);
        panel.setRange(bounds);
        panel.removeAllData();
        panel.addData(superposition, style);
    }

}
