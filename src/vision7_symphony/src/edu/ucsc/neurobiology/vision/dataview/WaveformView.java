package edu.ucsc.neurobiology.vision.dataview;

import java.io.*;

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
public class WaveformView
    extends JPanel {

    private PlotPanel panel;
    private ScatterPlotStyle style;
    private double[][] data;
    private int nElectrodes;
    int nSpikes;
    int nLeftPoints, nRightPoints, nPoints;
    double minAmplitude, maxAmplitude;
    SpikeFile spikes;
    RawDataFile rawData;
    JSpinner nLeftPointsBox, nRightPointsBox, nSpikesBox, minAmpBox, maxAmpBox;
    JCheckBox subtractMeanBox, goodSpikesBox, microvoltBox;
    JTextField electrodesBox;
    boolean subtractMean, goodSpikes;
    int[] electrodes;
    int arrayID;


    public WaveformView(RawDataFile rawData, SpikeFile spikes) throws IOException {
        super(new BorderLayout());
        this.spikes = spikes;
        this.rawData = rawData;
        arrayID = rawData.getHeader().getArrayID();

        style = new ScatterPlotStyle();
        style.setSymbolType(SymbolType.NONE);
        style.setConnectingPoints(true);

        ElectrodeMap electrodeMap = ElectrodeMapFactory.getElectrodeMap(spikes.getArrayID());
        nElectrodes = electrodeMap.getNumberOfElectrodes();

        JPanel controlPanel = new JPanel(new BorderLayout());
        controlPanel.add(makeUpperPanel(), BorderLayout.NORTH);
        controlPanel.add(makeLowerPanel(), BorderLayout.SOUTH);
        this.add(controlPanel, BorderLayout.NORTH);

        panel = new PlotPanel();
        this.add(panel, BorderLayout.CENTER);

        Vision.getInstance().createFrame(this, null, createMenu(), "Waveforms");
    }


    private JPanel makeLowerPanel() {
        JPanel controlPanel = new JPanel(new GridLayout(1, 0));

        controlPanel.add(new JLabel("Spikes", JLabel.RIGHT));
        nSpikesBox = new JSpinner(new SpinnerNumberModel(200, 0, 100000, 1));
        controlPanel.add(nSpikesBox);

        controlPanel.add(new JLabel("Points", JLabel.RIGHT));
        nLeftPointsBox = new JSpinner(new SpinnerNumberModel(10, 0, 500, 1));
        nRightPointsBox = new JSpinner(new SpinnerNumberModel(15, 0, 500, 1));
        controlPanel.add(nLeftPointsBox);
        controlPanel.add(nRightPointsBox);

        controlPanel.add(new JLabel("Amplitudes", JLabel.RIGHT));
        minAmpBox = new JSpinner(new SpinnerNumberModel(0, 0, 2048, 1));
        maxAmpBox = new JSpinner(new SpinnerNumberModel(1000, 0, 2048, 1));
        controlPanel.add(minAmpBox);
        controlPanel.add(maxAmpBox);

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
            public void actionPerformed(ActionEvent event) {
                nLeftPoints = ( (Integer) nLeftPointsBox.getValue()).intValue();
                nRightPoints = ( (Integer) nRightPointsBox.getValue()).intValue();
                nPoints = nLeftPoints + nRightPoints + 1;
                style.setConnectionPeriod(nPoints);
                nSpikes = ( (Integer) nSpikesBox.getValue()).intValue();
                minAmplitude = ( (Integer) minAmpBox.getValue()).intValue();
                maxAmplitude = ( (Integer) maxAmpBox.getValue()).intValue();
                subtractMean = subtractMeanBox.isSelected();
                goodSpikes = goodSpikesBox.isSelected();
                electrodes = getelectrodes();
                IOUtil.printArray(electrodes);
                try {
                    update();
                } catch (IOException e) {
                    Vision.reportException(e);
                }
            }
        });
        lowerPanel.add(updateButton, BorderLayout.EAST);

        return lowerPanel;
    }


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

        return ElectrodeMapFactory.getElectrodeMap(arrayID).getAdjacentsTo(el);
    }


    private JMenu[] createMenu() {
        JMenu menu = new JMenu("Waveforms");

        JMenuItem item = new JMenuItem("Save data");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                try {
                    PrintWriter pw = new PrintWriter(new FileWriter(
                        "C:\\MATLAB6p5\\work\\spikes.data"));

                    for (int s = 0; s < nSpikes; s++) {
                        for (int p = 0; p < nPoints * electrodes.length; p++) {
                            pw.print( (int) data[s][p] + " ");
                        }
                        pw.println();
                    }

                    pw.close();
                } catch (IOException e) {
                    Vision.reportException(e);
                }
            }
        });
        menu.add(item);

        return new JMenu[] {menu};
    }


    private void update() throws IOException {
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

        // load in all the required data
        data = new double[nSpikes][nPoints * electrodes.length];
        short[] samples = new short[nPoints];
        int plotIndex = 0, spikeIndex = -1;
        SpikeIterator iter = spikes.iterator();

        while (iter.hasNext() && (plotIndex < nSpikes)) {
            Spike s = (Spike) iter.next();
            if (s.electrode != electrodes[0]) {
                continue;
            }
            if (s.amplitude < minAmplitude || s.amplitude > maxAmplitude) {
                continue;
            }
            spikeIndex++;

            final int sampleNumber = s.time - nLeftPoints;
            if (sampleNumber < 0) {
                continue;
            }

            // load the data for this particular electrode
            for (int el = 0; el < electrodes.length; el++) {
                rawData.getData(electrodes[el], sampleNumber, samples);
//                    if (subtractMean) {
//                        for (int i = 0; i < nPoints; i++) {
//                            data[plotIndex][el * nPoints +
//                                i] = samples[i] - (s.amplitude + samples[nLeftPoints]);
//                        }
//                    } else {
                for (int i = 0; i < nPoints; i++) {
                    data[plotIndex][el * nPoints + i] = samples[i] * factor;
                }
//                    }
            }

            plotIndex++;
        }

        // create the supperposition plots
        ScatterPlot superposition = new ScatterPlot();
        for (plotIndex = 0; plotIndex < nSpikes; plotIndex++) {
            // add the waveform to the scatterplot
            for (int sample = 0; sample < nPoints * electrodes.length; sample++) {
                superposition.add(sample, data[plotIndex][sample]);
            }
        }

        double[] bounds = PlotUtil.getPlotBounds(superposition);
        panel.setRange(bounds);
        panel.removeAllData();
        panel.addData(superposition, style);
    }

}
