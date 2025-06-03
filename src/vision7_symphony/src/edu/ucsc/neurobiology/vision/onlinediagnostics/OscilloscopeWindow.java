package edu.ucsc.neurobiology.vision.onlinediagnostics;

import java.beans.*;

import java.awt.*;
import javax.swing.*;
import javax.swing.event.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.parameters.*;
import edu.ucsc.neurobiology.vision.plot.*;


/**
 * This class implements the Oscilloscope Window (which display the running raw data
 * in 1s chunks)
 * and the controller component used to controll the Oscilloscope.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class OscilloscopeWindow {
    private OscilloscopePlot[] oscilloscopePlot;
    private int nElectrodes;
    private JInternalFrame frame;


    public OscilloscopeWindow(
        final MultipleCompressedSampleInputStream sampleInputStream, RawDataHeader header,
        int n) {

        this.nElectrodes = header.getNumberOfElectrodes();

        JPanel p = new JPanel(new GridLayout(0, 1));
        oscilloscopePlot = new OscilloscopePlot[n];
        for (int i = 0; i < oscilloscopePlot.length; i++) {
            oscilloscopePlot[i] = new OscilloscopePlot("Oscilloscope", 20000, 20000);
            final OscilloscopePlot plot = oscilloscopePlot[i];

            sampleInputStream.addSampleListener(oscilloscopePlot[i]);

            PlotPanel oscilloscopePanel = new PlotPanel();
            oscilloscopePanel.setRange(0, 20000, -1000, 1000);
            oscilloscopePanel.setLabels("Time (samples)", "Amplitude (ADC Counts)");
            if (i != oscilloscopePlot.length - 1) {
                oscilloscopePanel.getXAxis().setShowlabel(false);
                oscilloscopePanel.getXAxis().setShowTicks(false);
            }
            oscilloscopePanel.addData(plot, new ScatterPlotStyle("Raw Data",
                SymbolType.NONE, 0, Color.black, true, Color.black, 1));

            p.add(oscilloscopePanel);
        }

        frame = Vision.getInstance().createFrame(
            p, getController(), null, "Oscilloscope");

        frame.addInternalFrameListener(new InternalFrameAdapter() {
            public void internalFrameClosed(InternalFrameEvent e) {
                for (int i = 0; i < oscilloscopePlot.length; i++) {
                    sampleInputStream.removeSampleListener(oscilloscopePlot[i]);
                }
            }
        });
    }


    public void dispose() {
        Vision.getInstance().removeFrame(frame);
    }


    private JComponent getController() {
        ParametersTable t = new ParametersTable();

        for (int i = 0; i < oscilloscopePlot.length; i++) {
            final OscilloscopePlot plot = oscilloscopePlot[i];
            t.addParameter(
                new IntegerParameter("Electrode", null, null, plot.getElectrode(),
                                     0, nElectrodes - 1),
                new PropertyChangeListener() {
                public void propertyChange(PropertyChangeEvent e) {
                    int electrode = ( (IntegerParameter) e.getNewValue()).getValue();
                    plot.setElectrode(electrode);
                }
            });
        }

        return new JScrollPane(t);
    }

}
