package edu.ucsc.neurobiology.vision.dataview;

import java.awt.Color;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.io.IOException;

import javax.swing.JComponent;
import javax.swing.JScrollPane;

import edu.ucsc.neurobiology.vision.Vision;
import edu.ucsc.neurobiology.vision.io.RawDataFile;
import edu.ucsc.neurobiology.vision.io.RawDataHeader;
import edu.ucsc.neurobiology.vision.onlinediagnostics.RawDataScatterPlot;
import edu.ucsc.neurobiology.vision.parameters.IntegerParameter;
import edu.ucsc.neurobiology.vision.parameters.ParametersTable;
import edu.ucsc.neurobiology.vision.plot.PlotPanel;
import edu.ucsc.neurobiology.vision.plot.ScatterPlotStyle;
import edu.ucsc.neurobiology.vision.plot.SymbolType;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class RawDataView {
    private RawDataScatterPlot rawData1 , rawData2;
    private PlotPanel rawDataPanel;
    private ScatterPlotStyle style1, style2;
    private RawDataFile rawDataFile;


    public RawDataView(File file) throws IOException {
        rawDataFile = new RawDataFile(file);
        this.rawData1 = new RawDataScatterPlot(rawDataFile, 0);
            this.rawData2 = new RawDataScatterPlot(rawDataFile, 0);

        rawDataPanel = new PlotPanel();
        rawDataPanel.setLabels("Time (seconds)", "Amplitude (ADC Counts)");
        style1 = new ScatterPlotStyle("Trace 1",
                                      SymbolType.NONE, 2, Color.red, true, Color.red, 1);
        style2 = new ScatterPlotStyle("Trace 2",
            SymbolType.NONE, 2, Color.blue, true, Color.blue, 1);
        rawDataPanel.addData(rawData2, style2);
        rawDataPanel.addData(rawData1, style1);
        rawDataPanel.setRange(new double[]{0, 1, -500, 500});

        rawDataPanel.setLimits(new double[]{0, rawDataFile.getHeader().getNumberOfSamples()*
                rawDataFile.getHeader().getTimePerSample()/1000, -4000, 4000});
//        rawDataPanel.addPopupAction(new AbstractAction("Fourier Transform") {
//            public void actionPerformed(ActionEvent e) {
//                double T = 1.0 / 20000.0;
//                int N = 4096;
//                double omega = 1.0 / (N * T);
//                int n = rawData1.getPointCount() / N;
//                System.out.println("FFT " + N);
//
//                double[] averageSpectrum = new double[N / 2];
//                double[] p = new double[3];
//                for (int i = 0, k = 0; i < n; i++) {
//                    double[] a = new double[N];
//                    double[] b = new double[N];
//                    for (int j = 0; j < N; j++) {
//                        rawData1.getDataPoint(k, p);
//                        k++;
//                        a[j] = p[1];
//                    }
//                    FFT.fft(a, b, +1);
//                    for (int j = 0; j < N / 2; j++) {
//                        averageSpectrum[j] += Math.sqrt(a[j] * a[j] + b[j] * b[j]);
//                    }
//                }
//                double max = Double.NEGATIVE_INFINITY;
//                for (int i = 1; i < averageSpectrum.length; i++) {
//                    if (averageSpectrum[i] > max) {
//                        max = averageSpectrum[i];
//                    }
//                }
//                MathUtil.divide(averageSpectrum, max);
//
//                ScatterPlot ss = new ScatterPlot();
//                for (int j = 1; j < N / 2; j++) {
//                    ss.add(j * omega, averageSpectrum[j]);
//                }
//
//                PlotPanel pp = new PlotPanel();
//                pp.addData(ss, style1);
//                pp.autoscale();
//                pp.setYRange(0, 1);
//                pp.setLabels("Frequency (Hz)", "Amplitude");
//                PlotUtil.showData("Average Spectrum (" + n + ")", pp);
//            }
//        });

        RawDataHeader h = rawDataFile.getHeader();
        rawDataPanel.addToLegend("Header Size: " + h.getHeaderSize());
        rawDataPanel.addToLegend("nElectrodes: " + h.getNumberOfElectrodes());
        rawDataPanel.addToLegend("nSamples: " + h.getNumberOfSamples());
        rawDataPanel.addToLegend("arrayID: " + h.extractArrayID(h.getArrayID()) +
                                 " (" + h.extractElectrodeMapPart(h.getArrayID()) +
                                 "/" +
                                 h.extractElectrodeMapPartsCount(h.getArrayID()) +
                                 ")");
        rawDataPanel.addToLegend("Sampling: " + h.getSamplingFrequency());
//            rawDataPanel.addToLegend("Comment: " + h.getComment());

        Vision app = Vision.getInstance();
        app.createFrame(rawDataPanel, getController(), null,
                        rawData1.getFileName() + " - Raw Data");
    }


    /*
        private JMenu createMenu() {
            JMenu menu = new JMenu("Raw Data");
//        menu.setMnemonic(KeyEvent.VK_N);
            JMenuItem item = new JMenuItem("Save Red");
//        item.setMnemonic(KeyEvent.VK_O);
            item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    rawData1.printASCII();
                }
            });
            menu.add(item);
            item = new JMenuItem("Save Blue");
//        item.setMnemonic(KeyEvent.VK_S);
            item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    rawData2.printASCII();
                }
            });
            menu.add(item);
            return menu;
        }
     */

    private JComponent getController() {
        ParametersTable t = new ParametersTable();

        t.addParameter(
            new IntegerParameter("Electrode (red)", null, null, 0, 0,
                                 rawData1.getElectrodesCount() - 1),
            new PropertyChangeListener() {
            public void propertyChange(PropertyChangeEvent e) {
                int electrode = ( (IntegerParameter) e.getNewValue()).getValue();
                rawData1.setElectrode(electrode);
              
            }
            
            
        });
        
//        t.addParameter(
//                new IntegerParameter("Time (s)", null, null, 0, 0,
//                                     rawData1.getPointCount()/20000),
//                new PropertyChangeListener() {
//                public void propertyChange(PropertyChangeEvent e) {
//                    int electrode = ( (IntegerParameter) e.getNewValue()).getValue();
//                    rawData1.setElectrode(electrode);
//                  
//                }
//                
//                
//            });
        
                t.addParameter(
                    new IntegerParameter("Electrode (blue)", null, null, 0, 0, rawData1.getElectrodesCount() - 1),
                    new PropertyChangeListener() {
                        public void propertyChange(PropertyChangeEvent e) {
             int electrode = ((IntegerParameter)e.getNewValue()).getValue();
                            rawData2.setElectrode(electrode);
                        }
                    }
                );
         
        return new JScrollPane(t);
    }

}
