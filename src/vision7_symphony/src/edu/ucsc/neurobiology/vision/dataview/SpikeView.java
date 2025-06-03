package edu.ucsc.neurobiology.vision.dataview;

import java.beans.*;
import java.io.*;

import java.awt.*;
import java.awt.geom.*;
import javax.swing.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.electrodemap.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.parameters.*;
import edu.ucsc.neurobiology.vision.plot.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class SpikeView {
    private JComponent controller;

    // the collection of all spikes read from a certain spike file (*.spikes)
    private SpikeFile spikeFile;

    // the electrode map associated with the spikes
    private ElectrodeMap electrodeMap;

    // the spike amplitude is multiplied with this factor to get the ellipse radius
    private double scale = 0.1;

    // the time window for which the spikes will bw displayed
    private int tMin = 0;
    private int tMax = 2000;

    // the PlotPanel instance on which the spikes will be plotted
    private PlotPanel spikeView;


    /**
     *
     */
    public SpikeView(SpikeFile spikeFile) {
        this.spikeFile = spikeFile;

        electrodeMap = ElectrodeMapFactory.getElectrodeMap(spikeFile.getArrayID());

        spikeView = new PlotPanel();
        spikeView.setRange( -300, 300, -300, 300);
        spikeView.setLabels("x (\u03BCm)", "y (\u03BCm)");
        // add the data plotter that knows how to draw electrode maps
        spikeView.addDataPlotter(new ElectrodeMapPlotter());
        // add the actual electrode map
        spikeView.addData(electrodeMap, new ElectrodeMapStyle("Electrode Map"));
        // add the spike iterator, which is the collection of all the spikes
        // from tMin to tMax. Each spike is represented as an ellipse
        spikeView.addData(new SpikeIter(), new FunctionStyle("Spikes", Color.blue, 1));

        spikeView.addToLegend("Electrodes: " + spikeFile.getNumberOfElectrodes());
        spikeView.addToLegend("Spikes: " + spikeFile.getSpikesCount());
        spikeView.addToLegend("Threshold: " + spikeFile.getThreshold());
        spikeView.addToLegend("Mean Time Constant: " + spikeFile.getMeanTimeConstant());
        spikeView.addToLegend("First TTL: " + spikeFile.getFirstTTL());
        spikeView.addToLegend("TTLs: " + spikeFile.getTTLTimes().length);

        // make an internal window to display the spikeView
        Vision app = Vision.getInstance();
        controller = getController();
        app.createFrame(spikeView, controller, null,
                        spikeFile.getFileNameNoExtension() + " - Spikes");
    }


    /**
     * This class implements the DataIterator interface to convert the spike collection
     * into a collection of ellipses which will can be drawn on any PlotPanel.
     */
    class SpikeIter
        implements DataIterator {
        // the ellipse instance used to draw the spikes (see plot.ParametricEllipse)
        private ParametricEllipse ellipse = new ParametricEllipse(0, 0, 0, 0, 0);
        // this iterator is obtained from the SpikeFile class, and it iterates through
        // all the spikes from tMin onwards.
        private SpikeIterator i;
        // the iterator above returns instances of the Spike class
        private Spike s;
        private Point2D p;


        public String getDescription() {
            return "Spikes";
        }


        // this method gets called before the PlotPanel wants to repaint this plot
        public void reset() {
            try {
                i = spikeFile.iterator(tMin);
            } catch (IOException e) {
                Vision.reportException(e);
            }
        }


        // asks you if you have more data to plot (ellipses in this case)
        public boolean hasNextData() {
            s = (Spike) i.next();
            return i.hasNext() && (s.time < tMax);
        }


        // asks you to return the next data (ellipses in this case)
        public Object nextData() {
            p = electrodeMap.getPosition(s.electrode, p);
            double a = s.amplitude * scale;
            // update the ellipse params with this spike params and return it
            ellipse.setParameters(p.getX(), p.getY(), a, a, 0, 1, 1);
            return ellipse;
        }
    }


    /**
     * this method creates and returns a panel of components used to manipulate
     * the spike plot.
     */
    private JComponent getController() {
        ParametersTable t = new ParametersTable();

        t.addParameter(
            new DoubleParameter("Scale (\u00b5m/mV)", null, null, scale),
            new PropertyChangeListener() {
            public void propertyChange(PropertyChangeEvent e) {
                scale = ( (DoubleParameter) e.getNewValue()).getValue();
                spikeView.replotAllData();
            }
        });

        t.addParameter(
            new RangeParameter("Range(ms)", null, null, tMin, tMax, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY),
            new PropertyChangeListener() {
            public void propertyChange(PropertyChangeEvent e) {
                RangeParameter p = (RangeParameter) e.getNewValue();
                tMin = (int) p.getMin();
                tMax = (int) p.getMax();
                spikeView.replotAllData();
            }
        });

        return new JScrollPane(t);
    }

}
