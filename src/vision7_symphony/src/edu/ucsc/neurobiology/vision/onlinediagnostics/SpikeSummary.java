package edu.ucsc.neurobiology.vision.onlinediagnostics;

import java.beans.*;
import java.util.*;

import javax.swing.*;
import javax.swing.event.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.electrodemap.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.parameters.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class SpikeSummary
    implements SpikeListener {

    public final static int HEXAGONAL = 0;
    public final static int RECTANGULAR = 1;

    final static int DEF_MAX_FREQ = 100;
    int maxSpikeFrequency = DEF_MAX_FREQ;
    static int DEF_REFRESH = 1000;

    final ElectrodeMap electrodeMap;
    final int nElectrodes;
    int numberOfSpikes[];
    int maxMagnatudes[];
    double startTime = 0; //beginning of current refresh
    int timeStep = DEF_REFRESH * 20; //length of refresh in milliseconds, converted to samples
    SpikeSummaryPlot spikeSummaryPlot;
    DetailedSpikeSummaryPlot detailedSpikeSummaryPlot;
    ArrayList<Spike> spikes = new ArrayList<Spike>(2000);

    JInternalFrame f1, f2;

//    Point2D.Double p0, p;
//    double a = 4*60;
//    double b = 2*a/Math.sqrt(3);
//    int mainElec = 65;


    public SpikeSummary(final SpikeFinder spikeFinder, ElectrodeMap electrodeMap) {
        this.electrodeMap = electrodeMap;
        this.nElectrodes = electrodeMap.getNumberOfElectrodes();
        numberOfSpikes = new int[nElectrodes];
        maxMagnatudes = new int[nElectrodes];
        spikeSummaryPlot = new SpikeSummaryPlot(electrodeMap);
        detailedSpikeSummaryPlot = new DetailedSpikeSummaryPlot(electrodeMap);

        f1 = Vision.getInstance().createFrame(
            spikeSummaryPlot, getController(), null, "Spike Summary");
        f2 = Vision.getInstance().createFrame(
            detailedSpikeSummaryPlot, null, null, "Detailed Spike Summary");

        InternalFrameListener ifl = new InternalFrameAdapter() {
            public void internalFrameClosed(InternalFrameEvent e) {
                try {
                    if (!f1.isClosed()) {
                        f1.setClosed(true);
                        spikeFinder.removeSpikeListener(SpikeSummary.this);
                    }
                    if (!f2.isClosed()) {
                        f2.setClosed(true);
                        spikeFinder.removeSpikeListener(SpikeSummary.this);
                    }
                } catch (PropertyVetoException ex) {
                    System.out.println("Close VETO");
                }
            }
        };
        f1.addInternalFrameListener(ifl);
        f2.addInternalFrameListener(ifl);

        spikeSummaryPlot.setElectrodeMap(HEXAGONAL);
        detailedSpikeSummaryPlot.setElectrodeMap(HEXAGONAL);

        spikeFinder.addSpikeListener(this);
//        p0 = new Point2D.Double(0, 0);
//        electrodeMap.getPosition(mainElec, p0);
//        p = new Point2D.Double(0, 0);
    }


    public void processSpike(Spike spike) {
        //This code assumes that the spikes are sorted by time.
        if (spike.time > (startTime + timeStep)) {
            spikeSummaryPlot.refresh(numberOfSpikes, maxMagnatudes,
                                     detailedSpikeSummaryPlot.getCenterElectrode());
            detailedSpikeSummaryPlot.refresh(spikes, startTime,
                                             spikeSummaryPlot.getCenterElectrode());
            startTime += timeStep; // output in 1 second blocks.  Time is measured in ms.

            spikes.clear();
            Arrays.fill(numberOfSpikes, 0);
            Arrays.fill(maxMagnatudes, 0);
        }

        // update the spike count and amplitude
        final int electrode = spike.electrode;
//        electrodeMap.getPosition(electrode, p);
//        if (Math.pow( (p.x - p0.x) / a, 2) + Math.pow( (p.y - p0.y) / b, 2) <= 1) {
        numberOfSpikes[electrode]++;
        maxMagnatudes[electrode] = (int) Math.max(maxMagnatudes[electrode],
                                                  spike.amplitude);
//        }

        //Copy spikes to holding array if they will be plotted
        if (detailedSpikeSummaryPlot.isElectrodeVisible(electrode)) {
            //I must copy spike, because otherwise I pass the reference, not the data and
            //the reference changes.
            spikes.add((Spike) spike.clone());

            /*     System.out.println("break2");
             for(int i=0; i<spikes.size(); i++) {
             if(((Spike) spikes.get(i)).getElectrode()==0) System.out.println(spikes.get(i));
             }*/
        }
    }


    public void finishSpikeProcessing() {}


    //user controlls for spikeSummaryPlot
    private JComponent getController() {
        ParametersTable t = new ParametersTable();
        EnumeratorParameter hexOrRect = new EnumeratorParameter("HexOrRect", null,
            "Electrade Map");
        hexOrRect.addChoice(SpikeSummary.HEXAGONAL, "Hexagonal");
        hexOrRect.addChoice(SpikeSummary.RECTANGULAR, "Rectangular");

        t.addParameter(
            hexOrRect,
            new PropertyChangeListener() {
            public void propertyChange(PropertyChangeEvent e) {
                int style = (int) ( (EnumeratorParameter) e.getNewValue()).getValue();
                spikeSummaryPlot.setElectrodeMap(style);
                detailedSpikeSummaryPlot.setElectrodeMap(style);
            }
        }
        );

        IntegerParameter maximumSpikeAmplitude = new IntegerParameter(
            "MaxSpikeAmplitude", null, null, spikeSummaryPlot.DEF_MAX_MAG_MAX, 1, 5000);

        t.addParameter(
            maximumSpikeAmplitude,
            new PropertyChangeListener() {
            public void propertyChange(PropertyChangeEvent e) {
                int newAmplitude = (int) ( (IntegerParameter) e.getNewValue()).getValue();
                spikeSummaryPlot.setMaxAmplitude(newAmplitude);
                detailedSpikeSummaryPlot.setMaxAmplitude(newAmplitude);
            }
        }
        );

        IntegerParameter maxFrequency = new IntegerParameter("maxFrequency (Hz)", null, null,
            DEF_MAX_FREQ, 1
            , 5000);
        spikeSummaryPlot.setSpikesMax( (DEF_MAX_FREQ * DEF_REFRESH * 20) / 20000);

        t.addParameter(
            maxFrequency,
            new PropertyChangeListener() {
            public void propertyChange(PropertyChangeEvent e) {
                maxSpikeFrequency = (int) ( (IntegerParameter) e.getNewValue()).
                                    getValue();
                spikeSummaryPlot.setSpikesMax( (maxSpikeFrequency * timeStep) / 20000);
            }
        }
        );

        IntegerParameter refreshRate = new IntegerParameter("Refresh Rate (ms)", null, null,
            DEF_REFRESH, 10, 100000);
        detailedSpikeSummaryPlot.setRefreshRate(timeStep);
        t.addParameter(
            refreshRate,
            new PropertyChangeListener() {
            public void propertyChange(PropertyChangeEvent e) {
                timeStep = 20 * ( (int) ( (IntegerParameter) e.getNewValue()).getValue());
                //*20 is because time is in samples, timeStep is in ms.
                 detailedSpikeSummaryPlot.setRefreshRate(timeStep);
                spikeSummaryPlot.setSpikesMax( (maxSpikeFrequency * timeStep) / 20000);

            }
        }
        );

        return new JScrollPane(t);
    }


    public void dispose() {
        Vision.getInstance().removeFrame(f1); //only f1 has to be elliminated
    }


}
