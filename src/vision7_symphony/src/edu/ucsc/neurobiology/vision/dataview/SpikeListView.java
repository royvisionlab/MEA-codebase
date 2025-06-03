package edu.ucsc.neurobiology.vision.dataview;

import java.beans.*;
import java.io.*;
import java.util.*;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.parameters.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class SpikeListView {
    private SpikeFile spikeFile;
    private SpikeTableModel tableModel;

    private int tMin = 0;
    private int tMax = 2000;
    private int electrode = 33;
    private int nElectrodes;


    public SpikeListView(SpikeFile spikeFile) throws IOException {
        this.spikeFile = spikeFile;
        nElectrodes = spikeFile.getNumberOfElectrodes();

        tableModel = new SpikeTableModel();
        JTable table = new JTable(tableModel);
        table.getColumnModel().getColumn(0).setHeaderValue("N");
        table.getColumnModel().getColumn(1).setHeaderValue("Electrode");
        table.getColumnModel().getColumn(2).setHeaderValue("Time(ms)");
        table.getColumnModel().getColumn(3).setHeaderValue("Amplitude(ADC)");
        tableModel.update();

        Vision app = Vision.getInstance();
        app.createFrame(new JScrollPane(table), getController(), null, "Spike List");
    }


    class SpikeTableModel
        extends AbstractTableModel {
        private ArrayList<Spike> spikes = new ArrayList<Spike> ();

        public void update() throws IOException {
            SpikeIterator i = spikeFile.iterator(tMin);
            spikes.clear();
            while (i.hasNext()) {
                Spike s = (Spike) i.next();
                if (s.electrode != electrode) {
                    continue;
                }
                if (s.time > tMax) {
                    break;
                }

                spikes.add( (Spike) s.clone());
            }
            this.fireTableChanged(new TableModelEvent(this));
        }


        public int getRowCount() {
            return spikes.size();
        }


        public int getColumnCount() {
            return 4;
        }


        public Object getValueAt(int rowIndex, int columnIndex) {
            switch (columnIndex) {
                case 0:
                    return "" + rowIndex;
                case 1:
                    return "" + spikes.get(rowIndex).electrode;
                case 2:
                    return "" + StringUtil.format(spikes.get(rowIndex).time);
                case 3:
                    return "" + StringUtil.format(spikes.get(rowIndex).amplitude);
                default:
                    return null;
            }
        }
    }


    private JComponent getController() {
        ParametersTable t = new ParametersTable();

        t.addParameter(
            new IntegerParameter("Electrode", null, null, electrode, 0, nElectrodes),
            new PropertyChangeListener() {
            public void propertyChange(PropertyChangeEvent event) {
                electrode = ( (IntegerParameter) event.getNewValue()).getValue();
                try {
                    tableModel.update();
                } catch (IOException e) {
                    Vision.reportException(e);
                }
            }
        }
        );

        t.addParameter(
            new RangeParameter("Range", null, null, tMin, tMax, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY),
            new PropertyChangeListener() {
            public void propertyChange(PropertyChangeEvent event) {
                RangeParameter p = (RangeParameter) event.getNewValue();
                tMin = (int) p.getMin();
                tMax = (int) p.getMax();
                try {
                    tableModel.update();
                } catch (IOException e) {
                    Vision.reportException(e);
                }
            }
        }
        );

        return new JScrollPane(t);
    }

}
