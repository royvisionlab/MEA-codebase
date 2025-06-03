package edu.ucsc.neurobiology.vision.analysis;

import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedList;

import edu.ucsc.neurobiology.vision.analysis.Imaging.DataContainer;
import edu.ucsc.neurobiology.vision.util.WaveformCalculator;

/**
 * Abstracted out of vision.analysis.Imaging
 */
class EICalculator extends Thread {
    private static int BUFFER_LENGTH = 1000;
    
    private LinkedList<DataContainer> fullList, emptyList;
    private boolean finished = false;
    HashMap<Integer, WaveformCalculator> calc;


    public EICalculator(int nlPoints, int nrPoints, int[] idList, int nElectrodes) {
        calc = new HashMap<Integer, WaveformCalculator>();

        for (int i = 0; i < idList.length; i++) {
            calc.put(idList[i], new WaveformCalculator(nElectrodes, nlPoints, nrPoints));
        }

        int nPoints = nlPoints + nrPoints + 1;
        emptyList = new LinkedList<DataContainer>();
        fullList = new LinkedList<DataContainer>();
        for (int i = 0; i < BUFFER_LENGTH; i++) {
            DataContainer c = new Imaging.DataContainer();
            c.spike = new short[nPoints + 2][nElectrodes];
            emptyList.add(c);
        }
    }


    /**
     * It is important that the implementation does not alter the input data.
     * @param data
     */
    public void add(DataContainer data) {
        DataContainer d;

        synchronized (emptyList) {
            while (emptyList.isEmpty()) {
                try {
                    emptyList.wait();
                } catch (InterruptedException e) {}
            }
            d = emptyList.removeFirst();
        }

        d.neuronID = data.neuronID;
        d.time = data.time;
        for (int i = 0; i < d.spike.length; i++) {
            System.arraycopy(data.spike[i], 0, d.spike[i], 0, data.spike[i].length);
        }

        synchronized (fullList) {
            fullList.addLast(d);
            fullList.notifyAll();
        }
    }


    public void finish() {
        synchronized (fullList) {
            finished = true;
            fullList.notifyAll();
        }
    }


    public void run() {
        DataContainer d;

        while (true) {
            synchronized (fullList) {
                while (fullList.isEmpty()) {
                    if (finished) return;
                    try {
                        fullList.wait();
                    } catch (InterruptedException e) {}
                }
                d = fullList.removeFirst();
            }

            try {
                calc.get(d.neuronID).addSpike(d.spike);
            } catch (IOException ex) {
                ex.printStackTrace();
            }

            synchronized (emptyList) {
                emptyList.addLast(d);
                emptyList.notifyAll();
            }
        }
    }

}
