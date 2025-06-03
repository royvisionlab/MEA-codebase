package edu.ucsc.neurobiology.vision.anf;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.electrodemap.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * This module plugs into the Raw Data Vision IO Pipeline (RawVIOP) and calculated the
 * covariance matrices on every electrode. At the end it save them to the CovarianceFile.
 *
 * @see CovarianceFile
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class CovariancesCalculator
    implements SampleListener, SpikeListener {

    private static final boolean DEBUG_FLOW = false;

    private final int nlPoints, nrPoints;
    private final int nPoints;
    private final int maxCovarianceSpikes;
    private final String covFileName;
    private final int nElectrodes, middleElectrode;
    private final float alpha;

    private boolean[] isDisconnected;
    private int maxAdjacentElectrodes;
    private int sampleIndex;
    private LinkedList<short[]> sampleList;
    private int[][] adjacent;
    private CovarianceCalculatorThread covThread1, covThread2;
    private short[] data;
    private float[] mean;
    private CovarianceMatrix[] covMatrix;

    private int[] nSpikes;
    private double[] step, remainder;
    VisionHeader header;


    public CovariancesCalculator(
        String covFileName, int nlPoints, int nrPoints, ElectrodeUsage electrodeUsage,
        int minCovarianceSpikes, int maxCovarianceSpikes, double minimizationError,
        VisionHeader spikeFileHeader, boolean doAlignment) {

        this.covFileName = covFileName;
        ElectrodeMap map = ElectrodeMapFactory.getElectrodeMap(spikeFileHeader.arrayID);
        this.nElectrodes = map.getNumberOfElectrodes();
        this.nlPoints = nlPoints;
        this.nrPoints = nrPoints;
        this.nPoints = nlPoints + nrPoints + 1;
        this.maxCovarianceSpikes = maxCovarianceSpikes;
        this.alpha = (float) (1.0 / (spikeFileHeader.meanTimeConstant * 20000.0));

        this.header = spikeFileHeader;
        header.electrodeUsage = electrodeUsage;
        header.minCovarianceSpikes = minCovarianceSpikes;
        header.maxCovarianceSpikes = maxCovarianceSpikes;
        header.nlPoints = nlPoints;
        header.nrPoints = nrPoints;
        header.minimizationError = minimizationError;
        header.covarianceType = VisionHeader.NOISE_UNWHITENED;

        //
        nSpikes = new int[nElectrodes];
        step = new double[nElectrodes];
        remainder = new double[nElectrodes];
        Arrays.fill(step, 1);
        Arrays.fill(remainder, 1);

        mean = new float[nElectrodes];

        isDisconnected = map.getDisconnectedElectrodesList();
        isDisconnected[0] = true;

        // create the adjacency matrix
        maxAdjacentElectrodes = 0;
        adjacent = new int[nElectrodes][];
        for (int electrode = 0; electrode < nElectrodes; electrode++) {
            adjacent[electrode] = map.getAdjacentsTo(electrode, electrodeUsage.ordinal());
            if (adjacent[electrode].length > maxAdjacentElectrodes) {
                maxAdjacentElectrodes = adjacent[electrode].length;
            }
        }

        data = new short[nPoints * maxAdjacentElectrodes + 1];

        // create the sample buffer
        sampleList = new LinkedList<short[]>();
        for (int i = 0; i < 5 * SpikeFinder.maxSpikeWidth; i++) {
        //for (int i = 0; i < 2 * SpikeFinder.maxSpikeWidth; i++) {
            sampleList.add(new short[nElectrodes]);
        }

        // create the CovarianceCalculator threads
        covMatrix = new CovarianceMatrix[nElectrodes];
        middleElectrode = nElectrodes / 2;
        covThread1 = new CovarianceCalculatorThread(
            nPoints, nlPoints, adjacent, 1, middleElectrode, minimizationError, doAlignment);
        covThread2 = new CovarianceCalculatorThread(
            nPoints, nlPoints, adjacent, middleElectrode + 1, nElectrodes - 1,
            minimizationError, doAlignment);
    }


    public void start() {
        covThread1.start();
        covThread2.start();
    }


//    int iii = 0;
    public void processSpike(Spike spike) {
        final int electrode = spike.electrode;

        if (isDisconnected[electrode] ||
            spike.time >= header.nSamples - nrPoints
            || spike.time < nlPoints) {
            return;
        }

        nSpikes[electrode]++;

//        remainder[electrode]++;
//        if (remainder[electrode] > step[electrode]) {
//            remainder[electrode] -= step[electrode];

        if ((nSpikes[electrode] < maxCovarianceSpikes) || (maxCovarianceSpikes < 0)) {

            // process the spike
            // assumes that the buffer is big enough.
            // will throw an exception if it is not.
            final int startTime =
                sampleList.size() - (sampleIndex - spike.time) - nlPoints;

            // copy the and string the spike data
            for (int i = 0; i < nPoints; i++) {
                short[] s = sampleList.get(startTime + i);
                for (int adj = 0; adj < adjacent[electrode].length; adj++) {
                    data[adj * nPoints + i] = s[adjacent[electrode][adj]];
                }
            }

//            if (iii == 0) {
//                PlotUtilities.showData("" + spike.electrode + ", " + spike.time,
//                                       new ScatterPlot(null, data, null), new ScatterPlotStyle());
//            }

            data[data.length - 1] = (short) electrode;
            if (sampleIndex >= nPoints) {
          //  	edu.ucsc.neurobiology.vision.plot.PlotPanel p = edu.ucsc.neurobiology.vision.plot.PlotUtil.showArray("Data", data);
           // 	try {
           // 	Thread.sleep(2000);
           // 	} catch(Exception e) {
          //  		
          //  	}
          //  	p.setVisible(false);

                if (electrode <= middleElectrode) {
                    covThread1.add(data);
                } else {
                    covThread2.add(data);
                }
            }
        }

//        iii++;
    }


    public void processSample(short[] sample) {
        // filter the sample and add it into the sample queue
        short[] s = sampleList.removeFirst();
        System.arraycopy(sample, 0, s, 0, nElectrodes);
        for (int i = 0; i < nElectrodes; i++) {
            mean[i] += alpha * (s[i] - mean[i]);
            s[i] -= mean[i];
        }
        sampleList.addLast(s);

        sampleIndex++;

        if (sampleIndex % 20000 == 0) { // recalculate the spike selection step
            double f = (double) header.nSamples / sampleIndex;
            for (int electrode = 0; electrode < nElectrodes; electrode++) {
                step[electrode] = f * nSpikes[electrode] / maxCovarianceSpikes;
            }
        }
    }


    public void finishSpikeProcessing() {
        if (DEBUG_FLOW) {
            System.out.println("finishProcessing() called on PCANFCovariances");
        }

        // finish the calculation
        covThread1.finish();
        covThread2.finish();
        while (covThread1.isAlive() || covThread2.isAlive()) {
            try {
                Thread.sleep(100);
            } catch (InterruptedException e) {}
        }

        if (DEBUG_FLOW) {
            System.out.println("The covariance threads died for " + covFileName + ".");
        }

        // save the covariance matrices
        try { //BEU
            CovarianceFile covFile = new CovarianceFile(covFileName, header);
            covThread1.save(covFile);
            covThread2.save(covFile);
            covFile.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    public void finishSampleProcessing() {
    }


    class CovarianceCalculatorThread
        extends Thread {

        private SpikeAligner aligner;
        private LinkedList<short[]> fullList, emptyList;
        private boolean finished = false;
        private int fromEl, toEl;
        private boolean doAlignment;


        public CovarianceCalculatorThread(
            int nPoints, int nlPoints, int[][] adjacent, int fromEl, int toEl,
            double minimizationError, boolean doAlignment) {

            aligner = new SpikeAligner(nPoints, nlPoints, adjacent, minimizationError);

            this.fromEl = fromEl;
            this.toEl = toEl;
            this.doAlignment = doAlignment;

            int maxAdjacents = -1;
            for (int electrode = 1; electrode < adjacent.length; electrode++) {
                if (adjacent[electrode].length > maxAdjacents) {
                    maxAdjacents = adjacent[electrode].length;
                }
            }

            for (int electrode = fromEl; electrode <= toEl; electrode++) {
                if (!isDisconnected[electrode]) {
                    covMatrix[electrode] = new CovarianceMatrix(
                        (nPoints - 2) * adjacent[electrode].length);
                }
            }

            emptyList = new LinkedList<short[]>();
            fullList = new LinkedList<short[]>();
            for (int i = 0; i < 10000; i++) {
                emptyList.add(new short[nPoints * maxAdjacents + 1]);
            }
        }


        public void add(short[] data) {
            short[] d;

            synchronized (emptyList) {
                while (emptyList.isEmpty()) {
                    try {
                        emptyList.wait();
                    } catch (InterruptedException e) {}
                }
                d = (short[]) emptyList.removeFirst();
            }

            System.arraycopy(data, 0, d, 0, data.length);

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
            short[] d;

            while (true) {
                synchronized (fullList) {
                    while (fullList.isEmpty()) {
                        if (finished) {
                            return;
                        }
                        try {
                            fullList.wait();
                        } catch (InterruptedException e) {}
                    }
                    d = (short[]) fullList.removeFirst();
                }

                final int electrode = (int) d[d.length - 1];
                try {
                    double[] alignedSpike = aligner.align(electrode, d, doAlignment);
                    covMatrix[electrode].addData(alignedSpike);
                } catch (CannotEvaluateException e) {
                    e.printStackTrace();
                }

                synchronized (emptyList) {
                    emptyList.addLast(d);
                    emptyList.notifyAll();
                }
            }
        }


        public CovarianceMatrix getCovariance(int electrode) {
            return covMatrix[electrode];
        }


        public void save(CovarianceFile dis) throws IOException {
            for (int el = fromEl; el <= toEl; el++) {
                dis.addCovariance(covMatrix[el], el);
            }
        }
    }


    public static void main(String[] args) {
        int N = 10;
        int n = 10;
        final double step = (double) N / n;
        double remainder = step;
        for (int i = 0; i < N; i++) {
            remainder++;
            if (remainder > step) {
                remainder -= step;
                System.out.println(i + " selected");
            }
        }
    }

}
