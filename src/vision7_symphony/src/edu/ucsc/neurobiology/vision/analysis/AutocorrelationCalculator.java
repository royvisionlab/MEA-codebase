package edu.ucsc.neurobiology.vision.analysis;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class AutocorrelationCalculator {
    static final double samplingPeriod = 0.05; //ms
    private final double timeWindow, timeBinning; //ms
    private LinkedList<Integer> list;
    private DoubleHistogram autocorrelation;

    public AutocorrelationCalculator(double timeWindow, double timeBinning) {
        this.list = new LinkedList<Integer>();
        list.clear();

        this.timeWindow = timeWindow / samplingPeriod;
        this.timeBinning = timeBinning / samplingPeriod;
//        this.t1 = t1;
//        this.t2 = t2;
//        this.timeBinning = timeBinning;

        this.autocorrelation =
            new DoubleHistogram("auto", 0, this.timeWindow, this.timeBinning);
    }


    public void addSpike(int time) {
        // calculate the autocorrelation
        if (list.isEmpty()) {
            list.addLast(new Integer(time));
        } else {
            int oldestTime = list.getFirst();

            while (time - oldestTime > timeWindow) {
                // remove oldest time
                list.removeFirst();
                if (list.isEmpty()) {
                    break;
                }

                for (Integer t : list) {
                    int dt = t - oldestTime;
                    autocorrelation.fill(dt, 1);
                }

                oldestTime = list.getFirst();
            }

            list.addLast(new Integer(time));
        }
    }


    public void finish() {
    }


    /**
     * Do not change the definition of the function. The time axis has to be in ms.
     *
     * @return DoubleHistogram
     */
    public DoubleHistogram getAutocorrelation() {
        DoubleHistogram h = new DoubleHistogram(
            "Auto", 0, timeWindow * samplingPeriod, timeBinning * samplingPeriod);

        for (int i = 0; i < h.getBinCount(); i++) {
            h.setBin(i, autocorrelation.getBin(i));
        }

        return h;
    }


    public static DoubleHistogram calculate(int[] times, double timeWindow, double timeBinning) {
        AutocorrelationCalculator c = new AutocorrelationCalculator(timeWindow, timeBinning);
        for (int t : times) c.addSpike(t);
        c.finish();
        return c.getAutocorrelation();
    }


    public static DoubleHistogram calculateISI(
        int[] time, int n, double timeWindow, double timeBinning) {

        DoubleHistogram h = new DoubleHistogram(
            "", 0, timeWindow / samplingPeriod,
            timeBinning / samplingPeriod);

        for (int i = 1; i < n; i++) {
            h.fill(time[i] - time[i - 1], 1);
        }

        DoubleHistogram h1 = new DoubleHistogram("Auto", 0, timeWindow, timeBinning);
        for (int i = 0; i < h.getBinCount(); i++) {
            h1.setBin(i, h.getBin(i));
        }

        return h1;
    }


    public static DoubleHistogram calculateISI(
        double[] time, int n, double timeWindow, double timeBinning) {

        DoubleHistogram h = new DoubleHistogram(
            "", 0, timeWindow / samplingPeriod,
            timeBinning / samplingPeriod);

        for (int i = 1; i < n; i++) {
            h.fill(time[i] - time[i - 1], 1);
        }

        DoubleHistogram h1 = new DoubleHistogram("Auto", 0, timeWindow, timeBinning);
        for (int i = 0; i < h.getBinCount(); i++) {
            h1.setBin(i, h.getBin(i));
        }

        return h1;
    }


    public static DoubleHistogram calculateISI(
        IntegerList time, double timeWindow, double timeBinning) {

        DoubleHistogram h = new DoubleHistogram(
            "", 0, timeWindow / samplingPeriod,
            timeBinning / samplingPeriod);

        int nSpikes = time.size();
        for (int i = 1; i < nSpikes; i++) {
            h.fill(time.get(i) - time.get(i - 1), 1);
        }

        DoubleHistogram h1 = new DoubleHistogram("Auto", 0, timeWindow, timeBinning);
        for (int i = 0; i < h.getBinCount(); i++) {
            h1.setBin(i, h.getBin(i));
        }

        return h1;
    }


    public static double getContamination(int id, NeuronFile neuronFile) throws
        IOException {
        int[] t = neuronFile.getSpikeTimes(id);
        return getContamination(t, neuronFile.getNumberOfSamples());
    }
    
    public static double getContamination(int[] spikeTimes, int nSamples) {
        DoubleHistogram auto = calculate(
                spikeTimes, 5, 1 / VisionParams.SAMPLES_PER_MILLISECOND);
        double N = spikeTimes.length;
        double Nb = auto.count(VisionParams.ACFT1, VisionParams.ACFT2);
        double T = nSamples / VisionParams.SAMPLES_PER_MILLISECOND;
        double dT = VisionParams.ACFT2 - VisionParams.ACFT1;
        return Nb * T / (dT * N * N);
    }


}
