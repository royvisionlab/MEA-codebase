package edu.ucsc.neurobiology.vision.util;

import java.io.*;

import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.math.*;


/**
 * An average waveform calculator that works by using the spike waveform data
 * provided by the caller.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class WaveformCalculator {
    private final int nLeftPoints, nRightPoints, nPoints, nElectrodes;
//    private final short[][] samples;
    private final float[][] average, error;
    private int nSpikes;
//    private final int nSamples;
    private short[] spike;
    private UniformSpline spline;


    public WaveformCalculator(int nElectrodes, int nLeftPoints, int nRightPoints) {
        this.nElectrodes = nElectrodes;
        this.nLeftPoints = nLeftPoints;
        this.nRightPoints = nRightPoints;
        this.nPoints = nLeftPoints + nRightPoints + 1;
//        this.nSamples = rawData.getHeader().getNumberOfSamples();
        average = new float[nElectrodes][nPoints];
        error = new float[nElectrodes][nPoints];

        spike = new short[nPoints + 2];
        spline = new UniformSpline(nPoints + 2);

        reset();

//        System.err.println("nPoints " + nPoints);
    }


    public short[][] createCompatibleBuffer() {
        return new short[nPoints + 2][nElectrodes];
    }


    public void reset() {
        nSpikes = 0;
        for (int electrode = 0; electrode < nElectrodes; electrode++) {
            for (int i = 0; i < nPoints; i++) {
                average[electrode][i] = 0;
                error[electrode][i] = 0;
            }
        }
    }


    /*
        public void addSpike(double time) throws IOException {
            final int startSample = (int) Math.round(time) - nLeftPoints - 1;
            if (startSample > 0 && startSample < nSamples - (samples.length)) {
                rawData.getData(startSample, samples);

                double t0 = time - startSample - nLeftPoints;

                // allign the spikes
                for (int electrode = 0; electrode < nElectrodes; electrode++) {
                    for (int i = 0; i < nPoints + 2; i++) {
                        spike[i] = samples[i][electrode];
                    }

                    spline.reSpline(spike);
                    for (int i = 0; i < nPoints; i++) {
                        double v = spline.getValueAt(i + t0);
                        average[electrode][i] += v;
                        error[electrode][i] += v * v;
                    }
                }

                n++;
            }
        }
     */

    public void addSpike(double time, short[][] samples) throws IOException {
        final int startSample = (int) Math.round(time) - nLeftPoints - 1;
        double t0 = time - startSample - nLeftPoints - 1;

        // align the spikes
        for (int electrode = 0; electrode < nElectrodes; electrode++) {
            for (int i = 0; i < nPoints + 2; i++)
                spike[i] = samples[i][electrode];

            spline.reSpline(spike);
            for (int i = 0; i < nPoints; i++) {
                double v = spline.getValueAt(i + t0);
                average[electrode][i] += v;
                error[electrode][i] += v * v;
            }
        }

        nSpikes++;
    }


    public void addSpike(short[][] samples) throws IOException {
        for (int electrode = 0; electrode < nElectrodes; electrode++) {
            for (int t = 0; t < nPoints; t++) {
                average[electrode][t] += samples[t][electrode];
                error[electrode][t] += samples[t][electrode] * samples[t][electrode];
            }
        }
        nSpikes++;
    }
    
    
    public void finish() {
        for (int electrode = 0; electrode < nElectrodes; electrode++) {
            for (int i = 0; i < nPoints; i++) {
                average[electrode][i] /= nSpikes;
                error[electrode][i] = (float) Math.sqrt(
                    (error[electrode][i] -
                     nSpikes * average[electrode][i] * average[electrode][i]) / nSpikes);
            }
        }
    }


    public int getCount() {
        return nSpikes;
    }


    public float[][] getAverage() {
        return average;
    }


    public float[][] getError() {
        return error;
    }
    
    public int getNSpikes() {
        return nSpikes;
    }


    public float[][][] getImage() {
        return new float[][][] {average, error};
    }


    public static double[][] getOverlap(
        int[] electrodes, int[] t, int nSpikes, int nlPoints, int nrPoints,
        RawDataFile rawData) throws IOException {

        final int nPoints = nlPoints + nrPoints + 1;
        short[][] samples = new short[nPoints][
                            rawData.getHeader().getNumberOfElectrodes()];
        double[][] overlap = new double[nSpikes][electrodes.length * nPoints];

        for (int s = 0; s < nSpikes; s++) {
            final int sampleNumber = t[s] - nlPoints;
            if (sampleNumber > 0) {
                rawData.getData(sampleNumber, samples);
                for (int el = 0, index = 0; el < electrodes.length; el++) {
                    for (int i = 0; i < nPoints; i++) {
                        overlap[s][index++] = samples[i][electrodes[el]];
                    }
                }
            }
        }

        return overlap;
    }


    public static float[][][] calculateEI(
        double[] times, int nSpikes, int nEISpikes,
        int nLeftPoints, int nRightPoints, RawDataFile rawDatFile) throws IOException {

        DoubleList noBurstingT = new DoubleList(10000, 2);
        for (int i = 1; i < nSpikes - 1; i++) {
            if ( (times[i] - times[i - 1] > nLeftPoints + 20) &&
                (times[i + 1] - times[i] > nRightPoints + 20)) {
                noBurstingT.add(times[i]);
            }
        }
        final int NN = noBurstingT.size();
        BinaryRandom r = new BinaryRandom( (double) Math.min(nEISpikes, NN) / NN);
        WaveformCalculator average = new WaveformCalculator(rawDatFile.getHeader().
            getNumberOfElectrodes(), nLeftPoints, nRightPoints);
        int nSamples = rawDatFile.getHeader().getNumberOfSamples();
        int nElectrodes = rawDatFile.getHeader().getNumberOfElectrodes();

        short[][] samples = new short[nElectrodes][nLeftPoints + nRightPoints + 2];
        for (int i = 0; i < NN; i++) {
            if (r.random()) {
                double time = noBurstingT.get(i);
                final int startSample = (int) Math.round(time) - nLeftPoints - 1;
                if (startSample > 0 && startSample < nSamples - (samples.length)) {
                    rawDatFile.getData(startSample, samples);
                    average.addSpike(time, samples);
                }
            }
        }
        average.finish();

        return average.getImage();
    }


    public static float[][] calculateEI(
        double[] times, int nSpikes, int nEISpikes, int nLeftPoints, int nRightPoints,
        short[] data) {

        DoubleList noBurstingT = new DoubleList(10000, 2);
        for (int i = 1; i < nSpikes - 1; i++) {
            if ( (times[i] - times[i - 1] > nLeftPoints + 20) &&
                (times[i + 1] - times[i] > nRightPoints + 20)) {
                noBurstingT.add(times[i]);
            }
        }
        final int NN = noBurstingT.size();
        BinaryRandom r = new BinaryRandom( (double) Math.min(nEISpikes, NN) / NN);

        int nPoints = nLeftPoints + nRightPoints + 1;
        float[] average = new float[nPoints];
        float[] error = new float[nPoints];
        short[] spike = new short[nPoints + 2];
        UniformSpline spline = new UniformSpline(nPoints + 2);
        int n = 0;

        for (int sIndex = 0; sIndex < NN; sIndex++) {
            if (r.random()) {
                double time = noBurstingT.get(sIndex);
                final int startSample = (int) Math.round(time) - nLeftPoints - 1;
                if (startSample > 0) {
                    // calculate the exact time of the spike
                    System.arraycopy(data, startSample, spike, 0, spike.length);
                    spline.reSpline(spike);
                    double t0 = time - startSample - nLeftPoints;

                    for (int i = 0; i < nPoints; i++) {
                        double v = spline.getValueAt(i + t0);
                        average[i] += v;
                        error[i] += v * v;
                    }

                    n++;
                }
            }
        }

        for (int i = 0; i < average.length; i++) {
            average[i] /= n;
            error[i] = (float) Math.sqrt(
                (error[i] - n * average[i] * average[i]) / n);
        }

        // remove the means
        double a = 0;
        for (int i = 0; i < 10; i++) {
            a += average[i];
        }
        a /= 10;
        for (int i = 0; i < average.length; i++) {
            average[i] -= a;
        }

        return new float[][] {average, error};
    }
}
