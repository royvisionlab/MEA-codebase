package edu.ucsc.neurobiology.vision.analysis;

import java.io.*;
import java.util.*;

import java.awt.*;
import javax.swing.*;

import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class DriftingContrastSinusoids {
    private HashMap<String, double[]> paramsMap;
    private int[] runList;
    private int[] nRuns;
    static int crRange = 500;
    static int ftRange = crRange / 2;
    NeuronFile neuronFile;

    final double dt = 8.34;
    final double T = dt / 1000;
    public final double[] contrasts = {0.03, 0.06, 0.12, 0.24, 0.48};
    public final int nContrasts = contrasts.length;
    final int N = 2 * 1024;
    final double omega = 1.0 / (N * T);

    double[][] a = new double[nContrasts][N];
    double[][] b = new double[nContrasts][N];
    int currentID;
    double timeToIgnore = 0;
    double contrast, temporalPeriod;


    public DriftingContrastSinusoids(String datasetFolder) throws IOException {
        String datasetName = new File(datasetFolder).getName();
        neuronFile = new NeuronFile(
            datasetFolder + File.separator + datasetName + ".neurons");
        loadStimulus(
            new File(datasetFolder).getParent() + File.separator + "stimuli" +
            File.separator + datasetName + ".txt");
    }


    public double getContrast() {
        return contrast;
    }


    public double getTemporalPeriod() {
        return temporalPeriod;
    }


    private void loadStimulus(String file) throws IOException {
        String[] fileName = {file};
        paramsMap = new HashMap<String, double[]>();

        InputStreamReader[] r = new InputStreamReader[fileName.length];
        StreamTokenizer[] st = new StreamTokenizer[fileName.length];
        IntegerList[] runs = new IntegerList[fileName.length];
        for (int i = 0; i < fileName.length; i++) {
            r[i] = new InputStreamReader(new FileInputStream(fileName[i]));
            runs[i] = new IntegerList();
            st[i] = new StreamTokenizer(r[i]);
            st[i].whitespaceChars('(', '(');
            st[i].whitespaceChars(')', ')');
            st[i].whitespaceChars(':', ':');
            st[i].whitespaceChars('#', '#');
            st[i].eolIsSignificant(true);
        }

        for (int fileIndex = 0; fileIndex < fileName.length; fileIndex++) {
            // read the stimulus type
            st[fileIndex].nextToken();
            if (st[fileIndex].ttype != StreamTokenizer.TT_WORD ||
                !st[fileIndex].sval.equals("type")) {
                throw new IOException("Missing type.");
            }
            st[fileIndex].nextToken();
            if (st[fileIndex].ttype != StreamTokenizer.TT_WORD) {
                throw new IOException("Wrong type.");
            }
//            String type = st[fileIndex].sval;

            // read the temporal period
            st[fileIndex].nextToken();
            if (st[fileIndex].ttype != StreamTokenizer.TT_WORD ||
                !st[fileIndex].sval.equals("temporal-period")) {
                throw new IOException("Missing temporal-period.");
            }
            st[fileIndex].nextToken();
            if (st[fileIndex].ttype != StreamTokenizer.TT_NUMBER) {
                throw new IOException("Wrong type.");
            }
            temporalPeriod = st[fileIndex].nval;

            // skip to the end of line
            do {
                st[fileIndex].nextToken();
            } while (st[fileIndex].ttype != StreamTokenizer.TT_EOL);

            // read the run lines
            while (true) {
                // read the contrast
                st[fileIndex].nextToken();
                if (st[fileIndex].ttype != StreamTokenizer.TT_WORD ||
                    !st[fileIndex].sval.equals("CONTRAST")) {
                    throw new IOException("Missing CONTRAST at line " +
                                          st[fileIndex].lineno());
                }
                st[fileIndex].nextToken();
                if (st[fileIndex].ttype != StreamTokenizer.TT_NUMBER) {
                    throw new IOException(
                        "Bad CONTRAST at line " + st[fileIndex].lineno());
                }
                double contrast = st[fileIndex].nval;

                // read the EOL
                st[fileIndex].nextToken();
                if (st[fileIndex].ttype != StreamTokenizer.TT_EOL) {
                    throw new IOException("Too many values at line " +
                                          st[fileIndex].lineno());
                }

                // figure out what kind of run we have
                int runNumber = -1;
                for (int i = 0; i < paramsMap.size(); i++) {
                    double[] params = (double[]) paramsMap.get("" + i);
                    if (contrast == params[0]) {
                        runNumber = i;
                        break;
                    }
                }

                // yet unknown run, register it
                if (runNumber == -1) {
                    runNumber = paramsMap.size();
                    double[] params = {contrast};
                    paramsMap.put("" + runNumber, params);
                    System.out.println(
                        "Run ID " + runNumber + ": " +
                        "Contrast " + contrast);
                }
                runs[fileIndex].add(runNumber);

                // see whether we reached the end of the file
                st[fileIndex].nextToken();
                if (st[fileIndex].ttype == StreamTokenizer.TT_EOF) {
                    break;
                } else {
                    st[fileIndex].pushBack();
                }
            }
        } // while

        r[0].close();
        runList = runs[0].toArray();

        nRuns = new int[paramsMap.size()];
        for (int i = 0; i < runList.length; i++) {
            nRuns[runList[i]]++;
        }
    }


    public DoubleHistogram[] getAverageResponse(double dt) throws IOException {
        int[] ttl = neuronFile.getTTLTimes();
        int[] times = neuronFile.getSpikeTimes(currentID);
        final double samplingRate = 20000.0; //Hz

        // check if the run numbers start with 0
        if (MathUtil.min(runList) != 0) {
            throw new IllegalStateException("Wrong run ID.");
        }
        final int nDifferentRuns = MathUtil.max(runList) + 1;

        // remove the intermediary ttl's
        IntegerList ttlList = new IntegerList();
        final int ttlsPerRepeat = 25;
        final int nRepeats = ttl.length / ttlsPerRepeat;

        double averageRunLength = 0;
        for (int run = 0; run < nRepeats; run++) {
            // add the TTL at which the sinusoid started
            ttlList.add(ttl[run * ttlsPerRepeat]);

            // calculate the average length of a sinusoid run plus the following gray time
            if (run != 0) {
                averageRunLength += ttlList.get(run) - ttlList.get(run - 1);
            }
        }
        // convert averageRunLength to seconds
        averageRunLength /= nRepeats * samplingRate;
//        System.out.println("averageRunLength " + averageRunLength);
        ttl = ttlList.toArray();

        // create the run histograms
        DoubleHistogram[] runHistograms = new DoubleHistogram[nDifferentRuns];
        for (int i = 0; i < nDifferentRuns; i++) {
            runHistograms[i] = new DoubleHistogram("" + i, 0, averageRunLength, dt);
        }

        // construct the spike stream
        SpikeStream stream = new SpikeStream(times, ttl);

        // process the spikes
        int currentRun = -1;
        int t;

        while ( (t = stream.getNext()) != Integer.MAX_VALUE) {
            if (t < 0) {
                currentRun++;
                if (currentRun >= runList.length) {
                    break;
                }
            } else {
                if (currentRun != -1) {
                    runHistograms[runList[currentRun]].fill(
                        (t - ttl[currentRun]) / samplingRate, 1);
                }
            }
        }

        // normalize the run histograms
        for (int i = 0; i < nDifferentRuns; i++) {
            runHistograms[i].scale(1.0 / (nRuns[i] * dt));
        }

        return runHistograms;
    }


    public DoubleHistogram[] getResponse(double dt, int nPeriods) throws IOException {
        DoubleHistogram[] h = getAverageResponse(dt);
        DoubleHistogram[] hh = new DoubleHistogram[h.length];
        double period = nPeriods * 30 * 8.34 / 1000;
        int binsPerPeriod = (int) Math.round(period / dt);
        System.out.println("period " + period);
        int shift = 0;
        int periodsToIgnore = (int) Math.ceil(timeToIgnore / period);
        for (int i = 0; i < h.length; i++) {
            hh[i] = new DoubleHistogram("", 0, period, dt);
            for (int j = periodsToIgnore; j < 24 * 30 * 8.34 / 1000 / dt; j++) {
                if (j >= shift) {
                    hh[i].fillBin( (int) ( (j - shift) % binsPerPeriod), h[i].getBin(j));
                }
            }
            hh[i].scale(1.0 / (24 / nPeriods - 1));
        }
        return hh;
    }


    /*
     public DoubleHistogram[] getResponse1(double dt, int nPeriods) throws IOException {
           DoubleHistogram[] h = getAverageReversingGratingsResponse(dt);
           DoubleHistogram[] hh = new DoubleHistogram[h.length];
           double period = nPeriods * 30 * 8.34 / 1000;
           int binsPerPeriod = (int) Math.round(period / dt);
           System.out.println("period " + period);
           int shift = 4;
           int periodsToIgnore = (int) Math.ceil(timeToIgnore / period);
           for (int i = 0; i < h.length; i++) {
               hh[i] = new DoubleHistogram("", 0, period, dt);
               for (int j = periodsToIgnore; j < 24 * 30 * 8.34 / 1000 / dt; j++) {
                   if (j >= shift) {
        hh[i].fillBin( (int) ( (j - shift) % binsPerPeriod), h[i].getBin(j));
                   }
               }
               hh[i].scale(1.0 / (24 / nPeriods - 1));
           }
           return hh;
       }
     */

    public String getRunInfo(int i) {
        double[] p = (double[]) paramsMap.get("" + i);
        return "contrast " + p[0];
    }


    private static final double getHarmonic(double frequency, double[] a, double omega) {
        int k = (int) Math.round(frequency / omega);
        double max = Double.NEGATIVE_INFINITY;
        for (int i = k - 4; i <= k + 4; i++) {
            if (a[i] > max) {
                max = a[i];
            }
        }
        return max;
    }


//    public void setCurrentNeuron(int id) throws IOException {
//        this.currentID = id;
//    }
//

    public void setCurrentNeuron(int id) throws IOException {
        this.currentID = id;
        for (int i = 0; i < a.length; i++) {
            Arrays.fill(a[i], 0);
            Arrays.fill(b[i], 0);
        }
        // load the reversing grating response histograms
        DoubleHistogram[] h = getAverageResponse(T);
        double[] noise = new double[h.length];
        for (int i = 0; i < h.length; i++) {
            // get the histogram
            double[] x = h[i].toArray();
            // subtract the mean
            double xm = MathUtil.mean(x);
            MathUtil.sub(x, xm);
            for (int j = (int) Math.round(1 / T); j < (int) Math.round(6 / T); j++) {
                a[i][j] = x[j];
            }
            FFT.fft(a[i], b[i], -1);
            for (int j = 0; j < N / 2; j++) {
                a[i][j] = Math.sqrt(a[i][j] * a[i][j] + b[i][j] * b[i][j]) * T;
            }
            noise[i] = MathUtil.mean(a[i], N / 4, N / 2);
        }
        /*
         for (int frequencyIndex = 0; frequencyIndex < periods.length; frequencyIndex++) {
                    Arrays.fill(maxHarmonic[frequencyIndex], Double.NEGATIVE_INFINITY);
                    Arrays.fill(meanHarmonic[frequencyIndex], 0);
                    for (int i = frequencyIndex * 4; i < (frequencyIndex + 1) * 4; i++) {
                        for (int j = 0; j < maxHarmonic.length; j++) {
                            double harmonic =
                                getHarmonic(4 * (j + 1), a[i], omega) - 1 * noise[i];
                            if (harmonic > maxHarmonic[frequencyIndex][j]) {
                                maxHarmonic[frequencyIndex][j] = harmonic;
                            }
                            meanHarmonic[frequencyIndex][j] += harmonic;
                        }
                    }
                    MathUtil.divide(meanHarmonic[frequencyIndex], 4);
                }
         */
    }


    /*
        public double getNonlinearityIndex(int frequencyIndex) {
            return meanHarmonic[frequencyIndex][1] / maxHarmonic[frequencyIndex][0];
//        double freq = 1 / (period[frequencyIndex] * 6.0 / 210.0);
        }
     */

    /*
        public void showOnePeriodPSTH() throws IOException {
            DoubleHistogram[] hh = getResponse(dt / 1000, 2);
            JFrame f = new JFrame();
            f.getContentPane().setLayout(new GridLayout(6, 4));
            HistogramStyle style = new HistogramStyle();
            for (int i = 0; i < hh.length; i++) {
                PlotPanel p = new PlotPanel();
                p.addData(hh[i], style);
                p.autoscale();
                p.setYRange(0, 100);
                p.setXAxisVisible(false);
                p.setYAxisVisible(false);
                ArrayList legend = new ArrayList();
                legend.add(getRunInfo(i));
                p.setAdditionalLegend(legend);
                f.add(p);
            }
            f.setBounds(0, 0, 1000, 800);
            f.setVisible(true);
        }
     */
    public void showPSTH() throws IOException {
        DoubleHistogram[] hh = getAverageResponse(dt / 1000);
        JFrame f = new JFrame();
        f.getContentPane().setLayout(new GridLayout(0, 1));
        HistogramStyle style = new HistogramStyle();
        for (int i = 0; i < hh.length; i++) {
            PlotPanel p = new PlotPanel();
            p.addData(hh[i], style);
            p.autoscale();
            p.setYRange(0, 150);
            p.setXAxisVisible(false);
            p.setYAxisVisible(false);
            p.addToLegend(getRunInfo(i));
            f.add(p);
        }
        f.setBounds(0, 0, 1000, 800);
        f.setVisible(true);
    }


    /*
             public void showPSTH(int periodIndex) throws IOException {
        DoubleHistogram[] hh = getAverageReversingGratingsResponse(dt / 1000);
        JFrame f = new JFrame();
        f.getContentPane().setLayout(new GridLayout(0, 1));
        f.setTitle("Neuron " + currentID + ", Period " + periods[periodIndex]);
        HistogramStyle style = new HistogramStyle();
        for (int i = periodIndex * nPhases; i < (periodIndex + 1) * nPhases; i++) {
            PlotPanel p = new PlotPanel();
            p.addData(hh[i], style);
            p.autoscale();
            p.setYRange(0, 150);
            p.setXAxisVisible(false);
            p.setYAxisVisible(false);
            ArrayList legend = new ArrayList();
            legend.add(getRunInfo(i));
            p.setAdditionalLegend(legend);
            f.add(p);
        }
        f.setBounds(0, 0, 1000, 800);
        f.setVisible(true);
             }
     */
    public void showSpectra() {
        ScatterPlotStyle style = new ScatterPlotStyle();
        style.setConnectingPoints(true);
        style.setSymbolType(SymbolType.NONE);
        JFrame f = new JFrame();
        f.setTitle(currentID + " spectra");
        f.getContentPane().setLayout(new GridLayout(0, 1));
        for (int i = 0; i < a.length; i++) {
            ScatterPlot s = new ScatterPlot();
            for (int j = 0; j < N / 2; j++) {
                s.add(j * omega, a[i][j]);
            }
            PlotPanel p = new PlotPanel();
            p.addData(s, style);
            p.autoscale();
            if (i < 8) {
                p.setYRange(0, 100);
            } else {
                p.setYRange(0, 10);
            }
            p.setXAxisVisible(false);
            p.setYAxisVisible(false);
            p.addToLegend(getRunInfo(i));
            f.add(p);
        }
        f.setBounds(0, 0, 1000, 800);
        f.setVisible(true);
    }


    /*
        public void showHarmonics(int[] id) throws IOException {
            ScatterPlotStyle st1 = new ScatterPlotStyle(ScatterPlotStyle.
                SYMBOLTYPE_FILLED_SQUARE, 4, Color.black, true, Color.black, 1);
            ScatterPlotStyle st2 = new ScatterPlotStyle(ScatterPlotStyle.
                SYMBOLTYPE_FILLED_SQUARE, 4, Color.red, true, Color.red, 1);
            PlotPanel p = new PlotPanel();
            for (int i = 0; i < id.length; i++) {
                setCurrentNeuron(id[i]);
                ScatterPlot tunungF1 = new ScatterPlot();
                ScatterPlot tunungF2 = new ScatterPlot();
//            double max1 = Double.NEGATIVE_INFINITY;
//            double max2 = Double.NEGATIVE_INFINITY;
//            for (int periodIndex = 0; periodIndex < period.length; periodIndex++) {
//                if (maxHarmonic[periodIndex][0] > max1) {
//                    max1 = maxHarmonic[periodIndex][0];
//                }
//                if (meanHarmonic[periodIndex][1] > max2) {
//                    max2 = meanHarmonic[periodIndex][1];
//                }
//            }
                double f = maxHarmonic[0][0];
                for (int periodIndex = 0; periodIndex < periods.length; periodIndex++) {
                    double freq = 1 / (periods[periodIndex] * 6.0 / 210.0);
                    tunungF1.add(freq, maxHarmonic[periodIndex][0] / f);
                    tunungF2.add(freq, meanHarmonic[periodIndex][1] / f);
                }
                p.addData(tunungF1, st1);
                p.addData(tunungF2, st2);
            }
            p.setRange(0, 6, 0, 5);
            p.setLabels("Spatial Frequency (cycles/deg)", "F1, F2 amplitude");
            PlotUtilities.showData("", p);
        }
     */
    /*
        public void showHarmonics() {
            ScatterPlot tunungF1 = new ScatterPlot();
            ScatterPlot tunungF2 = new ScatterPlot();
            ScatterPlotStyle st1 = new ScatterPlotStyle(ScatterPlotStyle.
                SYMBOLTYPE_FILLED_SQUARE, 4, Color.black, true, Color.black, 1);
            ScatterPlotStyle st2 = new ScatterPlotStyle(ScatterPlotStyle.
                SYMBOLTYPE_FILLED_SQUARE, 4, Color.red, true, Color.red, 1);
            for (int periodIndex = 0; periodIndex < periods.length; periodIndex++) {
                double freq = 1 / (periods[periodIndex] * 6.0 / 210.0);
                tunungF1.add(freq, maxHarmonic[periodIndex][0]);
                tunungF2.add(freq, meanHarmonic[periodIndex][1]);
            }
            PlotPanel p = new PlotPanel();
            p.addData(tunungF1, st1);
            p.addData(tunungF2, st2);
            p.setRange(0, 6, 0, 150);
            p.setLabels("Spatial Frequency (cycles/deg)", "F1, F2 amplitude");
            PlotUtilities.showData("", p);
        }
     */

}
