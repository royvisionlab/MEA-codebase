package edu.ucsc.neurobiology.vision.analysis;

import java.io.*;
import static java.lang.Math.*;
import java.util.*;

import static java.awt.Color.*;
import java.awt.*;
import javax.swing.*;

import static edu.ucsc.neurobiology.vision.io.IOUtil.*;
import edu.ucsc.neurobiology.vision.io.*;
import static edu.ucsc.neurobiology.vision.math.MathUtil.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.math.fitting.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.plot.HistogramStyle.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class ReversingGratings {
    // the elements of double[] are {spatialPeriod, temporalPeriod, spatialPhase}
    private HashMap<Integer, double[]> paramsMap;
    private int[] runIDList;
    private int[] nRuns;
    private NeuronFile neuronFile;
    private final int N;
    private double timeBinning, omega;
    private double[][] real, img;
    private int currentID;
    private double contrast, orientation;
    private double stimulusLength_fr;
    public double[] spatialPeriods, temporalPeriods, phases;
    public int nSpatialPeriods, nTemporalPeriods, nPhases;
    private double mmPerPixel = .0058;
    private String experimentName, datasetName;
    private int index1, index2;

    public static String xLabel = "Spatial Frequency (cyc/mm)";
    public static String yLabel = "<html>F<sub>1</sub> , F<sub>2</sub> (spikes/s)";
    static Color[] colors = {black, red, blue, green, magenta, orange};
    static String[] colors1 = {"black", "red", "blue", "green", "blue", "orange"};
    public static ScatterPlotStyle[] style1 = new ScatterPlotStyle[3];
    public static ScatterPlotStyle[] style2 = new ScatterPlotStyle[3];
    public static ScatterPlotStyle[] style1Thick = new ScatterPlotStyle[3];
    public static ScatterPlotStyle[] style2Thick = new ScatterPlotStyle[3];
    static {
        for (int tPeriod = 0, colorIndex = 0; tPeriod < 3; tPeriod++) {
            style1[tPeriod] = new ScatterPlotStyle("Fodd " + tPeriod,
                SymbolType.DISK, 1, colors[colorIndex], true,
                colors[colorIndex++], 0.25f);
            style2[tPeriod] = new ScatterPlotStyle("Feven " + tPeriod,
                SymbolType.DISK, 1, colors[colorIndex], true, colors[colorIndex++], 0.25f);
        }

        for (int tPeriod = 0, colorIndex = 0; tPeriod < 3; tPeriod++) {
            style1Thick[tPeriod] = new ScatterPlotStyle("Fodd " + tPeriod,
                SymbolType.NONE, 2, Color.black, true, Color.black, 1);
            style2Thick[tPeriod] = new ScatterPlotStyle("Feven " + tPeriod,
                SymbolType.NONE, 2, Color.red, true, Color.red, 1);
        }
        style1Thick[1].setDashPattern("3, 3");
        style2Thick[1].setDashPattern("3, 3");
        style1Thick[2].setDashPattern("1, 1");
        style2Thick[2].setDashPattern("1, 1");
    }


    public ReversingGratings(String datasetFolder, int nBinsPerPeriod) throws IOException {
        experimentName = new File(datasetFolder).getParentFile().getName();
        datasetName = new File(datasetFolder).getName();
        neuronFile = new NeuronFile(
            datasetFolder + File.separator + datasetName + ".neurons");
        loadStimulus(new File(datasetFolder).getParent() + File.separator + "stimuli" +
                     File.separator + datasetName + ".txt");

        double T = min(temporalPeriods) / 120; // temporal period in seconds
        this.timeBinning = T / nBinsPerPeriod;
        double nBins = (stimulusLength_fr / min(temporalPeriods) * nBinsPerPeriod);
        this.N = (int) Math.pow(2, Math.ceil(log2(nBins)));
//        System.err.println(N);
        this.omega = 1.0 / (N * timeBinning);

        real = new double[nPhases * nSpatialPeriods * nTemporalPeriods][N];
        img = new double[nPhases * nSpatialPeriods * nTemporalPeriods][N];
    }


    public void setRepeatSelection(int index1, int index2) {
        this.index1 = index1;
        this.index2 = index2;
    }


    public double getContrast() {
        return contrast;
    }


    private void loadStimulus(String fileName) throws IOException {
        paramsMap = new HashMap();
        InputStreamReader r = new InputStreamReader(new FileInputStream(fileName));
        StreamTokenizer st = new StreamTokenizer(r);
        IntegerList runs = new IntegerList();
        st.whitespaceChars('(', '(');
        st.whitespaceChars(')', ')');
        st.whitespaceChars(':', ':');
        st.whitespaceChars('#', '#');
        st.eolIsSignificant(false);

        // read the stimulus type
        st.nextToken();
        if (st.ttype != st.TT_WORD || !st.sval.toUpperCase().equals("TYPE")) {
            throw new IOException("Missing type.");
        }
        st.nextToken();
        if (st.ttype != st.TT_WORD) {
            throw new IOException("Wrong type.");
        }
        String type = st.sval;

        // read the rgb
        st.nextToken();
        if (st.ttype != st.TT_WORD || !st.sval.toUpperCase().equals("RGB")) {
            throw new IOException("Missing rgb.");
        }
        for (int i = 0; i < 3; i++) {
            st.nextToken();
            if (st.ttype != st.TT_NUMBER) {
                throw new IOException("Wrong type.");
            }
            contrast = st.nval;
        }

        orientation = readNumber(st, "ORIENTATION");
        stimulusLength_fr = readNumber(st, "FRAMES");
        
        double xStart = readNumber(st, "X-START");
        double xEnd = readNumber(st, "X-END");
        double yStart = readNumber(st, "Y-START");
        double yEND = readNumber(st, "Y-END");

        // skip to the end of line
  //      do {
  //          st.nextToken();
  //      } while (st.ttype != st.TT_EOL);
        
        DoubleList spatialPeriodsList = new DoubleList();
        DoubleList temporalPeriodsList = new DoubleList();
        DoubleList phaseList = new DoubleList();
        // read the run lines
        while (true) {
            double spatialPeriod = readNumber(st, "SPATIAL-PERIOD");
            double temporalPeriod = readNumber(st, "TEMPORAL-PERIOD");
            double spatialPhase = readNumber(st, "SPATIAL-PHASE");

            if (!spatialPeriodsList.contains(spatialPeriod)) {
                spatialPeriodsList.add(spatialPeriod);
            }
            if (!temporalPeriodsList.contains(temporalPeriod)) {
                temporalPeriodsList.add(temporalPeriod);
            }
            if (!phaseList.contains(spatialPhase / spatialPeriod)) {
                phaseList.add(spatialPhase / spatialPeriod);
            }

            // read the EOL
   //         st.nextToken();
   //         if (st.ttype != st.TT_EOL) {
   //             throw new IOException("Too many values at line " + st.lineno());
   //         }

            // figure out what kind of run we have
            int runNumber = -1;
            for (int i = 0; i < paramsMap.size(); i++) {
                double[] params = paramsMap.get(new Integer(i));
                if (spatialPeriod == params[0] && temporalPeriod == params[1] &&
                    spatialPhase == params[2]) {
                    runNumber = i;
                    break;
                }
            }

            // yet unknown run, register it
            if (runNumber == -1) {
                runNumber = paramsMap.size();
                double[] params = {spatialPeriod, temporalPeriod, spatialPhase};
                paramsMap.put(new Integer(runNumber), params);
            }
            runs.add(runNumber);

            // see whether we reached the end of the file
            st.nextToken();
            if (st.ttype == st.TT_EOF) {
                break;
            } else {
                st.pushBack();
            }
        }

        double[] _spatialPeriods = spatialPeriodsList.toArray();
        nSpatialPeriods = _spatialPeriods.length;
        Arrays.sort(_spatialPeriods);

        spatialPeriods = new double[nSpatialPeriods];
        for (int i = 0; i < nSpatialPeriods; i++) {
            spatialPeriods[i] = _spatialPeriods[nSpatialPeriods - i - 1];
        }

        temporalPeriods = temporalPeriodsList.toArray();
        Arrays.sort(temporalPeriods);
        nTemporalPeriods = temporalPeriods.length;

        phases = phaseList.toArray();
        Arrays.sort(phases);
        nPhases = phases.length;

        r.close();
        runIDList = runs.toArray();

        nRuns = new int[paramsMap.size()];
        for (int i = 0; i < runIDList.length; i++) {
            nRuns[runIDList[i]]++;
        }

        index1 = 0;
        index2 = max(nRuns) - 1;

//        System.err.println(index2);

//        System.out.println("runs");
//        for (int i = 0; i < runNumberList.length; i++) {
//            System.out.println(runNumberList[i] + ", ");
//        }
    }


    public DoubleHistogram[] getAverageResponseHistograms() throws
        IOException {

        int[] ttl = neuronFile.getTTLTimes();

        int[] times = neuronFile.getSpikeTimes(currentID);
        final double samplingRate = 20000.0; //Hz

        // check if the run numbers start with 0
        if (min(runIDList) != 0) {
            throw new IllegalStateException("Wrong run ID.");
        }
        final int nDifferentRuns = max(runIDList) + 1;
//        System.out.println("nDifferentRuns " + nDifferentRuns);

        // remove the intermediary ttl's
        IntegerList ttlList = new IntegerList();
        final int ttlsPerRepeat = ttl.length / runIDList.length;
        final int nRepeats = ttl.length / ttlsPerRepeat;
//        System.out.println(" ttl.length " +  ttl.length);

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
//        for (int i = 0; i < 10; i++) {
//            System.out.println(ttl[i]);
//        }

        // create the run histograms
        DoubleHistogram[] runHistograms = new DoubleHistogram[nDifferentRuns];
        for (int i = 0; i < nDifferentRuns; i++) {
            runHistograms[i] = new DoubleHistogram("" + i, 0, averageRunLength,
                timeBinning);
        }

        // construct the spike stream
        SpikeStream stream = new SpikeStream(times, ttl);

        // process the spikes
        int currentRunIndex = -1;
        int t;

        int[] n = new int[runIDList.length];
        Arrays.fill(n, -1);

        while ( (t = stream.getNext()) != Integer.MAX_VALUE) {
            if (t < 0) { // we are switching runs (stimulus conditions)
                currentRunIndex++;

                n[runIDList[currentRunIndex]]++;

                if (currentRunIndex >= runIDList.length) {
                    break;
                }
            } else {
                if (currentRunIndex != -1 &&
                    n[runIDList[currentRunIndex]] >= index1 &&
                    n[runIDList[currentRunIndex]] <= index2) {

                    runHistograms[runIDList[currentRunIndex]].fill(
                        (t - ttl[currentRunIndex]) / samplingRate, 1);
                }
            }
        }

        // normalize the run histograms
        for (int i = 0; i < nDifferentRuns; i++) {
            int nr = Math.min(nRuns[i], index2 - index1 + 1);
            runHistograms[i].scale(1.0 / timeBinning); // convert to spike rate
            runHistograms[i].scale(1.0 / Math.max(nr, 1));
        }

        return runHistograms;
    }


    public String getRunInfo(int i) {
        double[] p = paramsMap.get(new Integer(i));
        return
            "T" + StringUtil.format(p[1], 0) +
            " S" + StringUtil.format(p[0], 0) +
            " P" + StringUtil.format(p[2], 0)
            ;
    }


    /**
     *
     * @param frequency double The frequency in Hz at which to pick up the harmonic
     * @param n int
     * @return double
     */
    private final double getHarmonic(double frequency, int n) {
        int k = (int) Math.round(frequency / omega);
        int dk = 4;
//        int dk = (int) Math.round(1 / omega);

        double max = Double.NEGATIVE_INFINITY;
        for (int i = k - dk; i <= k + dk; i++) {
            if (i >= 0 && real[n][i] > max) {
                max = real[n][i];
            }
        }
        return max;
    }


    public double[] setCurrentNeuron(int id) throws IOException {
        this.currentID = id;
        for (int i = 0; i < real.length; i++) {
            /*Real FT vector [run][frequency], gets reset to magnitude*/
            Arrays.fill(real[i], 0);
            /*Im FT vector [run][freqency] */
            Arrays.fill(img[i], 0);
        }

        // load the reversing grating response histograms
        DoubleHistogram[] h = getAverageResponseHistograms();

        for (int i = 0; i < h.length; i++) {
            // get the histogram
            double[] x = h[i].toArray();

            // subtract the mean
            double xm = mean(x);
            sub(x, xm);

            // the old loop
            //            for (int j = (int) Math.round(1 / timeBinning); j < (int) Math.round(6 / timeBinning); j++) {
            final int n = (int) (stimulusLength_fr / 120.0 / timeBinning);
            for (int j = 0; j < n; j++) {
                real[i][j] = x[j];
            }

            FFT.fft(real[i], img[i], -1);

            for (int j = 0; j < N; j++) {
                real[i][j] = (2.0 / n) *
                             Math.sqrt(real[i][j] * real[i][j] + img[i][j] * img[i][j]);
//                img[i][j] = Math.atan(img[i][j] / real[i][j]);
            }

//            if (i == 56) {
//                PlotPanel  p = PlotUtil.showData(new ScatterPlot(x), style1[0]);
//                p.addToLegend(this.getRunInfo(i));
//
//                PlotUtil.showData(h[i], new HistogramStyle());
//                PlotUtil.showData(new ScatterPlot(real[i]), style1[0]);
//            }
        }

        return null;
    }


    public int getRunID(double spatialPeriod, double temporalPeriod, double phase) {
        for (Integer key : paramsMap.keySet()) {
            double[] p = paramsMap.get(key);

            if (p[0] == spatialPeriod && p[1] == temporalPeriod &&
                p[2] == phase * spatialPeriod) {
                return key.intValue();
            }
        }
        return -1;
    }


    /**
     * If level=1 then F1 and F2 are calculated.
     * If level=2 then F1+F3 and F2+F4 are calculated.
     * and so on
     *
     * @param level int
     * @return double[][][]
     * @throws IOException
     */
    public double[][][] getHarmonic(int level) throws IOException {
        int power = 2;
        double[][][] response = new double[nTemporalPeriods][3][nSpatialPeriods];

        for (int tPeriod = 0; tPeriod < nTemporalPeriods; tPeriod++) {
            for (int sPeriod = 0; sPeriod < nSpatialPeriods; sPeriod++) {
                double maxF1 = Double.NEGATIVE_INFINITY;
                double maxF2 = Double.NEGATIVE_INFINITY;
                double meanF1 = 0;
                double meanF2 = 0;

                for (int phaseIndex = 0; phaseIndex < nPhases; phaseIndex++) {
                    int runID = getRunID(spatialPeriods[sPeriod],
                                         temporalPeriods[tPeriod], phases[phaseIndex]);

                    double f1 = 0, f2 = 0;

                    GaussianMixture m = fit(runID);
                    for (int harmonic = 1; harmonic <= 2 * level; harmonic++) {
                        double f = m.getAmplitude(harmonic - 1);
                        if (harmonic % 2 == 0) { // even harmonics, average
                            f2 += Math.pow(f, power);
                        } else { // odd harmonics, max
                            f1 += Math.pow(f, power);
                        }
                    }

                    f1 = Math.pow(f1, 1.0 / power);
                    f2 = Math.pow(f2, 1.0 / power);

                    meanF1 += f1;
                    meanF2 += f2;
                    if (f1 > maxF1) {
                        maxF1 = f1;
                    }
                    if (f2 > maxF2) {
                        maxF2 = f2;
                    }
                }
                meanF1 /= nPhases;
                meanF2 /= nPhases;

                response[tPeriod][0][sPeriod] = getFrequency(sPeriod);
                response[tPeriod][1][sPeriod] = maxF1;
                response[tPeriod][2][sPeriod] = meanF2;
            }
        }

        return response;
    }


    public double getFrequency(int index) {
        return 1 / spatialPeriods[index] / mmPerPixel;
    }


    public void showPhaseDependence() throws IOException {
        int tPeriod = 0;
        JFrame f = new JFrame();
        f.setLayout(new GridLayout(0, 1));

        for (int sPeriod = 10 - 1; sPeriod >= 0; sPeriod--) {
            PlotPanel p = new PlotPanel();
            ScatterPlot f1 = new ScatterPlot();
            ScatterPlot f2 = new ScatterPlot();
            for (int phaseIndex = 0; phaseIndex < nPhases; phaseIndex++) {
                int runID = getRunID(spatialPeriods[sPeriod],
                                     temporalPeriods[tPeriod], phases[phaseIndex]);
                GaussianMixture m = fit(runID);
                f1.add(phaseIndex * 11.25, m.getAmplitude(0));
                f2.add(phaseIndex * 11.25, m.getAmplitude(1));
            }

            p.addData(f1, style1[0]);
            p.addData(f2, style2[0]);

            p.autoscale();
            p.setLabels("Spatial Phase (deg)", "Amplitude");
            p.setAxisVisible(false);
            f.add(p);
//            PlotUtil.showData("period " + sPeriod, p);

//            ScatterPlotStyle s = new ScatterPlotStyle();
//            s.setConnectingPoints(true);
//            s.setSymbolSize(4);
        }

        f.setBounds(0, 0, 250, 1000);
        f.setVisible(true);
    }


    public void showOnePeriodPSTH() throws IOException {
        boolean[] spatialPeriodState = new boolean[nSpatialPeriods];
         boolean[] temporalPeriodState = new boolean[nTemporalPeriods];
        boolean[] phaseState = new boolean[nPhases];
        Arrays.fill(spatialPeriodState, true);
        Arrays.fill(temporalPeriodState, true);
        Arrays.fill(phaseState, true);
        showOnePeriodPSTH(spatialPeriodState, temporalPeriodState, phaseState, -1);
    }


    public JPanel showOnePeriodPSTH(boolean[] spatialPeriodState, boolean[] temporalPeriodState, boolean[] phaseState,
                                    double max) throws IOException {

        DoubleHistogram[] hh = getOnePeriodHistograms(1);

        if (max == -1) {
            max = Double.MIN_VALUE;
            for (int i = 0; i < hh.length; i++) {
                double v = hh[i].getMaxValue();
                if (v > max) {
                    max = v;
                }
            }
        }

//        JPanel panel = new JPanel(new GridLayout(
//            0, countValues(true, phaseState) * nTemporalPeriods), false);

        JPanel panel = new JPanel(new GridLayout(
            countValues(true, phaseState) * nTemporalPeriods, 0), false);

        panel.setBackground(Color.white);
        HistogramStyle style = new HistogramStyle(
            "", OutlineType.RECTANGULAR, black, 1, true, black, false, black, 1);
        for (int tPeriod = 0; tPeriod < nTemporalPeriods; tPeriod++) {
            if(!temporalPeriodState[tPeriod]) {
                continue;
            }
            for (int sPhase = 0; sPhase < nPhases; sPhase++) {
                if (!phaseState[sPhase]) {
                    continue;
                }

                for (int sPeriod = 0; sPeriod < nSpatialPeriods; sPeriod++) {
                    if (!spatialPeriodState[sPeriod]) {
                        continue;
                    }

                    int n = getRunID(spatialPeriods[sPeriod], temporalPeriods[tPeriod],
                                     phases[sPhase]);

                    PlotPanel p = new PlotPanel();
                    p.addData(hh[n], style);
                    p.autoscale();
                    p.setYRange(0, max);
                    p.setXAxisVisible(false);
                    p.setYAxisVisible(false);
                    p.addToLegend(getRunInfo(n));

                    p.axesBorder.setPadding(1, 1, 1, 1);
                    panel.add(p);
                }
            }
        }
        JFrame f = new JFrame("" + currentID);
        f.setBounds(0, 0, 1500, 800);
//        f.setBounds(0, 0, 100, 100);

        f.add(panel);
        f.setVisible(true);

        return panel;
    }


    public PlotPanel getPSTH(int periodIndex, int phaseIndex, double max) throws
        IOException {
        DoubleHistogram[] hh = getOnePeriodHistograms(1);

        HistogramStyle style = new HistogramStyle(
            "", OutlineType.RECTANGULAR, black, 1 / 8f, true, black, false, black, 1 / 8f);

        int n = getRunID(spatialPeriods[periodIndex], temporalPeriods[0],
                         phases[phaseIndex]);

        PlotPanel p = new PlotPanel("", false, false, false, false);
        p.addData(hh[n], style);
        p.autoscale();
        p.setYRange(0, max);
        p.setXAxisVisible(false);
        p.setYAxisVisible(false);
//        p.axesBorder.setPadding(1, 1, 1, 1);

        return p;
    }


    private GaussianMixture fit(int runID) {
        // create the Gaussian Mixture
        double[] xi = new double[3 * 2];
        double[] ai = new double[3 * 2];
        double tPeriod = paramsMap.get(runID)[1];
        for (int f = 0; f < xi.length; f++) {
            xi[f] = (f + 1) * 120 / tPeriod;
            ai[f] = getHarmonic(xi[f], runID);
        }
        GaussianMixture mixture = new GaussianMixture(xi, ai, 0.1, 0);
        /*
                // fit
                int NN = N / 4;
                double[] x = new double[NN];
                double[] y = new double[NN];
                double[] err = new double[NN];
                for (int f = 0; f < NN; f++) {
                    x[f] = f * omega;
                    y[f] = real[runID][f];
                    err[f] = 1;
                }
                try {
                    Fitter.fit1D(mixture, x, y, err, NN);
                } catch (Exception ex) {
                    ex.printStackTrace();
                }
         */
        return mixture;
    }


    public void showSpectra() {
        ScatterPlotStyle style = new ScatterPlotStyle();
        style.setConnectingPoints(true);
        style.setSymbolType(SymbolType.NONE);
        style.setConnectionLineColor(gray);

        JFrame frame = new JFrame();
        frame.setTitle(currentID + " spectra");
        frame.getContentPane().setLayout(new GridLayout(
            nPhases * nTemporalPeriods, nSpatialPeriods));
        //        frame.getContentPane().setLayout(new GridLayout(
//            nSpatialPeriods, nPhases * nTemporalPeriods));

        double max = Double.NEGATIVE_INFINITY;
        for (int k = 0; k < nTemporalPeriods; k++) {
            for (int j = 0; j < nPhases; j++) {
                for (int i = 0; i < nSpatialPeriods; i++) {
                    int n = getRunID(spatialPeriods[i], temporalPeriods[k], phases[j]);

                    double m = max(real[n], 0, N - 1);
                    if (m > max) {
                        max = m;
                    }

                    ScatterPlot s = new ScatterPlot();
                    for (int ii = 0; ii < N / 2; ii++) {
                        s.add(ii * omega, real[n][ii]);
                    }

                    PlotPanel p = new PlotPanel();
                    p.addData(s, style);
//                    p.addData(fit(n), new FunctionStyle("", black, 1));
                    p.setXAxisVisible(false);
                    p.setYAxisVisible(false);
                    p.addToLegend(getRunInfo(n));
                    p.addToLegend("" + (int) m);
                    p.autoscale();

                    frame.add(p);
                }
            }
        }

        for (int i = 0; i < frame.getContentPane().getComponentCount(); i++) {
            ( (PlotPanel) frame.getContentPane().getComponent(i)).setYRange(0, max);
        }

        frame.setBounds(0, 0, 1000, 800);
        frame.setVisible(true);
    }


    public void showHarmonic(int level) throws IOException {
        double[][][] r = getHarmonic(level);
        PlotPanel p = new PlotPanel();

        for (int tPeriod = 0, cIndex = 0; tPeriod < nTemporalPeriods; tPeriod++) {
            String s1 = "SQUARE 3 " + colors1[cIndex] + ", SOLID 1 " + colors1[cIndex++];
            String s2 = "SQUARE 3 " + colors1[cIndex] + ", SOLID 1 " + colors1[cIndex++];
            p.addData(new ScatterPlot(r[tPeriod][0], r[tPeriod][1], null), s1);
            p.addData(new ScatterPlot(r[tPeriod][0], r[tPeriod][2], null), s2);

//            Spline s = new Spline(nSpatialPeriods);
//            s.reSpline(r[tPeriod][0], r[tPeriod][1]);
//            p.addData(s, new FunctionStyle("F1 spline", Color.black, 1));

//            s = new Spline(nSpatialPeriods);
//            s.reSpline(r[tPeriod][0], r[tPeriod][2]);
//            p.addData(s, new FunctionStyle("F2 spline", Color.red, 1));
        }

        p.setLabels(xLabel, yLabel);
        p.setAxesType(AxisType.LOG10, AxisType.LOG10);
        p.autoscale();
        PlotUtil.showData("", p, 600, 300);
    }


    public static class GaussianMixture
        extends FittableFunction implements FunctionData {

        final int n;
        double[] xi;
        double s;

        public GaussianMixture(double[] xi, double[] ai, double s, double B) {
            this.n = ai.length;
            this.xi = xi;
            this.s = s;

            double[] p = new double[n + 1];
            System.arraycopy(ai, 0, p, 0, n);
            p[n] = B;
            setParameters(p);
        }


        public double getValueAndDerivatives(
            final double[] c, final double[] parameters, final double[] derivatives) {

            double x = c[0];
            double v = 0;
            for (int i = 0; i < n; i++) {
                double d = exp( -0.5 * (x - xi[i]) * (x - xi[i]) / (s * s));
                derivatives[i] = d;
                v += parameters[i] * d;
            }
            derivatives[n] = 1;
            return parameters[n] + v;
        }


        public double getValueAt(double x) {
            double v = 0;
            for (int i = 0; i < n; i++) {
                double d = exp( -0.5 * (x - xi[i]) * (x - xi[i]) / (s * s));
                v += parameters[i] * d;
            }
            return parameters[n] + v;
        }


        public String getDescription() {
            return "";
        }


        public double getAmplitude(int i) {
            return parameters[i];
        }


        public double getBackground() {
            return parameters[n - 1];
        }
    }


    public DoubleHistogram[] getOnePeriodHistograms(int numberOfPeriods) throws
        IOException {

        DoubleHistogram[] h = getAverageResponseHistograms();
        DoubleHistogram[] hh = new DoubleHistogram[h.length];

        for (int i = 0; i < h.length; i++) {
            double temporalPeriod = paramsMap.get(i)[1]; // in frames
            double averagingPeriod_s = numberOfPeriods * (temporalPeriod / 120);
            int binsPerAveragingPeriod = (int) Math.round(averagingPeriod_s / timeBinning);
            int totalPeriods = (int) Math.floor(
                (stimulusLength_fr / 120) / averagingPeriod_s);

            hh[i] = new DoubleHistogram("", 0, averagingPeriod_s, timeBinning);
            for (int binID = 0; binID < totalPeriods * binsPerAveragingPeriod; binID++) {
                hh[i].fillBin(
                    (int) (binID % binsPerAveragingPeriod), h[i].getBin(binID));
            }
            hh[i].scale( (double) numberOfPeriods / totalPeriods);
        }

        return hh;
    }

}


/*
    public static Num get3dbCutoff(DOG1DFourierFunction ff, double wMax) {
        double fMax = ff.getValueAt(wMax);

        // get the 3db cutoff frequency
        double wCutoff = -1;
        for (double wi = wMax; wi < 160; wi += 0.01) {
            double db = 20 * Math.log10(ff.getValueAt(wi) / fMax);
            if (db < -3) {
                wCutoff = wi;
                break;
            }
        }

        return new Num(wCutoff, 0);
    }
 */

/*
    public double getNonlinearityIndex(int sPeriod) throws IOException {
        double[][][] response = getHarmonic(1);

        int tPeriod = 0;
        double f1 = response[tPeriod][1][sPeriod];
        double f2 = response[tPeriod][2][sPeriod];

        return f2 / f1;
    }
 */

/*
    public static double getNonlinearityIndex(int id, ParametersFile pFile) throws
        IOException {

        double[] f1 = pFile.geArrayCell(id, "T1reversingF1");
        double[] f2 = pFile.geArrayCell(id, "T1reversingF2");
        if (f1 == null && f2 == null) {
            return Double.NaN;
        }

        double maxN = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < f1.length; i++) {
            double n = f2[i] / f1[i];
            if (n > maxN) {
                maxN = n;
            }
        }
        return maxN;
    }
 */
