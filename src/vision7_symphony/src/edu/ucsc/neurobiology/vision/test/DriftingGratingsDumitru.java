package edu.ucsc.neurobiology.vision.test;

import java.io.*;

import static java.lang.Math.*;
import java.util.*;

import java.awt.*;
import javax.swing.*;
import javax.swing.border.*;

import edu.ucsc.neurobiology.vision.analysis.*;
import static edu.ucsc.neurobiology.vision.io.IOUtil.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.math.fitting.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.plot.HistogramStyle.*;
import edu.ucsc.neurobiology.vision.util.*;


public class DriftingGratingsDumitru {
    ParametersFile paramsFile;
    private HashMap<Integer, double[]> paramsMap;
    private int[] runNumberList;
    private int[] nRuns;
    private NeuronFile neuronFile;
    private final int N = 2 * 1024;
    private double dt, omega;
    private double[] noise;
    private double[][] real, img;
    private int currentID;
    private double contrast;
    private double stimulusLength_fr;
    public double[] spatialPeriods, temporalPeriods, directions;
    public int nSpatialPeriods, nTemporalPeriods, nDirections;
    private double mmPerPixel = .0058;
    String experimentName, datasetName;
    static String xLabel = "Spatial Frequency (cycles/mm)";
    static String yLabel = "<html>F<sub>1</sub> , F<sub>2</sub> amplitude";


    public DriftingGratingsDumitru(ParametersFile paramsFile) {
        this.paramsFile = paramsFile;
    }


    static Color[] colors = {
                            Color.black, Color.blue, Color.cyan, Color.darkGray,
                            Color.gray, Color.green, Color.lightGray, Color.magenta,
                            Color.orange, Color.pink, Color.red, Color.yellow};

    public DriftingGratingsDumitru(String datasetFolder) throws IOException {
        experimentName = new File(datasetFolder).getParentFile().getName();
        datasetName = new File(datasetFolder).getName();
        neuronFile = new NeuronFile(
            datasetFolder + File.separator + datasetName + ".neurons");
        loadStimulus(
            new File(datasetFolder).getParent() + File.separator + "stimuli" +
            File.separator + datasetName + ".txt");

        noise = new double[nDirections * nSpatialPeriods * nTemporalPeriods];
        real = new double[nDirections * nSpatialPeriods * nTemporalPeriods][N];
        img = new double[nDirections * nSpatialPeriods * nTemporalPeriods][N];
    }


    public double getContrast() {
        return contrast;
    }


    private void loadStimulus(String fileName) throws IOException {
        paramsMap = new HashMap<Integer, double[]>();
        InputStreamReader r = new InputStreamReader(new FileInputStream(fileName));
        StreamTokenizer st = new StreamTokenizer(r);
        IntegerList runs = new IntegerList();
        st.whitespaceChars('(', '(');
        st.whitespaceChars(')', ')');
        st.whitespaceChars(':', ':');
        st.whitespaceChars('#', '#');
        st.eolIsSignificant(true);

        // read the stimulus type
        st.nextToken();
        if (st.ttype != StreamTokenizer.TT_WORD || !st.sval.toUpperCase().equals("TYPE")) {
            throw new IOException("Missing type.");
        }
        st.nextToken();
        if (st.ttype != StreamTokenizer.TT_WORD) {
            throw new IOException("Wrong type.");
        }
//        String type = st.sval;

        // read the rgb
        st.nextToken();
        if (st.ttype != StreamTokenizer.TT_WORD || !st.sval.toUpperCase().equals("RGB")) {
            throw new IOException("Missing rgb.");
        }
        for (int i = 0; i < 3; i++) {
            st.nextToken();
            if (st.ttype != StreamTokenizer.TT_NUMBER) {
                throw new IOException("Wrong type.");
            }
            contrast = st.nval;
        }

        readNumber(st, "X-START");
        readNumber(st, "X-END");
        readNumber(st, "Y-START");
        readNumber(st, "Y-END");

        stimulusLength_fr = readNumber(st, "FRAMES");

        // skip to the end of line
        do {
            st.nextToken();
        } while (st.ttype != StreamTokenizer.TT_EOL);

        DoubleList spatialPeriodsList = new DoubleList();
        DoubleList temporalPeriodsList = new DoubleList();
        DoubleList directionList = new DoubleList();
        // read the run lines
        while (true) {
            double spatialPeriod = readNumber(st, "SPATIAL-PERIOD");
            double temporalPeriod = readNumber(st, "TEMPORAL-PERIOD");
            double direction = readNumber(st, "DIRECTION");

            if (!spatialPeriodsList.contains(spatialPeriod)) {
                spatialPeriodsList.add(spatialPeriod);
            }
            if (!temporalPeriodsList.contains(temporalPeriod)) {
                temporalPeriodsList.add(temporalPeriod);
            }
            if (!directionList.contains(direction)) {
                directionList.add(direction);
            }

            // read the EOL
            st.nextToken();
            if (st.ttype != StreamTokenizer.TT_EOL) {
                throw new IOException("Too many values at line " + st.lineno());
            }

            // figure out what kind of run we have
            int runNumber = -1;
            for (int i = 0; i < paramsMap.size(); i++) {
                double[] params = paramsMap.get(new Integer(i));
                if (spatialPeriod == params[0] && temporalPeriod == params[1] &&
                    direction == params[2]) {

                    runNumber = i;
                    break;
                }
            }

            // yet unknown run, register it
            if (runNumber == -1) {
                runNumber = paramsMap.size();
                double[] params = {spatialPeriod, temporalPeriod, direction};
                paramsMap.put(new Integer(runNumber), params);
            }
            runs.add(runNumber);

            // see whether we reached the end of the file
            st.nextToken();
            if (st.ttype == StreamTokenizer.TT_EOF) {
                break;
            } else {
                st.pushBack();
            }
        }

        spatialPeriods = spatialPeriodsList.toArray();
        Arrays.sort(spatialPeriods);
        nSpatialPeriods = spatialPeriods.length;
        System.out.println("sPeriod: ");
        IOUtil.printArray(spatialPeriods);

        temporalPeriods = temporalPeriodsList.toArray();
        Arrays.sort(temporalPeriods);
        nTemporalPeriods = temporalPeriods.length;
        System.out.println("tPeriod: ");
        IOUtil.printArray(temporalPeriods);

        directions = directionList.toArray();
        Arrays.sort(directions);
        nDirections = directions.length;
        System.out.println("directions: ");
        IOUtil.printArray(directions);

        r.close();
        runNumberList = runs.toArray();

        nRuns = new int[paramsMap.size()];
        for (int i = 0; i < runNumberList.length; i++) {
            nRuns[runNumberList[i]]++;
        }
    }


    public DoubleHistogram[] getAverageResponseHistograms(double dt) throws
        IOException {

        int[] ttl = neuronFile.getTTLTimes();

        int[] times = neuronFile.getSpikeTimes(currentID);
        final double samplingRate = 20000.0; //Hz

        // check if the run numbers start with 0
        if (MathUtil.min(runNumberList) != 0) {
            throw new IllegalStateException("Wrong run ID.");
        }
        final int nDifferentRuns = MathUtil.max(runNumberList) + 1;
//        System.out.println("nDifferentRuns " + nDifferentRuns);

        // remove the intermediary ttl's
        IntegerList ttlList = new IntegerList();
        final int ttlsPerRepeat = ttl.length / runNumberList.length;
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

                if (currentRun >= runNumberList.length) {
                    break;
                }

            } else {
                if (currentRun != -1) {
                    runHistograms[runNumberList[currentRun]].fill(
                        (t - ttl[currentRun]) / samplingRate, 1);
                }
            }
        }

        // normalize the run histograms
        for (int i = 0; i < nDifferentRuns; i++) {
            runHistograms[i].scale(1.0 / Math.max( (nRuns[i] * dt), 1 * dt));
        }

        return runHistograms;
    }


    public String getRunInfo(int i) {
        double[] p = paramsMap.get(new Integer(i));
        return
            "T " + StringUtil.format(120 / p[1], 0) +
            ", S " + StringUtil.format(p[0], 0)
//            + " P" + StringUtil.format(p[2], 0)
            /*+ " noise " + StringUtil.format(noise[i], 1)*/
            ;
    }


    private final double getHarmonic(double frequency, int runID) {
        int k = (int) Math.round(frequency / omega);
//        int dk = (int) Math.round(0.1 * 4 / omega);

//        double sum = 0;
//        for (int i = k - dk; i <= k + dk; i++) {
//            sum += a[n][i];
//        }
//        return sum - (2 * dk + 1) * noise[n];


        double max = Double.NEGATIVE_INFINITY;
        for (int i = k - 4; i <= k + 4; i++) {
            if (real[runID][i] > max) {
                max = real[runID][i];
            }
        }
        return max - noise[runID];
    }


    public void setBinning(double dt) {
        this.dt = dt;
        omega = 1.0 / (N * dt / 1000);
    }


    public double[] setCurrentNeuron(int id) throws IOException {
        this.currentID = id;
        for (int i = 0; i < real.length; i++) {
            /*Real FT vector [run][frequency], gets reset to magnitude*/
            Arrays.fill(real[i], 0);
            /*Im FT vector [run][freqency] */
            Arrays.fill(img[i], 0);
        }

        double T = dt / 1000;

        // load the reversing grating response histograms
        DoubleHistogram[] h = getAverageResponseHistograms(T);

//        double t = 0;

        for (int i = 0; i < h.length; i++) {
            // get the histogram
            double[] x = h[i].toArray();

            // subtract the mean
            double xm = MathUtil.mean(x);
            MathUtil.sub(x, xm);
            //I commented out the "use only six frames", not clear why it was here.  MG
            //Without it, changing binning to .5 causes problems.

            for (int j = (int) Math.round(1 / T); j < (int) Math.round(6 / T); j++) {
//            for (int j = 0; j < x.length; j++) {
                real[i][j] = x[j];
                img[i][j] = 0;
            }

//            long t1 = System.currentTimeMillis();
//            for (int k = 0; k < 1; k++) {
            FFT.fft(real[i], img[i], -1);
//            }
//            long t2 = System.currentTimeMillis();
//            t += (t2 - t1) / 1000.;

            for (int j = 0; j < N / 2; j++) {
                double amp = Math.sqrt(real[i][j] * real[i][j] + img[i][j] * img[i][j]) *
                             T;
                double phase = Math.atan(img[i][j] / real[i][j]);

                real[i][j] = amp;
                img[i][j] = phase;
//                real[i][j] = phase;
//                img[i][j] = amp;
            }
            noise[i] = MathUtil.mean(real[i], N / 4, N / 2);
        }

//        System.err.println(t);

        return null;
    }


    public int getRunID(double spatialPeriod, double temporalPeriod, double direction) {
        for (Integer key : paramsMap.keySet()) {
            double[] p = paramsMap.get(key);
            if (p[0] == spatialPeriod && p[1] == temporalPeriod && p[2] == direction) {
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
        double[][][] response = new double[nTemporalPeriods][3][nSpatialPeriods];

        for (int tPeriod = 0; tPeriod < nTemporalPeriods; tPeriod++) {
            for (int sPeriod = 0; sPeriod < nSpatialPeriods; sPeriod++) {
                double meanF1 = 0;
                double meanF2 = 0;

                for (int dirIndex = 0; dirIndex < nDirections; dirIndex++) {
                    int runID = getRunID(spatialPeriods[sPeriod],
                                         temporalPeriods[tPeriod], directions[dirIndex]);

                    double f1 = 0, f2 = 0;
                    ReversingGratings.GaussianMixture m = fit(runID);
                    for (int harmonic = 1; harmonic <= 2 * level; harmonic++) {
                        double f = m.getAmplitude(harmonic - 1);
                        if (harmonic % 2 == 0) { // even harmonics, average
                            f2 += Math.pow(f, 2);
                        } else { // odd harmonics, max
                            f1 += Math.pow(f, 2);
                        }
                    }
                    f1 = Math.pow(f1, 0.5);
                    f2 = Math.pow(f2, 0.5);

                    meanF1 += f1;
                    meanF2 += f2;
                }
                meanF1 /= nDirections;
                meanF2 /= nDirections;

                response[tPeriod][0][sPeriod] = 1 / spatialPeriods[sPeriod] / mmPerPixel;
                response[tPeriod][1][sPeriod] = meanF1;
                response[tPeriod][2][sPeriod] = meanF2;
            }
        }

        return response;
    }


    public static PlotPanel getHarmonicsPanel(int[] id, ParametersFile paramsFile) {
        double[] freqs = paramsFile.getArrayCell(id[0], "reversingFrequencies");
        PlotPanel p = new PlotPanel("F1F2");

        // make the styles
        ScatterPlotStyle[] style1 = new ScatterPlotStyle[3];
        ScatterPlotStyle[] style2 = new ScatterPlotStyle[3];
        for (int tPeriod = 0, colorIndex = 0; tPeriod < 3; tPeriod++) {
            style1[tPeriod] = new ScatterPlotStyle("Fodd " + tPeriod,
                SymbolType.FILLED_SQUARE, 4, colors[colorIndex], true,
                colors[colorIndex++], 1);
            style2[tPeriod] = new ScatterPlotStyle("Feven " + tPeriod,
                SymbolType.FILLED_SQUARE, 4, colors[colorIndex], true, colors[colorIndex++],
                1);
        }

        // calculate overall scaling factor
        double overallsum = 0;
        int nn = 0;
        for (int i = 0; i < id.length; i++) {
            for (int tPeriod = 0; tPeriod < 3; tPeriod++) {
                double[] _f1 = paramsFile.getArrayCell(id[i],
                    "T" + (tPeriod + 1) + "reversingF1");
                double[] _f2 = paramsFile.getArrayCell(id[i],
                    "T" + (tPeriod + 1) + "reversingF2");

                if (_f1 != null && _f2 != null) {
                    overallsum += MathUtil.sum(_f1) + MathUtil.sum(_f2);
                    nn++;
                }
            }
        }
        overallsum /= nn;

        // make the graphs
        for (int i = 0; i < id.length; i++) {
            double sum = 0.0;

            for (int tPeriod = 0; tPeriod < 3; tPeriod++) {
                double[] _f1 = paramsFile.getArrayCell(id[
                    i], "T" + (tPeriod + 1) + "reversingF1");
                double[] _f2 = paramsFile.getArrayCell(id[
                    i], "T" + (tPeriod + 1) + "reversingF2");

                if (_f1 != null && _f2 != null) {
                    sum += MathUtil.sum(_f1);
                    sum += MathUtil.sum(_f2);

                    MathUtil.multiply(_f1, overallsum / sum);
                    MathUtil.multiply(_f2, overallsum / sum);

                    p.addData(new ScatterPlot(freqs, _f1, null), style1[tPeriod]);
                    p.addData(new ScatterPlot(freqs, _f2, null), style2[tPeriod]);
                }
            }
        }

        p.setRange(0, 3, 0, 150);
        p.setLabels(xLabel, yLabel);
        p.setAxesType(AxisType.LOG10, AxisType.LOG10);
        p.autoscale();
        return p;
    }


    public void showOnePeriodPSTH() throws IOException {
        DoubleHistogram[] hh = getResponse(dt / 1000, 1);
        double max = Double.MIN_VALUE;
        for (int i = 0; i < hh.length; i++) {
            double v = hh[i].getMaxValue();
            if (v > max) {
                max = v;
            }
        }

        JPanel panel = new JPanel(
            new GridLayout(0, nDirections * nTemporalPeriods /* / 2*/));
        HistogramStyle style = new HistogramStyle(
            "", OutlineType.RECTANGULAR, Color.black, 1, true, Color.black, false,
            Color.black, 1);
        for (int i = 0; i < nSpatialPeriods; i++) {
//            if (i != 3 && i != 9) {
//                continue;
//            }
            for (int k = 0; k < nTemporalPeriods; k++) {
                for (int j = 0; j < nDirections; j++) {
//                    if (j % 2 != 0) {
//                        continue;
//                    }

                    int n = getRunID(spatialPeriods[i], temporalPeriods[k], directions[j]);

                    PlotPanel p = new PlotPanel();
                    p.addData(hh[n], style);
                    p.autoscale();
                    p.setYRange(0, max);
                    p.setXAxisVisible(false);
                    p.setYAxisVisible(false);
//                    p.addToLegend(getRunInfo(n));

                    Border b = BorderFactory.createEmptyBorder(2, 2, 2, 2);
                    p.setBorder(BorderFactory.createCompoundBorder(b, p.getBorder()));
                    panel.add(p);
                }
            }
        }
        JFrame f = new JFrame("" + currentID);
        f.setBounds(0, 0, 1500, 800);
        f.add(panel);
        f.setVisible(true);

//        saveAsEPS(panel, "" + currentID, 6, 2);
    }


    private ReversingGratings.GaussianMixture fit(int runID) {
        double[] xi = new double[3 * 2];
        double[] ai = new double[3 * 2];
        for (int f = 0; f < xi.length; f++) {
            double tPer = 120 / paramsMap.get(new Integer(runID))[1];
            xi[f] = tPer * (f + 1);
            ai[f] = getHarmonic(tPer * (f + 1), runID);
        }
        int NN = N / 2;
        double[] x = new double[NN];
        double[] y = new double[NN];
        double[] err = new double[NN];
        for (int f = 0; f < NN; f++) {
            x[f] = f * omega;
            y[f] = real[runID][f];
            err[f] = 1;
        }

        ReversingGratings.GaussianMixture mixture = new ReversingGratings.GaussianMixture(
            xi, ai, 0.1, noise[runID]);
        try {
//            Fitter.fit1D(mixture, x, y, err, NN);
        } catch (Exception ex) {
            System.err.println("ex");
//            ex.printStackTrace();
        }

        return mixture;
    }


    public void showSpectra() {
        ScatterPlotStyle style = new ScatterPlotStyle();
        style.setConnectingPoints(true);
        style.setSymbolType(SymbolType.NONE);
        style.setConnectionLineColor(Color.gray);

        JFrame frame = new JFrame();
        frame.setTitle(currentID + " spectra");
        frame.getContentPane().setLayout(new GridLayout(nSpatialPeriods,
            nDirections * nTemporalPeriods));

        double max = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < nSpatialPeriods; i++) {
            for (int k = 0; k < nTemporalPeriods; k++) {
                for (int j = 0; j < nDirections; j++) {
                    int n = getRunID(spatialPeriods[i], temporalPeriods[k], directions[j]);

                    double m = MathUtil.max(real[n], 0, N / 2 - 1);
                    if (m > max) {
                        max = m;
                    }

                    ScatterPlot s = new ScatterPlot();
                    for (int ii = 0; ii < N / 2; ii++) {
                        s.add(ii * omega, real[n][ii]);
                    }

                    PlotPanel p = new PlotPanel();
                    p.addData(s, style);
                    p.addData(fit(n), new FunctionStyle("", Color.black, 1));

                    p.setXAxisVisible(true);
                    p.setYAxisVisible(false);
                    p.addToLegend(getRunInfo(n));
//                    p.addToLegend("" + (int) m);
                    p.autoscale();
                    p.setXRange(0, 5 * 120 / temporalPeriods[k]);

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


    public DoubleHistogram[] getResponse(double dt_sec, int numberOfPeriods) throws
        IOException {
        DoubleHistogram[] h = getAverageResponseHistograms(dt_sec);
        DoubleHistogram[] hh = new DoubleHistogram[h.length];

        for (int i = 0; i < h.length; i++) {
            double temporalPeriod = paramsMap.get(new Integer(i))[1];

            double averagingPeriod_s = numberOfPeriods * temporalPeriod * 8.34 / 1000;
            int binsPerAveragingPeriod = (int) Math.round(averagingPeriod_s / dt_sec);
            int totalPeriods = (int) Math.ceil(
                stimulusLength_fr * 8.34 / 1000 / averagingPeriod_s);

            hh[i] = new DoubleHistogram("", 0, averagingPeriod_s, dt_sec);
            for (int binID = 0; binID < totalPeriods * temporalPeriod; binID++) {
                hh[i].fillBin( (int) (binID % binsPerAveragingPeriod),
                              h[i].getBin(binID));
            }
            hh[i].scale( (double) numberOfPeriods / totalPeriods);
        }

        return hh;
    }


    public void showHarmonics(int[] id, String name) throws IOException {
        int nPeriods = spatialPeriods.length;

        ScatterPlotStyle st1 = new ScatterPlotStyle(
            SymbolType.FILLED_SQUARE, 4, Color.black, true, Color.black,
            1);
        ScatterPlotStyle st2 = new ScatterPlotStyle(
            SymbolType.FILLED_SQUARE, 4, Color.red, true, Color.red, 1);
        PlotPanel p = new PlotPanel();
        double[] _f1 = new double[nPeriods];
        double[] _f2 = new double[nPeriods];

        for (int neuron = 0; neuron < id.length; neuron++) {
            setCurrentNeuron(id[neuron]);
            for (int periodIndex = 0; periodIndex < nPeriods; periodIndex++) {
                double maxF1 = Double.NEGATIVE_INFINITY;
                double meanF2 = 0;

                for (int dirIndex = 0; dirIndex < nDirections; dirIndex++) {
                    int runID = getRunID(spatialPeriods[periodIndex], temporalPeriods[0],
                                         directions[dirIndex]);
                    double f = 120 / paramsMap.get(new Integer(runID))[1];
                    double f1 = getHarmonic(f * 1, runID);
                    if (f1 > maxF1) {
                        maxF1 = f1;
                    }
                    double f2 = getHarmonic(f * 2, runID);
                    meanF2 += f2;
                }
                meanF2 /= nDirections;

                _f1[periodIndex] = maxF1;
                _f2[periodIndex] = meanF2;
            }

            double norm = MathUtil.sum(_f1);
            MathUtil.divide(_f1, norm);
            MathUtil.divide(_f2, norm);

//            double maxNonlinearity = Double.NEGATIVE_INFINITY;
            ScatterPlot tunungF1 = new ScatterPlot();
            ScatterPlot tunungF2 = new ScatterPlot();
            for (int periodIndex = 0; periodIndex < nPeriods; periodIndex++) {
                double freq = Math.log( (1 /
                                         (spatialPeriods[periodIndex] * 6.0 / 210.0)));
                tunungF1.add(freq, _f1[periodIndex]);
                tunungF2.add(freq, _f2[periodIndex]);

//                double n = _f2[periodIndex] / _f1[periodIndex];
//                if (n > maxNonlinearity) {
//                    maxNonlinearity = n;
//                }
            }

            p.addData(tunungF1, st1);
            p.addData(tunungF2, st2);

//            System.out.println(id[neuron] + " : " + maxNonlinearity);
        }

        p.setRange(0, 3, 0, 150);
        p.setLabels(xLabel, yLabel);
        p.autoscale();
        PlotUtil.showData(name, p);
    }


    public void showHarmonics(String classID) throws IOException {
        int[] id = paramsFile.getNeuronsInClass(classID);
        showHarmonics(id, classID);
    }


    public void showHarmonic(int level) throws IOException {
        double[][][] response = getHarmonic(level);

        PlotPanel p = new PlotPanel();
//        int cIndex = 0;
        for (int tPeriod = 0; tPeriod < nTemporalPeriods; tPeriod++) {
            double area = 0;
            for (int sPeriod = 7; sPeriod < nSpatialPeriods; sPeriod++) {
                area += response[tPeriod][1][sPeriod];
            }

            ScatterPlot tunungF1 = new ScatterPlot();
            ScatterPlot tunungF2 = new ScatterPlot();
            for (int sPeriod = 0; sPeriod < nSpatialPeriods; sPeriod++) {
                tunungF1.add(response[tPeriod][0][sPeriod],
                             response[tPeriod][1][sPeriod] / area);
                tunungF2.add(response[tPeriod][0][sPeriod],
                             response[tPeriod][2][sPeriod] / area);
            }
//            p.addData(tunungF1, new ScatterPlotStyle("" + 120 / temporalPeriods[tPeriod],
//                SymbolType.FILLED_SQUARE, 4, colors[cIndex], true, colors[cIndex], 1));
//            cIndex++;
            p.addData(tunungF1, new ScatterPlotStyle("" + 120 / temporalPeriods[tPeriod],
                SymbolType.FILLED_SQUARE, 4, Color.black, true, Color.black, 1));
            p.addData(tunungF2, new ScatterPlotStyle("" + 120 / temporalPeriods[tPeriod],
                SymbolType.FILLED_SQUARE, 4, Color.red, true, Color.red, 1));
        }
        p.setLabels("Spatial Period (syc/mm)", "Response");
        p.setAxesType(AxisType.LOG10, AxisType.LOG10);
        p.autoscale();
        PlotUtil.showData("", p, 600, 400);
        /*
                p = new PlotPanel();
                cIndex = 0;
                for (int sPeriod = 0; sPeriod < nSpatialPeriods; sPeriod++) {
                    ScatterPlot tunungF1 = new ScatterPlot();
                    ScatterPlot tunungF2 = new ScatterPlot();
                    for (int tPeriod = 0; tPeriod < nTemporalPeriods; tPeriod++) {
                        double f = 120 / temporalPeriods[tPeriod];
                        tunungF1.add(f, response[tPeriod][1][sPeriod]);
                        tunungF2.add(f, response[tPeriod][2][sPeriod]);
                    }
                    p.addData(tunungF1, new ScatterPlotStyle("" + StringUtil.format(1 / spatialPeriods[sPeriod] / mmPerPixel, 2),
         SymbolType.FILLED_SQUARE, 4, colors[cIndex], true, colors[cIndex], 1));
                    cIndex++;
                }
                p.setLabels("Temporal Period (Hz)", "Response");
                p.setAxesType(AxisType.LOG10, AxisType.LOG10);
                p.autoscale();
                PlotUtil.showData("", p, 600, 400);
         */
    }


    public static Num get3dbCutoff(DOG1DFourierFunction ff, double wMax) {
        double fMax = ff.getValueAt(wMax);

        // get the 3db cutoff frequency
        double wCutoff = -1;
        for (double wi = wMax; wi < 160; wi += 0.01) {
            double db = 20 * log10(ff.getValueAt(wi) / fMax);
            if (db < -3) {
                wCutoff = wi;
                break;
            }
        }

        return new Num(wCutoff, 0);
    }


    public static DOG1DFourierFunction fit(
        double[] freqs, double[] f, double[] fErr,
        PlotPanel p1, String name, HashMap<String, Num> params) {

        DOG1DFourierFunction ff = null;
        for (double s = 0.01; s < 1; s += 0.005) {
            DOG1DFourierFunction fun;
            double max = MathUtil.max(f);
//            double s = 0.100;
//            if (name.equals("F1")) {
//                fun = new DOG1DFourierFunction(max, s, 0);
//                fun.setParameterState(2, false);
//            } else {
            fun = new DOG1DFourierFunction(max, -max / 5, s, s * 10, f[0]);
//            }

            try {
                Fitter.fit1D(fun, freqs, f, fErr, freqs.length);
                if (ff != null) {
                    if (fun.getChiSquared() <
                        ff.getChiSquared() /*&&  fun.getS2() > 1.5 * fun.getS1()*/) {
                        ff = fun;
                    }
                } else {
                    ff = fun;
                }
            } catch (Exception ex) {
            }
        }
        /*
                // get the frequency of maximum
                Num wMax;
                if (ff.getA2() < 0) {
//            double v = -a2 * s2 * s2 / (a1 * s1 * s1);
//            w = sqrt(2 * log(v) / (s2 * s2 - s1 * s1));
                    wMax = ff.a2().mul(ff.s2().sqr()).inv();
                    wMax = wMax.div(ff.a1()).div(ff.s1().sqr());

                    Num h = Num.zero.add(ff.s2().sqr()).sub(ff.s1().sqr());
                    wMax = wMax.ln().mul(Num.two).div(h).sqrt();
                } else {
                    wMax = Num.zero;
                }

                // get the amaplitude at max
                Num fMax = ff.getValueAt(wMax);

                // get the 3db cutoff frequency and error
                DOG1DFourierFunction f1 = new DOG1DFourierFunction(0, 0, 0, 0, 0);
                DOG1DFourierFunction f2 = new DOG1DFourierFunction(0, 0, 0, 0, 0);
                double error = Double.NEGATIVE_INFINITY;
                int param = -1;
                for (int i = 0; i < ff.getParametersCount(); i++) {
                    f1.setParameters(ff.getParameters());
                    f2.setParameters(ff.getParameters());
                    f1.setParameter(i, ff.getParameter(i) + ff.getParameterError(i));
                    f2.setParameter(i, ff.getParameter(i) - ff.getParameterError(i));

                    Num cutoff1 = get3dbCutoff(f1, wMax.x);
                    Num cutoff2 = get3dbCutoff(f2, wMax.x);
                    if (Math.abs(cutoff1.x - cutoff2.x) > error) {
                        error = Math.abs(cutoff1.x - cutoff2.x);
                        param = i;
                    }
                }

                Num wCutoff = get3dbCutoff(ff, wMax.x);
                wCutoff = new Num(wCutoff.x, error);
         */
        p1.addToLegend(
            "a1=" + ff.a1().toString(1) + ", " +
            "a2=" + ff.a2().toString(1) + ", " +
            "s1=" + ff.s1().toString(3) + ", " +
            "s2=" + ff.s2().toString(3) + ", " +
            "b=" + ff.b().toString(3) + ", " +
            "\u03C7\u00B2=" + StringUtil.format(ff.getChiSquared(), 1) + ", "
            );
        /*
                p1.addToLegend("wc = " + wCutoff.toString(1) + " @ " +
                               ff.getParameterNames()[param]);

                params.put(name + "wMax", wMax);
                params.put(name + "wCutoff", wCutoff);
                params.put(name + "fMax", fMax);
         */
        return ff;
    }


    public static PlotPanel getAverageHarmonicsPanel(
        int[] id, ParametersFile pFile, boolean db, ObjectOutputStream out) throws
        IOException {

        double[] _freqs = pFile.getArrayCell(id[0], "reversingFrequencies");
        final int nFreq = _freqs.length;
//        double[] freqs = new double[nFreq];
        for (int i = 0; i < nFreq; i++) {
            _freqs[i] *= 2 * Math.PI;
        }

        MeanVarianceCalculator[] mvcF1 = new MeanVarianceCalculator[nFreq];
        MeanVarianceCalculator[] mvcF2 = new MeanVarianceCalculator[nFreq];
        for (int periodIndex = 0; periodIndex < nFreq; periodIndex++) {
            mvcF1[periodIndex] = new MeanVarianceCalculator();
            mvcF2[periodIndex] = new MeanVarianceCalculator();
        }

        // calculate overall scaling factor
        double overallsum = 0;
        int nn = 0;
        for (int index = 0; index < id.length; index++) {
            double[] _f1 = pFile.getArrayCell(id[index], "T1reversingF1");
            double[] _f2 = pFile.getArrayCell(id[index], "T1reversingF2");
            if (_f1 != null) {
                overallsum += MathUtil.sum(_f1) + MathUtil.sum(_f2);
                nn++;
            }
        }
        overallsum /= nn;

        for (int index = 0; index < id.length; index++) {
            double[] _f1 = pFile.getArrayCell(id[index], "T1reversingF1");
            double[] _f2 = pFile.getArrayCell(id[index], "T1reversingF2");

            if (_f1 != null) {
                double sum = MathUtil.sum(_f1) + MathUtil.sum(_f2);
                if (sum > 0.0) {
                    MathUtil.multiply(_f1, overallsum / sum);
                    MathUtil.multiply(_f2, overallsum / sum);
                }

                for (int periodIndex = 0; periodIndex < nFreq; periodIndex++) {
                    mvcF1[periodIndex].add(_f1[periodIndex]);
                    mvcF2[periodIndex].add(_f2[periodIndex]);
                }
            }
        }

        double[] f1 = new double[nFreq];
        double[] f2 = new double[nFreq];
        double[] f1Err = new double[nFreq];
        double[] f2Err = new double[nFreq];
        for (int sPeriod = 0; sPeriod < nFreq; sPeriod++) {
            f1[sPeriod] = mvcF1[sPeriod].getMean();
            f1Err[sPeriod] = mvcF1[sPeriod].getMeanVariance();
            f2[sPeriod] = mvcF2[sPeriod].getMean();
            f2Err[sPeriod] = mvcF1[sPeriod].getMeanVariance();
        }

        PlotPanel p1 = new PlotPanel("F1F2db");

        // fit the F1 and F2 curves
        HashMap<String, Num> params = new HashMap<String, Num>();
        int nnn = 10;
        double[] _f1 = new double[nnn];
        double[] _f1Err = new double[nnn];
        double[] _fr = new double[nnn];
        for (int i = 0; i < _f1.length; i++) {
            _f1[i] = f1[i + nFreq - nnn];
            _f1Err[i] = f1Err[i + nFreq - nnn];
            _fr[i] = _freqs[i + nFreq - nnn];
        }

        DOG1DFourierFunction bestF1 = fit(_fr, _f1, _f1Err, p1, "F1", params);
        DOG1DFourierFunction bestF2 = fit(_freqs, f2, f2Err, p1, "F2", params);
        if (out != null) {
            out.writeObject(bestF1);
            out.writeObject(bestF2);
        }
//            Num f1Max = params.get("F1fMax");
//            Num f2Max = params.get("F2fMax");
//            p1.addToLegend("F2max/F1max = " + f2Max.div(f1Max).toString(2));

        /*
                // claculate "c"
                Num c1 = params.get("F1wCutoff");
                Num c2 = params.get("F2wCutoff");
                Num c = c2.div(c1);
                p1.addToLegend("C2/C1 = " + c.toString(1));
                if (out != null) {
                    out.writeObject(c);
                }

                // calculate "n"
                Num maxN = new Num(Double.NEGATIVE_INFINITY, 0);
                int sMax = -1;
                for (int sPeriod = 0; sPeriod < nFreq; sPeriod++) {
         Num n = new Num(f2[sPeriod], f2Err[sPeriod]).div(f1[sPeriod], f1Err[sPeriod]);
                    if (n.x > maxN.x) {
                        maxN = n;
                        sMax = sPeriod;
                    }
                }
                p1.addToLegend("n = " + maxN.toString(1) + " @ " + freqs[sMax]);
                if (out != null) {
                    out.writeObject(maxN);
                }
         */
        // convert the dB
        if (db) {
            double _f1Max = MathUtil.max(f1);
            double _f2Max = MathUtil.max(f2);
            for (int periodIndex = 0; periodIndex < nFreq; periodIndex++) {
                f1[periodIndex] = 20 * Math.log10(f1[periodIndex] / _f1Max);
                f2[periodIndex] = 20 * Math.log10(f2[periodIndex] / _f2Max);
            }
            p1.setLabels(xLabel, yLabel + " dB");
        } else {
            p1.setLabels(xLabel, yLabel);
        }

//        bestF1.setParameter(1, bestF1.getParameter(1) * (2 * Math.PI));
//        try {
//            bestF1.setParameter(4, bestF1.getParameter(4) * (2 * Math.PI));
//        } catch (Exception e) {}
//
//        bestF2.setParameter(1, bestF2.getParameter(1) * (2 * Math.PI));
//        try {
//            bestF2.setParameter(4, bestF2.getParameter(4) * (2 * Math.PI));
//        } catch (Exception e) {}

        p1.addData(new ScatterPlot(_freqs, f1, f1Err), new ScatterPlotStyle(
            SymbolType.NONE, 4, Color.black, true, Color.black, 1));
        p1.addData(new ScatterPlot(_freqs, f2, f2Err), new ScatterPlotStyle(
            SymbolType.NONE, 4, Color.red, true, Color.red, 1));
        p1.addData(bestF1, new FunctionStyle("Fodd", Color.black, 1));
        p1.addData(bestF2, new FunctionStyle("Feven", Color.red, 1));
        p1.setAxesType(AxisType.LOG10, AxisType.LOG10);
        p1.setLegendLocation(PlotPanel.Location.DOWNLEFT);
        p1.autoscale();
//        p1.setYRange(4, 200);
        p1.setLegendVisible(true);

        p1.addData(new DOG1DFourierFunction(
            400, 400 * -0.84,
            0.160,
            0.160, 0
                   ), new FunctionStyle("RF", Color.gray, 5));

        return p1;
    }


    public static double[][] getAverageF1F2(int[] id, ParametersFile pFile) throws
        IOException {
        final int nFreq = pFile.getArrayCell(id[0], "reversingFrequencies").length;

        MeanVarianceCalculator[] mvcF1 = new MeanVarianceCalculator[nFreq];
        MeanVarianceCalculator[] mvcF2 = new MeanVarianceCalculator[nFreq];
        for (int periodIndex = 0; periodIndex < nFreq; periodIndex++) {
            mvcF1[periodIndex] = new MeanVarianceCalculator();
            mvcF2[periodIndex] = new MeanVarianceCalculator();
        }

        // calculate overall scaling factor
        double overallsum = 0;
        int nn = 0;
        for (int index = 0; index < id.length; index++) {
            double[] _f1 = pFile.getArrayCell(id[index], "T1reversingF1");
            double[] _f2 = pFile.getArrayCell(id[index], "T1reversingF2");
            if (_f1 != null) {
                overallsum += MathUtil.sum(_f1) + MathUtil.sum(_f2);
                nn++;
            }
        }
        overallsum /= nn;

        for (int index = 0; index < id.length; index++) {
            double[] _f1 = pFile.getArrayCell(id[index], "T1reversingF1");
            double[] _f2 = pFile.getArrayCell(id[index], "T1reversingF2");

            if (_f1 != null) {
                double sum = MathUtil.sum(_f1) + MathUtil.sum(_f2);
                if (sum > 0.0) {
                    MathUtil.multiply(_f1, overallsum / sum);
                    MathUtil.multiply(_f2, overallsum / sum);
                }

                for (int periodIndex = 0; periodIndex < nFreq; periodIndex++) {
                    mvcF1[periodIndex].add(_f1[periodIndex]);
                    mvcF2[periodIndex].add(_f2[periodIndex]);
                }
            }
        }

        double[] f1 = new double[nFreq];
        double[] f2 = new double[nFreq];
        double[] f1Err = new double[nFreq];
        double[] f2Err = new double[nFreq];
        for (int sPeriod = 0; sPeriod < nFreq; sPeriod++) {
            f1[sPeriod] = mvcF1[sPeriod].getMean();
            f1Err[sPeriod] = mvcF1[sPeriod].getMeanVariance();
            f2[sPeriod] = mvcF2[sPeriod].getMean();
            f2Err[sPeriod] = mvcF1[sPeriod].getMeanVariance();
        }

        return new double[][] {f1, f2};
    }

}
