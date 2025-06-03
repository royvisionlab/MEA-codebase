package edu.ucsc.neurobiology.vision.analysis;

import java.io.*;
import java.util.*;

import static java.awt.Color.*;
import java.awt.*;
import javax.swing.*;

import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.math.*;
import static edu.ucsc.neurobiology.vision.plot.SymbolType.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class MovingBar {
    final double samplingRate = 20000.0;
    final String datasetName;
    final String baseName;
    final NeuronFile nf;
    final int nDifferentRuns;
    final double[] runLength;
    final int totalRuns;
    final int[] nRuns;

    int[] runNumberList;
    HashMap<Integer, double[]> paramsMap;
    int[] ttl;
    static int nb = 4;


    public MovingBar(String baseName, String datasetName) throws IOException {
        this.baseName = baseName;
        this.datasetName = datasetName;

        paramsMap = new HashMap();
        loadMovingFilterStimuli();

        // rearrange the ID's
        HashMap m = new HashMap();
        double minColor;
        Integer minID;
        int newID = 0;
        while (!paramsMap.isEmpty()) {
            minColor = Double.POSITIVE_INFINITY;
            minID = null;
            for (Integer id : paramsMap.keySet()) {
                double[] p = paramsMap.get(id);
                if (p[0] < minColor) {
                    minColor = p[0];
                    minID = id;
                }
            }
            double[] p = paramsMap.remove(minID);
            m.put(newID, p);
            newID++;
        }
        paramsMap = m;
        loadMovingFilterStimuli();
        nDifferentRuns = paramsMap.size();

        nf = new NeuronFile(baseName + File.separator + datasetName +
                            File.separator + datasetName + ".neurons");

        // initialize all structures
        ttl = nf.getTTLTimes();

        // check if the run numbers start with 0
        if (MathUtil.min(runNumberList) != 0) {
            throw new IllegalStateException("Wrong run ID.");
        }
        final int nDifferentRuns = MathUtil.max(runNumberList) + 1;

        // remove the odd ttl's
        IntegerList ttlList = new IntegerList();
        for (int i = 0; i < ttl.length; i++) {
            if (i % 2 == 0) {
                ttlList.add(ttl[i]);
            }
        }
        ttl = ttlList.toArray();
        totalRuns = ttl.length - 1; // FIXME maybe remove -1

        // calculate the average run length
        runLength = new double[nDifferentRuns];
        nRuns = new int[nDifferentRuns];
        for (int i = 0; i < ttl.length - 1; i++) {
            int d = ttl[i + 1] - ttl[i];
            runLength[runNumberList[i]] += d;
            nRuns[runNumberList[i]]++;
        }
        for (int i = 0; i < runLength.length; i++) {
            runLength[i] /= nRuns[i];
//            System.err.println(runLength[i]);
        }
    }


    private void loadMovingFilterStimuli() throws IOException {
        InputStreamReader r = new InputStreamReader(new FileInputStream(
            baseName + File.separator + "stimuli" + File.separator + datasetName + ".txt"));
        StreamTokenizer st = new StreamTokenizer(r);
        IntegerList _runs = new IntegerList();
        st.whitespaceChars('(', '(');
        st.whitespaceChars(')', ')');
        st.whitespaceChars(':', ':');
        st.whitespaceChars('#', '#');
        st.eolIsSignificant(true);

        // read the stimulus type
        st.nextToken();
        if (st.ttype != st.TT_WORD || !st.sval.equals("TYPE")) {
            throw new IOException("Missing TYPE.");
        }
        st.nextToken();
        if (st.ttype != st.TT_WORD) {
            throw new IOException("Wrong TYPE.");
        }
        String type = st.sval;
//            System.out.println("TYPE : " + type);

        // skip to the end of line
        do {
            st.nextToken();
        } while (st.ttype != st.TT_EOL);

        // read the run lines
        while (true) {
            // read the color of the bar
            st.nextToken();
            if (st.ttype != st.TT_WORD || !st.sval.equals("RGB")) {
                throw new IOException("Missing RGB at line " + st.lineno());
            }
            double[] rgb = new double[3];
            for (int i = 0; i < 3; i++) {
                st.nextToken();
                if (st.ttype != st.TT_NUMBER) {
                    throw new IOException("Bad RGB at line " + st.lineno());
                }
                rgb[i] = st.nval;
            }

            // read the direction of motion of the bar
            st.nextToken();
            if (st.ttype != st.TT_WORD) {
                if (!st.sval.equals("X-DELTA") || !st.sval.equals("Y-DELTA")) {
                    throw new IOException("Missing X-DELTA/Y-DELTA at line " + st.lineno());
                }
            }
            st.nextToken();
            if (st.ttype != st.TT_NUMBER) {
                throw new IOException("Missing X-DELTA at line " + st.lineno());
            }
            double speed = st.nval;

            // read the EOL
            st.nextToken();
            if (st.ttype != st.TT_EOL) {
                throw new IOException("Too many values at line " + st.lineno());
            }

            // figure out what kind of run we have
            int runNumber = -1;
            for (int i = 0; i < paramsMap.size(); i++) {
                double[] params = paramsMap.get(i);
                if (rgb[0] == params[0] &&
                    rgb[1] == params[1] &&
                    rgb[2] == params[2] &&
                    speed == params[3]) {

                    runNumber = i;
                    break;
                }
            }
            if (runNumber == -1) {
                runNumber = paramsMap.size();
                double[] params = {rgb[0], rgb[1], rgb[2], speed};
                paramsMap.put(runNumber, params);
//                System.out.println(
//                    "Run number " + runNumber + " associated to params: Color(" +
//                    rgb[0] + " " + rgb[1] + " " + rgb[2] + ") Speed " + speed);
            }
            _runs.add(runNumber);

            // see whether we reached the end of the file
            st.nextToken();
            if (st.ttype == st.TT_EOF) {
                break;
            } else {
                st.pushBack();
            }
        }

//        System.out.println("nDifferentRuns " + nDifferentRuns);
//        System.out.println("runs " + runs.size());

        runNumberList = _runs.toArray();
        r.close();
    }


    /**
     *
     * @param id int
     * @param dt double in seconds
     * @throws IOException
     */
    public void showResponse(int id, double dt) throws IOException {
        DoubleHistogram[] barResponseUD = getMovingFilterResponse(id, dt);

        double max = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < barResponseUD.length; i++) {
            if (barResponseUD[i].getMaxValue() > max) {
                max = barResponseUD[i].getMaxValue();
            }
        }

        JFrame f = new JFrame("" + id);
        f.getContentPane().setLayout(new GridLayout(2, 1));

        for (int i = 0; i < barResponseUD.length; i++) {
            PlotPanel p = new PlotPanel();
            double[] par = paramsMap.get(i);
            p.addToLegend("Color: " + par[0] + ", Speed: " + par[3]);
            p.addData(barResponseUD[i], new HistogramStyle());
            p.autoscale();
            //            p.setYRange(0, max);
            p.setYRange(0, 250);
            f.add(p);
        }

        f.setBounds(10, 10, 300, 300);
        f.setVisible(true);
    }


    public static void showResponse1(int id, MovingBar[] b) throws IOException {
        JFrame f = new JFrame("" + id);
        f.getContentPane().setLayout(new GridLayout(0, 2));

        double max = Double.NEGATIVE_INFINITY;
        for (int k = 0; k < b.length; k++) {
            double dt = b[k].runLength[0] / 20000 / 500 /*bins*/;
            DoubleHistogram[] barResponseUD = b[k].getMovingFilterResponse(id, dt);
            for (int i = 0; i < barResponseUD.length; i++) {
                if (barResponseUD[i].getMaxValue() > max) {
                    max = barResponseUD[i].getMaxValue();
                }
            }
        }

        for (int k = 0; k < b.length; k++) {
            double dt = b[k].runLength[0] / 20000 / 500 /*bins*/;
            DoubleHistogram[] barResponseUD = b[k].getMovingFilterResponse(id, dt);

            for (int i = 0; i < barResponseUD.length; i++) {
                PlotPanel p = new PlotPanel();
                double[] par = b[k].paramsMap.get(i);
                p.addToLegend("Color: " + par[0] + ", Speed: " + par[3]);
                p.addData(barResponseUD[i], new HistogramStyle());
                p.autoscale();
                p.setYRange(0, max);
                f.add(p);
            }
        }

        f.setBounds(10, 10, 600, 800);
        f.setVisible(true);
    }


    public static double getPeak(DoubleHistogram h, int n) {
        int k = h.getMaxValueBin();
        double v = 0;
        for (int i = k - n; i <= k + n; i++) {
            v += h.getBin(i);
        }
        return v / (2 * n + 1);
    }


    public static void showTunning(int[] id, MovingBar[] b, String name) throws
        IOException {
        JFrame f = new JFrame("" + id);
        f.getContentPane().setLayout(new GridLayout(0, 2));
        PlotPanel p = new PlotPanel();

        for (int i = 0; i < id.length; i++) {
            ScatterPlot spOFF = new ScatterPlot();
            ScatterPlot spON = new ScatterPlot();
            for (int k = 0; k < b.length; k++) {
                double dt = b[k].runLength[0] / 20000 / 500 /*bins*/;
                DoubleHistogram[] barResponseUD = b[k].getMovingFilterResponse(id[i], dt);
                spON.add(Math.pow(2, k + 1), getPeak(barResponseUD[0], nb));
                spOFF.add(Math.pow(2, k + 1), getPeak(barResponseUD[1], nb));
            }

            p.addData(spOFF, new ScatterPlotStyle(SQUARE, 4, black, true, black, 2));
            p.addData(spON, new ScatterPlotStyle(SQUARE, 4, red, true, red, 2));
        }

        p.autoscale();
        p.padY();
        p.setXRange(1.8, 35);
        p.setAxesType(AxisType.LOG10, AxisType.LINEAR);
        PlotUtil.showData(name, p);
    }


    public static void showAverageTunning(int[] id, MovingBar[] b, String name) throws
        IOException {
        JFrame f = new JFrame("" + id);
        f.getContentPane().setLayout(new GridLayout(0, 2));
        PlotPanel p = new PlotPanel();

        ScatterPlot spOFF = new ScatterPlot();
        ScatterPlot spON = new ScatterPlot();
        for (int k = 0; k < b.length; k++) {
            double on = 0, off = 0;
            for (int i = 0; i < id.length; i++) {
                double dt = b[k].runLength[0] / 20000 / 500 /*bins*/;
                DoubleHistogram[] barResponseUD = b[k].getMovingFilterResponse(id[i], dt);
                on += getPeak(barResponseUD[0], nb);
                off += getPeak(barResponseUD[1], nb);
            }
            spON.add(Math.pow(2, k + 1), on / id.length);
            spOFF.add(Math.pow(2, k + 1), off / id.length);
        }

        p.addData(spOFF, new ScatterPlotStyle(SQUARE, 4, black, true, black, 2));
        p.addData(spON, new ScatterPlotStyle(SQUARE, 4, red, true, red, 2));
        p.autoscale();
        p.padY();
        p.setXRange(1.8, 35);
        p.setAxesType(AxisType.LOG10, AxisType.LINEAR);
        PlotUtil.showData(name, p);
    }


    /**
     *
     * @param id int
     * @param dt double time binning in seconds
     * @return DoubleHistogram[]
     * @throws IOException
     */
    public DoubleHistogram[] getMovingFilterResponse(int id, double dt) throws
        IOException {

        // create the run histograms
        DoubleHistogram[] runHistograms = new DoubleHistogram[nDifferentRuns];
        for (int i = 0; i < nDifferentRuns; i++) {
            runHistograms[i] = new DoubleHistogram("" + i, 0, runLength[i] / samplingRate,
                dt);
        }

        // process the spikes
        SpikeStream stream = new SpikeStream(nf.getSpikeTimes(id), ttl);
        boolean firstTTL = false;
        int currentTTL = -1;
        int currentRun = -1;
        int t;
        while ( (t = stream.getNext()) != Integer.MAX_VALUE) {
            if (t < 0) {
                currentTTL = -t;

                if (firstTTL) {
                    currentRun = 0;
                } else {
                    currentRun++;
                    if (currentRun >= totalRuns) {
                        break;
                    }
                }
            } else {
                if (currentRun != -1) {
                    runHistograms[runNumberList[currentRun]].fill(
                        (t - currentTTL) / samplingRate, 1);
                }
            }
        }

        // normalize the run histograms
        for (int i = 0; i < nDifferentRuns; i++) {
            runHistograms[i].scale(1.0 / (nRuns[i] * dt));
        }

        return runHistograms;
    }


    public static int[][] loadMovingFilterStimuli(String[] fileName, HashMap paramsMap) throws
        IOException {

        final boolean automatic;
        if (paramsMap == null) {
            paramsMap = new HashMap();
            automatic = true;
        } else {
            automatic = false;
        }

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
            if (st[fileIndex].ttype != st[fileIndex].TT_WORD ||
                !st[fileIndex].sval.equals("TYPE")) {
                throw new IOException("Missing TYPE.");
            }
            st[fileIndex].nextToken();
            if (st[fileIndex].ttype != st[fileIndex].TT_WORD) {
                throw new IOException("Wrong TYPE.");
            }
            String type = st[fileIndex].sval;
//            System.out.println("TYPE : " + type);

            // skip to the end of line
            do {
                st[fileIndex].nextToken();
            } while (st[fileIndex].ttype != st[fileIndex].TT_EOL);

            // read the run lines
            while (true) {
                // read the color of the bar
                st[fileIndex].nextToken();
                if (st[fileIndex].ttype != st[fileIndex].TT_WORD ||
                    !st[fileIndex].sval.equals("RGB")) {
                    throw new IOException("Missing RGB at line " + st[fileIndex].lineno());
                }
                double[] rgb = new double[3];
                for (int i = 0; i < 3; i++) {
                    st[fileIndex].nextToken();
                    if (st[fileIndex].ttype != st[fileIndex].TT_NUMBER) {
                        throw new IOException("Bad RGB at line " + st[fileIndex].lineno());
                    }
                    rgb[i] = st[fileIndex].nval;
                }

                // read the direction of motion of the bar
                st[fileIndex].nextToken();
                if (st[fileIndex].ttype != st[fileIndex].TT_WORD) {
                    if (!st[fileIndex].sval.equals("X-DELTA") ||
                        !st[fileIndex].sval.equals("Y-DELTA")) {
                        throw new IOException(
                            "Missing X-DELTA/Y-DELTA at line " + st[fileIndex].lineno());
                    }
                }
                st[fileIndex].nextToken();
                if (st[fileIndex].ttype != st[fileIndex].TT_NUMBER) {
                    throw new IOException(
                        "Missing X-DELTA at line " + st[fileIndex].lineno());
                }
                double speed = st[fileIndex].nval;

                // read the EOL
                st[fileIndex].nextToken();
                if (st[fileIndex].ttype != st[fileIndex].TT_EOL) {
                    throw new IOException("Too many values at line " +
                                          st[fileIndex].lineno());
                }

                // figure out what kind of run we have
                int runNumber = -1;
                for (int i = 0; i < paramsMap.size(); i++) {
                    double[] params = (double[]) paramsMap.get("" + i);
                    if (rgb[0] == params[0] &&
                        rgb[1] == params[1] &&
                        rgb[2] == params[2] &&
                        speed == params[3]) {

                        runNumber = i;
                        break;
                    }
                }
                if (runNumber == -1) {
                    if (automatic) {
                        runNumber = paramsMap.size();
                        double[] params = {rgb[0], rgb[1], rgb[2], speed};
                        paramsMap.put("" + runNumber, params);
                        System.out.println(
                            "Run number " + runNumber +
                            " associated to params: Color(" +
                            rgb[0] + " " + rgb[1] + " " + rgb[2] + ") Speed " +
                            speed);
                    } else {
                        throw new IllegalArgumentException("Unknown params");
                    }
                }
                runs[fileIndex].add(runNumber);

                // see whether we reached the end of the file
                st[fileIndex].nextToken();
                if (st[fileIndex].ttype == st[fileIndex].TT_EOF) {
                    break;
                } else {
                    st[fileIndex].pushBack();
                }
            }
        }

//        System.out.println("nDifferentRuns " + nDifferentRuns);
//        System.out.println("runs " + runs.size());

        int[][] allRuns = new int[fileName.length][];
        for (int fileIndex = 0; fileIndex < fileName.length; fileIndex++) {
            allRuns[fileIndex] = runs[fileIndex].toArray();
            r[fileIndex].close();
        }

        return allRuns;
    }

}
