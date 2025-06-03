package edu.ucsc.neurobiology.vision.analysis;

import java.io.*;

import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class DriftingDirectionSinusoids {
//    private HashMap paramsMap;
    NeuronFile neuronFile;
    DoubleList directions = new DoubleList();

    final double dt = 8.34;
    final double T = dt / 1000;
    public final double[] contrasts = {0.03, 0.06, 0.12, 0.24, 0.48};
    public final int nContrasts = contrasts.length;
    final int N = 2 * 1024;
    final double omega = 1.0 / (N * T);

//    double[][] maxHarmonic = new double[periods.length][8];
//    double[][] meanHarmonic = new double[periods.length][8];
//    double[][] a = new double[nContrasts][N];
//    double[][] b = new double[nContrasts][N];
    int currentID;
    double timeToIgnore = 0;
    double contrast, temporalPeriod;
    int[] times, ttl;


    public DriftingDirectionSinusoids(String datasetFolder) throws IOException {
        String datasetName = new File(datasetFolder).getName();
        neuronFile = new NeuronFile(
            datasetFolder + File.separator + datasetName + ".neurons");
        loadStimulus(
            new File(datasetFolder).getParent() + File.separator + "stimuli" +
            File.separator + datasetName + ".txt");
    }


    private void loadStimulus(String fileName) throws
        IOException {

        InputStreamReader r = new InputStreamReader(new FileInputStream(fileName));
        StreamTokenizer st = new StreamTokenizer(r);

        st.whitespaceChars('(', '(');
        st.whitespaceChars(')', ')');
        st.whitespaceChars(':', ':');
        st.whitespaceChars('#', '#');
        st.eolIsSignificant(true);

//        while (true) {
//            st.nextToken();
//            System.out.println(st.ttype);
//        }

        // read the stimulus type
        st.nextToken();
        if (st.ttype != StreamTokenizer.TT_WORD || !st.sval.equals("TYPE")) {
            throw new IOException("Missing TYPE.");
        }
        st.nextToken();
        if (st.ttype != StreamTokenizer.TT_WORD) {
            throw new IOException("Wrong TYPE.");
        }
        String type = st.sval;
        System.out.println("TYPE : " + type);

        // skip to the end of line
        do {
            st.nextToken();
        } while (st.ttype != StreamTokenizer.TT_EOL);

        // read the run lines
        while (true) {
            // skip spatial period and temporal period
            st.nextToken();
            st.nextToken();
            st.nextToken();
            st.nextToken();

            // read the direction of the drift
            st.nextToken();
            if (st.ttype != StreamTokenizer.TT_WORD || !st.sval.equals("DIRECTION")) {
                throw new IOException("Missing DIRECTION at line " + st.lineno());
            }
            st.nextToken();
            if (st.ttype != StreamTokenizer.TT_NUMBER) {
                throw new IOException("Missing DIRECTION at line " + st.lineno());
            }
            double direction = st.nval;

            // read the EOL
            st.nextToken();
            if (st.ttype != StreamTokenizer.TT_EOL) {
                throw new IOException("Too many values at line " + st.lineno());
            }

            directions.add(direction);
            System.out.println("direction " + direction);

            // see whether we reached the end of the file
            st.nextToken();
            if (st.ttype == StreamTokenizer.TT_EOF) {
                break;
            } else {
                st.pushBack();
            }
        }

        r.close();
    }


    public void setCurrentNeruon(int id) throws IOException {
        times = neuronFile.getSpikeTimes(id);
        ttl = neuronFile.getTTLTimes();
    }


    /**
     * The drifting sinusoid stimulus consists of a period of sinusoid floowed by a
     * period of gray. During the sinusoid period there are TTLs generated after every
     * small time interval

     * @param binning double
     * @param runIndex int
     * @return DoubleHistogram
     */
    public DoubleHistogram getResponse(double binning, int runIndex) {
        // remove the intermediary ttl's
        final int ttlsPerDirection = ttl.length / directions.size();
        System.out.println("ttlsPerDirection " + ttlsPerDirection);

//        double averageRunLength = 0;
//        // calculate the average length of a sinusoid run plus the following gray time
//        for (int run = 0; run < directions.size(); run++) {
//            if (run != 0) {
//                averageRunLength += ttl[run] - ttl[run - 1];
//            }
//        }
//        System.out.println(averageRunLength);
        // convert averageRunLength to seconds
//        averageRunLength /= directions.size() * 20000.0;

        int ttl1 = ttl[ (runIndex + 0) * ttlsPerDirection];
        int ttl2 = ttl[ (runIndex + 1) * ttlsPerDirection];

        // process the spikes
        DoubleHistogram h = new DoubleHistogram("", 0, (ttl2 - ttl1) / 20000.0, binning);
        for (int i = 0; i < times.length; i++) {
            if (times[i] >= ttl1 && times[i] < ttl2) {
                h.fill( (times[i] - ttl1) / 20000.0, 1);
            }
        }
        // average and convert to spike rate
        h.scale(1.0 / binning);

        return h;
    }


    public DoubleHistogram getResponse(double binning) {
        DoubleHistogram hh = null;
        for (int i = 0; i < directions.size() - 1; i++) {
            DoubleHistogram h = getResponse(binning, i);
            if (hh == null) {
                hh = h;
            } else {
                for (int j = 0; j < hh.getBinCount(); j++) {
                    hh.fillBin(j, h.getBin(j));
                }
            }
        }
        return hh;
    }


    public double getContrast() {
        return contrast;
    }


    public double getTemporalPeriod() {
        return temporalPeriod;
    }


    /*
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
     */

    public String getRunInfo(int i) {
        return "direction " + directions.get(i);
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


    public void showTTLs() throws IOException {
        int[] ttl = neuronFile.getTTLTimes();
        ScatterPlot s = new ScatterPlot();
        for (int i = 0; i < ttl.length; i++) {
            s.add(ttl[i], 1);
        }
        ScatterPlotStyle style = new ScatterPlotStyle();
        style.setSymbolType(SymbolType.VERTICAL_LINE);
        style.setSymbolSize(100);
        PlotUtil.showData(ttl.length + " TTLs", s, style).setRange(0,
            ttl[ttl.length - 1] * 1.1, 0, 1);

        int[] ttlDistance = new int[ttl.length - 1];
        for (int i = 0; i < ttl.length - 1; i++) {
            ttlDistance[i] = ttl[i + 1] - ttl[i];
        }
        DoubleHistogram h = new DoubleHistogram("", 0, MathUtil.max(ttlDistance) * 1.1, 1);
        for (int i = 0; i < ttlDistance.length; i++) {
            h.fill(ttlDistance[i], 1);
        }
        PlotUtil.showData("TTL spacing hisogram", h, new HistogramStyle());
    }

}
