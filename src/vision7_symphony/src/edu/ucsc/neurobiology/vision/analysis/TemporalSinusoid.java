package edu.ucsc.neurobiology.vision.analysis;

import static java.lang.Math.exp;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import edu.ucsc.neurobiology.vision.io.NeuronFile;
import edu.ucsc.neurobiology.vision.io.ParametersFile;
import edu.ucsc.neurobiology.vision.math.FFT;
import edu.ucsc.neurobiology.vision.math.MathUtil;
import edu.ucsc.neurobiology.vision.math.fitting.FittableFunction;
import edu.ucsc.neurobiology.vision.neuronviewer.SpikeRatePlotMaker;
import edu.ucsc.neurobiology.vision.plot.DoubleHistogram;
import edu.ucsc.neurobiology.vision.plot.FunctionData;
import edu.ucsc.neurobiology.vision.plot.PlotPanel;
import edu.ucsc.neurobiology.vision.plot.PlotUtil;
import edu.ucsc.neurobiology.vision.plot.ScatterPlot;


/**
 *
 * @author Dumitru Petrusca Grivich, University of California, Santa Cruz
 */
public class TemporalSinusoid {
    double temporalPeriod = 15;
    double t2 = 60;

    private NeuronFile neuronFile;
    private final int N = 8 * 4096;
    private double dt, omega;
    private double[] real, img;
    private int[] currentID;
    private double contrast;
    private double stimulusLength_fr;
    String experimentName, datasetName;


    public TemporalSinusoid(String datasetFolder) throws IOException {
        experimentName = new File(datasetFolder).getParentFile().getName();
        datasetName = new File(datasetFolder).getName();
        neuronFile = new NeuronFile(
            datasetFolder + File.separator + datasetName + ".neurons");

        real = new double[N];
        img = new double[N];
    }


    public double getContrast() {
        return contrast;
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
            if (i >= 0 && real[i] > max) {
                max = real[i];
            }
        }
        return max;
    }


    /**
     *
     * @param nBinsPerPeriod int number of bins in one period
     */
    public void setBinning(int nBinsPerPeriod) {
        double T = 1000 * temporalPeriod / 120; // temporal period in ms
        this.dt = T / nBinsPerPeriod;
        System.err.println("dt = " + dt + " ms");
        omega = 1.0 / (N * dt / 1000);
    }


    public void setCurrentNeurons(int ...id) throws IOException {
        this.currentID = id;
        Arrays.fill(real, 0);
        Arrays.fill(img, 0);

        final double T = dt / 1000; // temporal period in seconds

        // get the histogram
        DoubleHistogram h = SpikeRatePlotMaker.getSpikeRateHistogram(
            neuronFile, currentID, 0, t2, T);
        double[] x = h.toArray();

        // subtract the mean
        double xm = MathUtil.mean(x);
        MathUtil.sub(x, xm);

        for (int j = 0; j < x.length; j++) {
            real[j] = x[j];
            img[j] = 0;
        }

        FFT.fft(real, img, -1);

        for (int j = 0; j < N / 2; j++) {
            double amp =
                Math.sqrt(real[j] * real[j] + img[j] * img[j]) * T;
            double phase = Math.atan(img[j] / real[j]);
            real[j] = amp;
            img[j] = phase;
        }
    }


    /*
        public void showOnePeriodPSTH() throws IOException {
            DoubleHistogram hh = getOnePeriodHistograms(1);
            double max = hh.getMaxValue();
//        System.err.println(max);

            PlotPanel p = new PlotPanel();
            p.addData(hh, new HistogramStyle(
                "", OutlineType.RECTANGULAR, black, 1, true, black, false, black, 1));
            p.autoscale();
            p.setYRange(0, max);

            PlotUtil.showData("", p);
        }
     */



    private GaussianMixture fit(int runID) {
        // create the Gaussian Mixture
        double[] xi = new double[3 * 2];
        double[] ai = new double[3 * 2];
        for (int f = 0; f < xi.length; f++) {
            xi[f] = (f + 1) * 120 / temporalPeriod;
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
        ScatterPlot s = new ScatterPlot();
        for (int ii = 0; ii < N / 2; ii++) {
            s.add(ii * omega, real[ii]);
        }

        PlotPanel p = new PlotPanel();
        p.addData(s, "NONE 0 black, SOLID 1 black");
        p.autoscale();
        p.padY();

        PlotUtil.showData("", p);
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


    /*
        public DoubleHistogram getOnePeriodHistograms(int numberOfPeriods) throws
            IOException {

            final double T = dt / 1000; // temporal period in seconds
            DoubleHistogram h = SpikeRatePlotMaker.getSpikeRateHistogram(
                neuronFile, new int[] {currentID}, 0, t2, T);
                System.err.println(h.getBinCount());

            double averagingPeriod_s = numberOfPeriods * (temporalPeriod / 120);
     int binsPerAveragingPeriod = (int) Math.round(averagingPeriod_s / (dt / 1000));
            int totalPeriods = (int) Math.floor(
                (stimulusLength_fr / 120) / averagingPeriod_s);

     DoubleHistogram hh = new DoubleHistogram("", 0, averagingPeriod_s, dt / 1000);
            for (int binID = 0; binID < totalPeriods * binsPerAveragingPeriod; binID++) {
                hh.fillBin(
                    (int) (binID % binsPerAveragingPeriod), h.getBin(binID));
            }
//        hh.scale( (double) numberOfPeriods / totalPeriods);

            return hh;
        }
     */


    public void simulateSpectra() {
        Arrays.fill(real, 0);
        Arrays.fill(img, 0);
        final double T = dt / 1000; // temporal period in seconds

        // get the histogram
        DoubleHistogram h = new DoubleHistogram("", 0, t2, T);
        for (int i = 0; i < h.getBinCount(); i++) {
            double v = Math.sin(2 * Math.PI * 8 * (i + 0.5) * T);
            double c = 0.75;
            if (v < 0) {
                v = 0;
            } else if (v > c) {
                v = c;
            }
            h.setBin(i, v);
        }
        double[] x = h.toArray();

        // subtract the mean
        MathUtil.sub(x, MathUtil.mean(x));

        for (int j = 0; j < x.length; j++) {
            real[j] = x[j];
            img[j] = 0;
        }

        FFT.fft(real, img, -1);

        for (int j = 0; j < N / 2; j++) {
            double amp = Math.sqrt(real[j] * real[j] + img[j] * img[j]) * T;
            double phase = Math.atan(img[j] / real[j]);
            real[j] = amp;
            img[j] = phase;
        }
    }


    public static void main(String[] args) throws Exception {
        ParametersFile f = new ParametersFile(
            "F:\\Good Data\\2006-02-10-0\\data010\\data010.params");
        TemporalSinusoid t = new TemporalSinusoid("F:\\Good Data\\2006-02-10-0\\data010");
        t.setBinning(50);

//        t.setCurrentNeurons(f.getNeuronsInClass("All/ON/Midget"));
//        t.showSpectra();
//        t.setCurrentNeurons(f.getNeuronsInClass("All/OFF/Midget"));
//        t.showSpectra();
//        t.setCurrentNeurons(f.getNeuronsInClass("All/ON/Parasol"));
//        t.showSpectra();
        t.setCurrentNeurons(f.getNeuronsInClass("All/OFF/Parasol"));
        t.showSpectra();
//        t.setCurrentNeurons(f.getNeuronsInClass("All/OFF/Large"));
//        t.showSpectra();
//        t.setCurrentNeurons(f.getNeuronsInClass("All/ON/nc4"));
//        t.showSpectra();

        t.simulateSpectra();
        t.showSpectra();
    }
}
