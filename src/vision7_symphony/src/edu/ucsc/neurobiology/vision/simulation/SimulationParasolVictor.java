package edu.ucsc.neurobiology.vision.simulation;

import java.io.*;
import java.util.*;

import java.awt.*;
import javax.swing.*;

import edu.ucsc.neurobiology.vision.analysis.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.math.fitting.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.stimulus.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class SimulationParasolVictor
    extends Function implements FunctionMinimization.IterationInfo {

    boolean DEBUG = false;
    boolean USE_CNTRAST_GAIN = false;

    AdvancedMovie[] movie;
    double[] refreshTime;
    double stixelWidth, stixelHeight, stixelArea;
    double pSize;
    int width, height;
    int nData;
    double[][] dataON /*, dataOFF*/;
    int[] dataLength;
    DoubleHistogram[] dataHistON /*, dataHistOFF*/;
    double[][] currentON /*, currentOFF*/;
    double[][] simulatedDataON /*, simulatedDataOFF*/;
    STATimeFunction1 subunitCenterTCFunction, subunitSurroundTCFunction;
    final double staLength = 400; // in ms
    double[] parameters;
    double[] freeParameters;
    boolean[] freeParameter;

    double dtReversingGratings;
    double cgcT;
    double cgcB;
    double nlfA, nlfI0, nlfB;
    double x0, y0, dyReversingGratings, dyDriftingSinusoinds;
    double centerA;
    double surroundA;
    double centerSigmaX;
    double centerSigmaY;
    double surroundSigmaX;
    double surroundSigmaY;

    double centerTCt1;
    double centerTCt2;
    double centerTCa;
    double centerTCn1;
    double centerTCn2;

    double surroundTCt1;
    double surroundTCt2;
    double surroundTCa;
    double surroundTCn1;
    double surroundTCn2;

    JFrame window = new JFrame();
    PlotPanel[] currentPanel;
    PlotPanel layoutPanel = new PlotPanel();
    FunctionStyle st = new FunctionStyle("f");

    ScatterPlotStyle blackThin = new ScatterPlotStyle(
        SymbolType.NONE, 1, Color.black, true, Color.black, 1);

    ScatterPlotStyle redThin = new ScatterPlotStyle(
        SymbolType.NONE, 1, Color.black, true, Color.red, 1);
    ScatterPlotStyle redThick = new ScatterPlotStyle(
        SymbolType.NONE, 1, Color.black, true, Color.red, 2);

    ScatterPlotStyle blueThin = new ScatterPlotStyle(
        SymbolType.NONE, 1, Color.black, true, Color.blue, 1);
    ScatterPlotStyle blueThick = new ScatterPlotStyle(
        SymbolType.NONE, 1, Color.black, true, Color.blue, 2);

    ScatterPlotStyle style3 = new ScatterPlotStyle(
        SymbolType.NONE, 1, Color.black, true, Color.green, 2);


    public SimulationParasolVictor(Object[][] input, double[][] x) throws IOException {
        movie = new AdvancedMovie[input.length];
        DoubleHistogram[] dON = new DoubleHistogram[input.length];
//        DoubleHistogram[] dOFF = new DoubleHistogram[input.length];

        for (int i = 0; i < input.length; i++) {
            movie[i] = (AdvancedMovie) input[i][0];
            dON[i] = (DoubleHistogram) input[i][1];
//            dOFF[i] = (DoubleHistogram) input[i][2];
        }

        initialize(movie, dON);

        parameters = new double[x.length];
        freeParameter = new boolean[x.length];
        DoubleList changeP = new DoubleList();
        for (int i = 0; i < parameters.length; i++) {
            parameters[i] = x[i][1];
            if (x[i][0] == 1) {
                freeParameter[i] = true;
                changeP.add(x[i][1]);
            } else {
                freeParameter[i] = false;
            }
        }
        freeParameters = changeP.toArray();
    }


    public void simulate() throws CannotEvaluateException {
        final double[][] xi = new double[freeParameters.length][freeParameters.length];
        for (int i = 0; i < freeParameters.length; i++) {
            xi[i][i] = 1;
        }

        FunctionMinimization.powell(freeParameters, xi, 1e-6, this, this);
    }


    private void initialize(AdvancedMovie[] movie, DoubleHistogram[] dON) throws
        IOException {
        window.getContentPane().setLayout(new GridLayout(0, 2));
        layoutPanel.setRange(0, 3.2, 0, 1.6);

        this.movie = movie;
        this.dataHistON = dON;
//        this.dataHistOFF = dOFF;
        this.nData = movie.length;
        this.dataON = new double[nData][];
//        this.dataOFF = new double[nData][];
        this.currentON = new double[nData][];
//        this.currentOFF = new double[nData][];
//        this.simulatedDataOFF = new double[nData][];
        this.simulatedDataON = new double[nData][];
        this.refreshTime = new double[nData];
        this.currentPanel = new PlotPanel[nData];
        this.width = movie[0].getWidth();
        this.height = movie[0].getHeight();
        this.pSize = movie[0].getFrame(0).getStixelWidth();
        this.stixelWidth = movie[0].getFrame(0).getStixelWidth();
        this.stixelHeight = movie[0].getFrame(0).getStixelHeight();
        this.stixelArea = stixelWidth * stixelHeight;
        this.dataLength = new int[nData];

        for (int n = 0; n < nData; n++) {
            refreshTime[n] = movie[n].getRefreshTime();
            currentPanel[n] = new PlotPanel();

            this.dataON[n] = dON[n].toArray();
//            if (dOFF[n] != null) {
//                this.dataOFF[n] = dOFF[n].toArray();
//            }
            dataLength[n] = dataON[n].length;
            System.out.println("dataLength[n] " + dataLength[n] + ", " +
                               dON[n].getBinCount());

//            currentOFF[n] = new double[dataLength[n]];
            currentON[n] = new double[dataLength[n]];
            simulatedDataON[n] = new double[dataLength[n]];
//            simulatedDataOFF[n] = new double[dataLength[n]];
        }
    }


    public void receiveInfo(String info) {
        System.out.println(info);
        for (int k = 0; k < parameters.length; k++) {
            int free = freeParameter[k] ? 1 : 0;
            System.out.println("{" + free + ", " + parameters[k] + "},");
        }
    }


    public double getValue() throws CannotEvaluateException {
        return getValue(freeParameters);
    }


    public double getValue(final double[] p) throws CannotEvaluateException {
        for (int i = 0, j = 0; i < parameters.length; i++) {
            if (freeParameter[i]) {
                parameters[i] = p[j];
                j++;
            }
        }
        readParameters(parameters);

        // construct the needed functions
        subunitCenterTCFunction = new STATimeFunction1(
            1, centerTCt1, centerTCa, centerTCt2, centerTCn1, centerTCn2);
        subunitSurroundTCFunction = new STATimeFunction1(
            1, surroundTCt1, surroundTCa, surroundTCt2, surroundTCn1, surroundTCn2);

        // construct the STA of the middle subunit
//        constructSTA(x0, y0, subunitSTA);

        // calculate the currents
        for (int n = 0; n < nData; n++) {
            // clear the currents
            for (int i = 0; i < dataLength[n]; i++) {
                currentON[n][i] = 0;
//                currentOFF[n][i] = 0;
            }

            if (movie[n] instanceof MovingBarMovie) {
                calculateMovingFilterCurrent(n);
//            } else if (movie[n] instanceof FullScreenFlashMovie) {
//                calculateFullScreenFlashCurrent(n);
            } else if (movie[n] instanceof DriftingSinusoidMovie) {
                calculateDriftingSinusoidCurrent(n);
//            } else if (movie[n] instanceof ReversingGratingMovie) {
//                calculateReversingGratingsCurrent(n);
            } else {
                throw new IllegalArgumentException("Unknown movie: " +
                    movie[n].getClass());
            }
        }

        // generate the spike rate
        for (int n = 0; n < nData; n++) {
            for (int i = 0; i < dataLength[n]; i++) {
//                simulatedDataON[n][i] =
//                    nlfA / (1 + StrictMath.exp( - (currentON[n][i] - nlfI0) / nlfSigma));
//                simulatedDataOFF[n][i] =
//                    nlfA / (1 + StrictMath.exp( - (currentOFF[n][i] - nlfI0) / nlfSigma));

                double dx = currentON[n][i] - nlfI0;
                simulatedDataON[n][i] = (dx > 0) ? nlfB + nlfA * Math.pow(dx, 1) : nlfB;
            }
        }

        if (DEBUG) {
            showResult();
        }

        // calculate the likelihood
        double l = 0;
        for (int n = 0; n < nData; n++) {
            for (int i = 0; i < dataLength[n]; i++) {
                l += simulatedDataON[n][i] -
                    dataON[n][i] * Math.log(simulatedDataON[n][i]);
            }
        }
        return l;
    }


    /*
        void constructSTA(double x, double y, float[][] sta) {
            double a, b, t, xx, yy, subunitCenterG, subunitSurroundG;
            Gaussian2DFunction subunitCenterGaussian = new Gaussian2DFunction(
                centerA, x, y, centerSigmaX, centerSigmaY, 0, 0);
            Gaussian2DFunction subunitSurroundGaussian = new Gaussian2DFunction(
                surroundA, x, y, surroundSigmaX, surroundSigmaY, 0, 0);
            for (int f = 0; f < staDepth; f++) {
                a = - (staDepth - 1 - f + 1) * refreshTime;
                b = - (staDepth - 1 - f + 0) * refreshTime;
                t = - (staDepth - 1 - f + 0.5) * refreshTime;
                    subunitCenterTC[f] = Integrator.qsimp(subunitCenterTCFunction, a, b,
                        1e-5);
         subunitSurroundTC[f] = Integrator.qsimp(subunitSurroundTCFunction, a, b,
                        1e-5);
                for (int i = 0; i < width; i++) {
                    sta[f][3 * i + 0] = 0;
                    sta[f][3 * i + 1] = 0;
                    sta[f][3 * i + 2] = 0;
                    for (int j = 0; j < height; j++) {
                        xx = (i + 0.5) * stixelSize;
                        yy = (j + 0.5) * stixelSize;
                        subunitCenterG = subunitCenterGaussian.getValue(xx, yy);
                        subunitSurroundG = subunitSurroundGaussian.getValue(xx, yy);
                        double c = (subunitCenterG * subunitCenterTC[f] +
         subunitSurroundG * subunitSurroundTC[f]) * stixelArea;
                        final int index = (j * width + i) * 3;
                        staFrames[f][index + 0] = (float) c;
                        staFrames[f][index + 1] = (float) c;
                        staFrames[f][index + 2] = (float) c;
                        sta[f][3 * i + 0] += (float) c;
                        sta[f][3 * i + 1] += (float) c;
                        sta[f][3 * i + 2] += (float) c;
                    }
                }
            }
        }
     */

    double t0 = 0;
    int i0 = 1;

    void calculateMovingFilterCurrent(int n) throws CannotEvaluateException {
        MovingBarMovie m = ( (MovingBarMovie) movie[n]);
        final int angle = (int) Math.round(m.getAngle() * 180 / Math.PI);
        final double velocity = m.getVelocity() * pSize / refreshTime[n];
        final double w = m.getThickness() * pSize;
        final double I0 = m.getContrast();
        final double ac = centerSigmaX * centerSigmaX + w * w;
        final double as = surroundSigmaX * surroundSigmaX + w * w;
        final double x1;
        if (angle == 0) {
            x1 = m.getX1() * pSize;
        } else if (angle == 90 || angle == -90) {
            x1 = m.getY1() * pSize;
        } else {
            System.out.println(angle);
            throw new Error("only 0 and 90 deg bars accepted");
        }
        final double Ac = 2 * Math.PI * I0 * centerA *
                          centerSigmaX * centerSigmaY * w / Math.sqrt(ac);
        final double As = 2 * Math.PI * I0 * surroundA *
                          surroundSigmaX * surroundSigmaY * w / Math.sqrt(as);

        FunctionDataAdapter f = new FunctionDataAdapter() {
            public double getValueAt(double t) {
                double xb = x1 + velocity * t;
                double c = 0;
                if (angle == 0) {
                    c += Ac * StrictMath.exp( - (x0 - xb) * (x0 - xb) / (2 * ac)) *
                        subunitCenterTCFunction.getValueAt(t - t0);
                    c += As * StrictMath.exp( - (x0 - xb) * (x0 - xb) / (2 * as)) *
                        subunitSurroundTCFunction.getValueAt(t - t0);
                } else {
                    c += Ac * StrictMath.exp( - (y0 - xb) * (y0 - xb) / (2 * ac)) *
                        subunitCenterTCFunction.getValueAt(t - t0);
                    c += As * StrictMath.exp( - (y0 - xb) * (y0 - xb) / (2 * as)) *
                        subunitSurroundTCFunction.getValueAt(t - t0);
                }

                return c;
            }
        };

        double[] filter = new double[100];
        for (int i = 0; i < filter.length; i++) {
            final double a = - (filter.length - 1 - i + 1) * refreshTime[n];
            final double b = - (filter.length - 1 - i + 0) * refreshTime[n];
            filter[i] = Integrator.gauss8(f, a, b, 1e-4);
        }
//        PlotUtilities.showData(new ScatterPlot(filter), new ScatterPlotStyle());

        double[] c = new double[dataLength[n]];
        double[] g = new double[dataLength[n] + 1];
        double[] y = new double[dataLength[n]];
        Arrays.fill(g, 1);
        for (int i = 0; i < dataLength[n]; i++) {
            t0 = (i + 0) * refreshTime[n];

            if (USE_CNTRAST_GAIN) {
                y[i] = g[i] * Integrator.gauss8(f, t0 - staLength, t0, 1e-4);
                for (int j = 0; j <= i; j++) {
                    c[i] += Math.abs(y[j]) * Math.exp( - (i - j) * refreshTime[n] / cgcT);
                }
                c[i] *= (cgcB / cgcT) * refreshTime[n];
                g[i + 1] = 1 / (1 + Math.pow(c[i], 4));
            } else {
                y[i] = Integrator.gauss8(f, t0 - staLength, t0, 1e-4);
            }

            currentON[n][i] = y[i];
//            currentOFF[n][i] = -currentON[n][i];
        }

        if (DEBUG) {
//            PlotUtilities.showData("c", new ScatterPlot(c), blackThin).autoscale();
//            PlotUtilities.showData("g", new ScatterPlot(g), blackThin).autoscale().
//                setYRange(0, 1);
            currentPanel[n].addData(new ScatterPlot(currentON[n]), redThin);
//            currentPanel[n].addData(new ScatterPlot(currentOFF[n]), blueThin);
        }
    }


    void calculateDriftingSinusoidCurrent(int n) throws CannotEvaluateException {
        DriftingSinusoidMovie m = ( (DriftingSinusoidMovie) movie[n]);
        final double Ty = m.getSpatialPeriod() * pSize;
        final double T = m.getTemporalPeriod();
        final double wy = 2 * Math.PI / Ty;
        final double w = 2 * Math.PI / T;
        final double I0 = m.getContrast();
        final double cA = -2 * centerA * I0 * Math.PI *
                          centerSigmaX * centerSigmaY *
                          StrictMath.exp( -centerSigmaY * centerSigmaY * wy * wy / 2);
        final double sA = -2 * surroundA * I0 * Math.PI *
                          surroundSigmaX * surroundSigmaY *
                          StrictMath.exp( -surroundSigmaY * surroundSigmaY * wy * wy / 2);

        FunctionDataAdapter tc = new FunctionDataAdapter() {
            public double getValueAt(double t) {
                return (
                    cA * subunitCenterTCFunction.getValueAt(t - t0)
                    +
                    sA * subunitSurroundTCFunction.getValueAt(t - t0)
                    ) * Math.cos(w * t - wy * (y0 - dyDriftingSinusoinds));
            }
        };

        for (int i = 0; i < dataLength[n]; i++) {
            t0 = i * refreshTime[n];
            currentON[n][i] = Integrator.gauss8(tc, t0 - staLength, t0, 1e-5);
        }

        if (DEBUG) {
            currentPanel[n].addData(new ScatterPlot(currentON[n]), redThin);
        }
    }


    public static double heaviside(double x) {
        return (x < 0) ? 0 : 1;
    }


    void calculateFullScreenFlashCurrent(int n) throws CannotEvaluateException {
        double Io = 0.48;
        final double C = 2 * Math.PI * Io * centerA * centerSigmaX * centerSigmaY;
        final double S = 2 * Math.PI * Io * surroundA * surroundSigmaX * surroundSigmaY;
        FunctionDataAdapter f = new FunctionDataAdapter() {
            public double getValueAt(double t) {
                return (heaviside(t) - heaviside(t - 1000)) *
                    (C * subunitCenterTCFunction.getValueAt(t - t0) +
                     S * subunitSurroundTCFunction.getValueAt(t - t0));
            }
        };

        double[] c = new double[dataLength[n]];
        double[] g = new double[dataLength[n] + 1];
        Arrays.fill(g, 1);
        double[] y = new double[dataLength[n]];
        for (int i = 0; i < dataLength[n]; i++) {
            t0 = (i + 0) * refreshTime[n];

            if (USE_CNTRAST_GAIN) {
                y[i] = g[i] * Integrator.gauss8(f, t0 - staLength, t0, 1e-4);
                for (int j = 0; j <= i; j++) {
                    c[i] += Math.abs(y[j]) * Math.exp( - (i - j) * refreshTime[n] / cgcT);
                }
                c[i] *= (cgcB / cgcT) * refreshTime[n];
                g[i + 1] = 1 / (1 + Math.pow(c[i], 4));
            } else {
                y[i] = Integrator.gauss8(f, t0 - staLength, t0, 1e-4);
            }

            currentON[n][i] = y[i];
        }

        if (DEBUG) {
//            PlotUtilities.showData("c", new ScatterPlot(c), blackThin).autoscale();
//            PlotUtilities.showData("g", new ScatterPlot(g), blackThin).autoscale().
//                setYRange(0, 1);
            currentPanel[n].addData(new ScatterPlot(currentON[n]), redThin);
        }
    }


    /*
        void calculateReversingGratingsCurrent(int n) throws CannotEvaluateException {
            ReversingGratingMovie m = ( (ReversingGratingMovie) movie[n]);
            final double Ty = m.getSpatialPeriod() * pSize;
            final double phase = m.getPhase() * pSize;
            final double T = m.getTemporalPeriod();
            final double I0 = m.getContrast();
            final double wy = 2 * Math.PI / Ty;
//        System.out.println("wy " + wy);
            final double w = 2 * Math.PI / T;
            final double cA = 2 * centerA * I0 * Math.PI * centerSigmaX * centerSigmaY *
                              Math.exp( -centerSigmaX * centerSigmaY * wy * wy / 2);
     final double sA = 2 * surroundA * I0 * Math.PI * surroundSigmaX * surroundSigmaY *
                              Math.exp( -surroundSigmaX * surroundSigmaY * wy * wy / 2);
            FunctionDataAdapter C_SIN = new FunctionDataAdapter() {
                public double getValueAt(double t) {
                    return Math.sin(w * t) * subunitCenterTCFunction.getValueAt(t - t0);
                }
            };
            FunctionDataAdapter S_SIN = new FunctionDataAdapter() {
                public double getValueAt(double t) {
     return Math.sin(w * t) * subunitSurroundTCFunction.getValueAt(t - t0);
                }
            };
            double[] c_sin = new double[dataLength[n]];
            double[] s_sin = new double[dataLength[n]];
            for (int i = 0; i < dataLength[n]; i++) {
                t0 = i * refreshTime[n] - dtReversingGratings;
                c_sin[i] = Integrator.gauss8(C_SIN, t0 - staLength, t0, 1e-5);
                s_sin[i] = Integrator.gauss8(S_SIN, t0 - staLength, t0, 1e-5);
            }

            double sin = Math.sin(wy * (y0 - dyReversingGratings + phase));
            double[] c = new double[dataLength[n]];
            double[] g = new double[dataLength[n] + 1];
            double[] y = new double[dataLength[n]];
            g[0] = 1;
            for (int i = 0; i < dataLength[n]; i++) {
                t0 = (i + 0) * refreshTime[n];

                if (USE_CNTRAST_GAIN) {
                    y[i] = g[i] * (cA * c_sin[i] + sA * s_sin[i]) * sin;
                    for (int j = 0; j <= i; j++) {
     c[i] += Math.abs(y[j]) * Math.exp( - (i - j) * refreshTime[n] / cgcT);
                    }
                    c[i] *= (cgcB / cgcT) * refreshTime[n];
                    g[i + 1] = 1 / (1 + Math.pow(c[i], 4));
                } else {
                    y[i] = (cA * c_sin[i] + sA * s_sin[i]) * sin;
                }

                currentON[n][i] = y[i];
            }

            if (DEBUG) {
//            PlotUtilities.showData("c", new ScatterPlot(c), blackThin).autoscale();
//            PlotUtilities.showData("g", new ScatterPlot(g), blackThin).autoscale().
//                setYRange(0, 1);
                currentPanel[n].addData(new ScatterPlot(currentON[n]), blackThin);
            }
        }
     */

    void showResult() {
        /*
                 // show the STA
                 ArrayList frameList = new ArrayList();
                 for (int i = 0; i < staDepth; i++) {
            frameList.add(new STAFrame(width, height, 1, staFrames[i]));
                 }
                 VisionUtilities.showSTA(
            "STA", new STA(frameList, refreshTime), 24, false);
         */
//        SigmoidFunction nlf = new SigmoidFunction(nlfA, nlfI0, nlfSigma);
//        PlotUtilities.showData("NLF", nlf, new FunctionStyle()).
//            setRange( -2, 5, 0, 200);

        System.out.println("C/S = " +
                           (centerA * centerSigmaX * centerSigmaY) /
                           (surroundA * surroundSigmaX * surroundSigmaY));
        System.out.println("s2/s1 = " + (surroundSigmaX / centerSigmaX));
        System.out.println("surroundA = " + surroundA);

        // show timecourses
        PlotPanel p = new PlotPanel();
        p.addData(subunitCenterTCFunction, new FunctionStyle("", Color.black, 2));
        p.addData(subunitSurroundTCFunction, new FunctionStyle("", Color.black, 1));
        p.setRange( -300, 0, -2, 2);
        PlotUtil.showData("Subunit TC", p);

        // show spacial filters
        PlotPanel spacePanel = new PlotPanel();
        spacePanel.addData(
            new DOG1DFunction(centerA, surroundA, x0, centerSigmaX,
                              surroundSigmaX), new FunctionStyle("", Color.red, 1));
        spacePanel.addData(
            new DOG1DFunction(centerA, surroundA, x0, centerSigmaY,
                              surroundSigmaY), new FunctionStyle("", Color.blue, 1));
        spacePanel.setRange(0, 3.2, -0.5, 5);
        PlotUtil.showData("Receptive Field", spacePanel);

        // show results
        for (int n = 0; n < nData; n++) {
//            currentPanel[n].addData(new ScatterPlot( currentON[n]), redThick);
//            currentPanel[n].addData(new ScatterPlot( currentOFF[n]), blueThick);
            currentPanel[n].autoscale();
            currentPanel[n].setYRange( -5, 5);
            currentPanel[n].setAxisVisible(false);
            window.add(currentPanel[n]);

            PlotPanel hPanelON = new PlotPanel();
            hPanelON.addData(dataHistON[n], new HistogramStyle());
            ScatterPlot spON = new ScatterPlot("");
            for (int j = 0; j < dataLength[n]; j++) {
                spON.add(j * refreshTime[n] / 1000, simulatedDataON[n][j]);
            }
            hPanelON.addData(spON, redThin);
            hPanelON.autoscale();
            hPanelON.setAxisVisible(false);
            hPanelON.setYRange(0, 100);
            window.add(hPanelON);
            /*
                        if (dataHistOFF[n] != null) {
                            PlotPanel hPanelOFF = new PlotPanel();
                            hPanelOFF.addData(dataHistOFF[n], new HistogramStyle());
                            ScatterPlot spOFF = new ScatterPlot("");
                            for (int j = 0; j < dataLength[n]; j++) {
                 spOFF.add(j * refreshTime[n] / 1000, simulatedDataOFF[n][j]);
                            }
                            hPanelOFF.addData(spOFF, blackThin);
                            hPanelOFF.autoscale();
                            hPanelOFF.setYRange(0, 170);
                            window.add(hPanelOFF);
                        } else {
                            window.add(new JLabel("Empty"));
                        }
             */
        }
//        window.setBounds(0, 0, 1280, 1000);
        window.setBounds(0, 0, 800, 1000);
        window.setVisible(true);
    }


    /*
        public double[] secondFilter(double[] f, double[] x, double dt) {
            double[] y = new double[x.length];
            double[] c = new double[x.length];
            double[] g = new double[x.length];
            g[0] = 1;
//       for (int i = 0; i < x.length; i++) {
//            int j0 = i - f.length + 1;
//            for (int j = Math.max(j0, 0); j < i; j++) {
//                y[i] += x[j] * f[j - j0];
//            }
//        }
            for (int i = 0; i < x.length - 1; i++) {
                int j0 = i - f.length + 1;
                for (int j = Math.max(j0, 0); j < i; j++) {
                    y[i] += x[j] * f[j - j0];
                }
                y[i] *= g[i];
                for (int j = 0; j < i; j++) {
                    c[i] += Math.abs(y[j]) * Math.exp( - (i - j) * dt / cgcT);
                }
                c[i] *= (cgcB / cgcT) * dt;
                g[i + 1] = 1 / (1 + Math.pow(c[i], 4));
            }
            PlotUtilities.showData("c", new ScatterPlot(c), blackThin).autoscale();
            PlotUtilities.showData("g", new ScatterPlot(g), blackThin).autoscale();
            return y;
        }
     */

    public static double[] dotProduct(double[] filter, double[] data) {
        double[] result = new double[data.length];

        for (int n = 0; n < data.length; n++) {
            int firstIndex = n - filter.length + 1;
            for (int i = 0; i < filter.length; i++) {
                if (firstIndex + i >= 0) {
                    result[n] += filter[i] * data[firstIndex + i];
                }
            }
        }

        return result;
    }


    public static void dotProduct1(float[][] staFrames, float[][] movie,
                                   double[] current) {
        for (int n = 0; n < movie.length; n++) {
            current[n] = 0;

            // calculate the instantaneous current
            int firstFrame = n - staFrames.length + 1;
            for (int f = 0; f < staFrames.length; f++) {
                if (firstFrame + f >= 0) {
                    for (int i = 0; i < staFrames[f].length; i++) {
                        current[n] += staFrames[f][i] * movie[firstFrame + f][i];
                    }
                }
            }
        }
    }


    public static DoubleHistogram[][] getMovingBarResponses(
        int id, double dt, String baseName, int[][] dataSetNummers) throws
        IOException {

        DoubleHistogram[][] d = new DoubleHistogram[dataSetNummers.length][];
        for (int vIndex = 0; vIndex < dataSetNummers.length; vIndex++) {
            String[] sName = new String[dataSetNummers[vIndex].length];
            for (int i = 0; i < sName.length; i++) {
                sName[i] = baseName + File.separator + "stimuli" + File.separator +
                           "s" + StringUtil.format(dataSetNummers[vIndex][i], 0, 2) +
                           ".txt";
//                System.out.println("stimulus: " + sName[i]);
            }
//            int[][] stimulus = VisionUtil.loadMovingFilterStimuli(sName, null);
//            final int maxRunNumber = MathUtil.max(stimulus[0]) + 1;
//            d[vIndex] = new DoubleHistogram[maxRunNumber];
//            int[] totalRepeats = new int[maxRunNumber];

            for (int n = 0; n < dataSetNummers[vIndex].length; n++) {
                String setName =
                    "data" + StringUtil.format(dataSetNummers[vIndex][n], 0, 3);
                String nName = baseName + File.separator + setName +
                               File.separator + setName +
                               VisionParams.NEURON_FILE_EXTENSION;
                NeuronFile nf = new NeuronFile(nName);
//                DoubleHistogram[] h = MovingBar.getMovingFilterResponse(
//                    nf.getSpikeTimes(id), nf.getTTLTimes(), stimulus[n],
//                    dt / /*(vIndex + 1) / */ 1000);
//Math.pow(2, vIndex)
//                for (int runIndex = 0; runIndex < maxRunNumber; runIndex++) {
//                    final int nRepeats = MathUtil.countValues(runIndex, stimulus[n]);
//                    totalRepeats[runIndex] += nRepeats;

//                    PlotUtilities.showData(h[runIndex], new HistogramStyle());

//                    h[runIndex].scale(nRepeats);
//                    if (d[vIndex][runIndex] == null) {
//                        d[vIndex][runIndex] = h[runIndex];
//                    } else {
//                        d[vIndex][runIndex].fill(h[runIndex]);
//                    }
//                }
            }

//            for (int runIndex = 0; runIndex < maxRunNumber; runIndex++) {
//                d[vIndex][runIndex].scale(1.0 / totalRepeats[runIndex]);
//                System.out.println("vIndex " + vIndex + ", run " + runIndex + ", n " +
//                                   totalRepeats[runIndex]);
//            }
        }

        return d;
    }


    public static void main(String[] args) throws IOException {
        int id = 899;
        String baseName = "c:\\data\\2003-09-19-3";
        int[][] LR = { {9, 17}
                     , {11, 19}
                     , {13, 21}
                     , {15, 23}
        };
        int[][] UD = { {10}
                     , {12}
                     , {14}
                     , {16}
        };

        double precisionFactor = 1;
        final double dt = 8.34 / precisionFactor;
        final double f = 8;
        final int w = (int) (640 / f), h = (int) (320 / f);
        final double pSize = 5e-3 * f;

        // MOVING BARS
        double x1 = -3 * 16 / f;
        double x2 = 2000;
        double y1 = h / 2;
        double y2 = h / 2;
        MovingBarMovie m2lr = new MovingBarMovie(
            500, w, h, pSize, dt / 1, 1000, 16 / f, x1, y1, x2, y2,
            2 / f / precisionFactor / 1, 0.5, 0.48);
        MovingBarMovie m4lr = new MovingBarMovie(
            500, w, h, pSize, dt / 1, 1000, 16 / f, x1, y1, x2, y2,
            4 / f / precisionFactor / 1, 0.5, 0.48);
        MovingBarMovie m8lr = new MovingBarMovie(
            500, w, h, pSize, dt / 1, 1000, 16 / f, x1, y1, x2, y2,
            8 / f / precisionFactor / 1, 0.5, 0.48);
        MovingBarMovie m16lr = new MovingBarMovie(
            500, w, h, pSize, dt / 1, 1000, 16 / f, x1, y1, x2, y2,
            16 / f / precisionFactor / 1, 0.5, 0.48);
        MovingBarMovie m2lrB = new MovingBarMovie(
            500, w, h, pSize, dt / 1, 1000, 16 / f, x1, y1, x2, y2,
            2 / f / precisionFactor / 1, 0.5, -0.48);
        MovingBarMovie m4lrB = new MovingBarMovie(
            500, w, h, pSize, dt / 1, 1000, 16 / f, x1, y1, x2, y2,
            4 / f / precisionFactor / 1, 0.5, -0.48);
        MovingBarMovie m8lrB = new MovingBarMovie(
            500, w, h, pSize, dt / 1, 1000, 16 / f, x1, y1, x2, y2,
            8 / f / precisionFactor / 1, 0.5, -0.48);
        MovingBarMovie m16lrB = new MovingBarMovie(
            500, w, h, pSize, dt / 1, 1000, 16 / f, x1, y1, x2, y2,
            16 / f / precisionFactor / 1, 0.5, -0.48);

        x1 = w / 2;
        x2 = w / 2;
        y1 = -3 * 16 / f;
        y2 = 1000;
        MovingBarMovie m2udW = new MovingBarMovie(
            500, w, h, pSize, dt / 1, 1000, 16 / f, x1, y1, x2, y2,
            2 / f / precisionFactor / 1, 0.5, 0.48);
        MovingBarMovie m4udW = new MovingBarMovie(
            500, w, h, pSize, dt / 1, 1000, 16 / f, x1, y1, x2, y2,
            4 / f / precisionFactor / 1, 0.5, 0.48);
        MovingBarMovie m8udW = new MovingBarMovie(
            500, w, h, pSize, dt / 1, 1000, 16 / f, x1, y1, x2, y2,
            8 / f / precisionFactor / 1, 0.5, 0.48);
        MovingBarMovie m16udW = new MovingBarMovie(
            500, w, h, pSize, dt / 1, 1000, 16 / f, x1, y1, x2, y2,
            16 / f / precisionFactor / 1, 0.5, 0.48);
        MovingBarMovie m2udB = new MovingBarMovie(
            500, w, h, pSize, dt / 1, 1000, 16 / f, x1, y1, x2, y2,
            2 / f / precisionFactor / 1, 0.5, -0.48);
        MovingBarMovie m4udB = new MovingBarMovie(
            500, w, h, pSize, dt / 1, 1000, 16 / f, x1, y1, x2, y2,
            4 / f / precisionFactor / 1, 0.5, -0.48);
        MovingBarMovie m8udB = new MovingBarMovie(
            500, w, h, pSize, dt / 1, 1000, 16 / f, x1, y1, x2, y2,
            8 / f / precisionFactor / 1, 0.5, -0.48);
        MovingBarMovie m16udB = new MovingBarMovie(
            500, w, h, pSize, dt / 1, 1000, 16 / f, x1, y1, x2, y2,
            16 / f / precisionFactor / 1, 0.5, -0.48);

        x1 = w / 2;
        x2 = w / 2;
        y1 = h + 3 * 16 / f;
        y2 = -1000;
        MovingBarMovie m2du = new MovingBarMovie(
            500, w, h, pSize, dt / 1, 1000, 16 / f, x1, y1, x2, y2,
            2 / f / precisionFactor / 1, 0.5, 0.48);
        MovingBarMovie m4du = new MovingBarMovie(
            500, w, h, pSize, dt / 1, 1000, 16 / f, x1, y1, x2, y2,
            4 / f / precisionFactor / 1, 0.5, 0.48);
        MovingBarMovie m8du = new MovingBarMovie(
            500, w, h, pSize, dt / 1, 1000, 16 / f, x1, y1, x2, y2,
            8 / f / precisionFactor / 1, 0.5, 0.48);
        MovingBarMovie m16du = new MovingBarMovie(
            500, w, h, pSize, dt / 1, 1000, 16 / f, x1, y1, x2, y2,
            16 / f / precisionFactor / 1, 0.5, 0.48);

        DoubleHistogram[][] barResponseLR = getMovingBarResponses(id, dt, baseName, LR);
        DoubleHistogram[][] barResponseUD = getMovingBarResponses(id, dt, baseName, UD);

        // FULL SCREEN FLASHES
//        FullScreenFlashMovie fm = new FullScreenFlashMovie(
//            w, h, pSize, dt, 1000, 0.48f, 1000, 0f);
        //MG Replaced code, but didn't test.  Check nFlashBlocks in FlashesClassification Constructor before using
//        FlashesClassification flashesClassification = new FlashesClassification(new
//            String[] {"data005"}, "c:\\data\\2003-09-19-3\\", 150);
//        DoubleHistogram hFlash1 = flashesClassification.getFlashResponseHistogram(id);
//        NeuronFile nfFlash1 = new NeuronFile(
//            "c:\\data\\2003-09-19-3\\data005\\data005.neurons");
//        DoubleHistogram hFlash1 = VisionUtilities.getFlashResponse(
//            nfFlash1.getSpikeTimes(id), nfFlash1.getTTLTimes(), dt / 1000.0);
        /*
                // REVERSING GRATINGS
                ReversingGratings[] m = new ReversingGratings[2];
                m[0] = new ReversingGratings("c:\\data\\2003-09-19-3\\data030");
                m[1] = new ReversingGratings("c:\\data\\2003-09-19-3\\data031");
                ReversingGratingMovie[][][] rGratingMovie =
                    new ReversingGratingMovie[m.length][m[0].nPeriods][m[0].nPhases];
         DoubleHistogram[][] revGratingsResponse = new DoubleHistogram[m.length][];
                for (int n = 0; n < m.length; n++) {
                    m[n].setCurrentNeuron(id);
                    revGratingsResponse[n] = m[n].getResponse(dt / 1000, 2);
//            revGratingsResponse[n] = m[n].getAverageReversingGratingsResponse(dt / 1000);
                    for (int period = 0; period < m[n].nPeriods; period++) {
                        for (int i = 0; i < m[n].nPhases; i++) {
                            double phase = i * m[n].periods[period] / 8;
         rGratingMovie[n][period][i] = new ReversingGratingMovie(w, h, pSize,
         dt, m[n].getContrast(), m[n].periods[period] / f, phase / f,
                                m[n].getTemporalPeriod() * 8.34);
                        }
                    }
                }
         */
        // DRIFTING CONTRAST SINUSOIDS
        DriftingContrastSinusoids m1 = new DriftingContrastSinusoids(
            "c:\\data\\2003-09-19-3\\data032");
        DriftingSinusoidMovie[] driftingMovie = new DriftingSinusoidMovie[m1.nContrasts];
        DoubleHistogram[] driftingResponse = new DoubleHistogram[m1.nContrasts];
        m1.setCurrentNeuron(id);
        driftingResponse = m1.getAverageResponse(dt / 1000);
//        driftingResponse = m1.getResponse(dt / 1000, 4);
        for (int contrast = 0; contrast < m1.nContrasts; contrast++) {
            driftingMovie[contrast] = new DriftingSinusoidMovie(
                1000, w, h, pSize, dt, 96 / f, m1.getTemporalPeriod() * 8.34,
                m1.contrasts[contrast], 1000);
        }

        // DRIFTING SINUSOID
//        NeuronFile dsnf1 = new NeuronFile(
//            "c:\\data\\2003-09-19-3\\data026\\data026.neurons");
//        double[] dir1 = VisionUtilities.loadDriftingSinusoidStimulus(
//            "c:\\data\\2003-09-19-3\\stimuli\\data026.txt");
//        DriftingSinusoidMovie sinusoid1 = new DriftingSinusoidMovie(
//            7194, w, h, pSize, dt, 64 / f, 32 * 8.34, 0.48, 3597);
//        DoubleHistogram sinData1 = VisionUtilities.getAverageDriftingSinusoidsResponse(
//            dsnf1.getSpikeTimes(id), dsnf1.getTTLTimes(), dir1, dt / 1000);
//
//        NeuronFile dsnf2 = new NeuronFile(
//            "c:\\data\\2003-09-19-3\\data027\\data027.neurons");
//        double[] dir2 = VisionUtilities.loadDriftingSinusoidStimulus(
//            "c:\\data\\2003-09-19-3\\stimuli\\data027.txt");
//        DriftingSinusoidMovie sinusoid2 = new DriftingSinusoidMovie(
//            7194, w, h, pSize, dt, 64 / f, 32 * 8.34, 0.48, 3597);
//        DoubleHistogram sinData2 = VisionUtilities.getAverageDriftingSinusoidsResponse(
//            dsnf2.getSpikeTimes(id), dsnf2.getTTLTimes(), dir2, dt / 1000);

        NeuronFile dsnf3 = new NeuronFile(
            "c:\\data\\2003-09-19-3\\data028\\data028.neurons");
//        double[] dir3 = VisionUtilities.loadDriftingSinusoidStimulus(
//            "c:\\data\\2003-09-19-3\\stimuli\\data028.txt");
        DriftingSinusoidMovie sinusoid3 = new DriftingSinusoidMovie(
            7194, w, h, pSize, dt, 64 / f, 32 * 8.34, 0.48, 7194);
//        DoubleHistogram sinData3 = VisionUtilities.getAverageDriftingSinusoidsResponse(
//            dsnf3.getSpikeTimes(id), dsnf3.getTTLTimes(), dir3, dt / 1000);

        /*  31
                {0, 1.7816674722399253},
                {1, 0.6395285032432203},
                {0, 0.64 - 0.1147},
                {0, 0.64 - 0.04},
                {0, 16.900794766943584},
                {1, 8.592948519626237},
                {1, 1.7178975051925183E11},
                {1, 142.71335839815202},
                {1, 0.2977389853150494},
                {1, 2.735528972719633},
                {0, 1.27},
                {0, 0.048},
                {0, -37.7895},
                {0, -94.2682},
                {0, -0.30027562094580296},
                {0, 13.0},
                {0, 8.0},
                {0, -51.26},
                {0, -64.51},
                {0, -0.231},
                {0, 15.0},
                {0, 1.0},
         */
        /* 37
                             {0, 1.7816674722399253}
                             , {1, 0.5}
                             , {0, 0.5 - 0.1147}
                             , {0, 0.5 - 0.04}
                             , {0, 16.900794766943584}
                             ,

                             {1, 13.347665048744094}
                             , {1, 0.5434004130641261}
                             , {1, 114.5278971752148}
                             , {1, 0.25122928532113836}
                             , {1, 1.2702031604478015}
                             ,

                             {0, 1.67}
                             , {0, 0.043}
                             , {0, -37.81}
                             , {0, -103.07}
                             , {0, -0.313}
                             , {0, 12}
                             , {0, 10}

                             , {0, -53.47}
                             , {0, -134.22}
                             , {0, -0.256}
                             , {0, 23}
                             , {0, 9}
                             ,
         */

        final double[][] p = { {0, 1.7816674722399253}, {0, 1.225019116982974}, {1, 0.072},
                             {0, 0.04}, {0, 16.900794766943584}, {0, 11.190447956366672},
                             {0, 1.7178975051935184E11}, {0, 125.14368746637771}, {0,
                             0.22594780836306633}, {0, 2.077401556818934}, {0, 2.73}, {0,
                             0.0495}, {0, -39.5202}, {0, -90.443}, {0,
                             -0.24372261456062377}, {0, 14.0}, {0, 5.0}, {0, -51.5616},
                             {0, -116.8916}, {0, -0.15619865430674443}, {0, 30.0}, {0,
                             39.0},
        };

// cgc
//        {0, 0.76283233527005},
//        {0, 41.39186576981667},

//        final SimulationParasolVictor s = new SimulationParasolVictor(new Object[][] {
        // contrast 0.24
//            {rGratingMovie[0][0][0], revGratingsResponse[0][0 * 4 + 0]},
//             {rGratingMovie[0][0][1], revGratingsResponse[0][0 * 4 + 1]},
//             {rGratingMovie[0][0][2], revGratingsResponse[0][0 * 4 + 2]},
//             {rGratingMovie[0][0][3], revGratingsResponse[0][0 * 4 + 3]},
//             {rGratingMovie[0][1][0], revGratingsResponse[0][1 * 4 + 0]},
//             {rGratingMovie[0][1][1], revGratingsResponse[0][1 * 4 + 1]},
//             {rGratingMovie[0][1][2], revGratingsResponse[0][1 * 4 + 2]},
//             {rGratingMovie[0][1][3], revGratingsResponse[0][1 * 4 + 3]},
        // contrast 0.48
//            {rGratingMovie[1][0][0], revGratingsResponse[1][0 * 4 + 0]},
//            {rGratingMovie[1][0][1], revGratingsResponse[1][0 * 4 + 1]},
//            {rGratingMovie[1][0][2], revGratingsResponse[1][0 * 4 + 2]},
//            {rGratingMovie[1][0][3], revGratingsResponse[1][0 * 4 + 3]},
//            {rGratingMovie[1][1][0], revGratingsResponse[1][1 * 4 + 0]},
//            {rGratingMovie[1][1][1], revGratingsResponse[1][1 * 4 + 1]},
//            {rGratingMovie[1][1][2], revGratingsResponse[1][1 * 4 + 2]},
//            {rGratingMovie[1][1][3], revGratingsResponse[1][1 * 4 + 3]},
//            {rGratingMovie[1][2][0], revGratingsResponse[1][2 * 4 + 0]},
//            {rGratingMovie[1][2][1], revGratingsResponse[1][2 * 4 + 1]},
//            {rGratingMovie[1][2][2], revGratingsResponse[1][2 * 4 + 2]},
//            {rGratingMovie[1][2][3], revGratingsResponse[1][2 * 4 + 3]},
//            // DRIFTING CONTRAST SINUSOIDS
//            {driftingMovie[0], driftingResponse[0]},
//            {driftingMovie[1], driftingResponse[1]},
//            {driftingMovie[2], driftingResponse[2]},
//            {driftingMovie[3], driftingResponse[3]},
//            {driftingMovie[4], driftingResponse[4]},

        // DRIFTING SINUSOIDS
//            {sinusoid1, sinData1},
//            {sinusoid2, sinData2},
//            {sinusoid3, sinData3},

        // U -> D bars
//            {m2udW, barResponseUD[0][0]},
//            {m4udW, barResponseUD[1][1]},
//            {m8udW, barResponseUD[2][3]},
//            {m16udW, barResponseUD[3][2]},
//            {m2udB, barResponseUD[0][1]},
//            {m4udB, barResponseUD[1][2]},
//            {m8udB, barResponseUD[2][0]},
//            {m16udB, barResponseUD[3][0]},
//            {fm, hFlash1},

        // D -> U bars
//            {m2du, barResponseUD[0][2]},
//            {m4du, barResponseUD[1][3]},
//            {m8du, barResponseUD[2][1]},
//            {m16du, barResponseUD[3][1]},
//
//            ,

        // L->R bars
//            {m2lr, barResponseLR[0][2]},
//            {m4lr, barResponseLR[1][2]},
//            {m8lr, barResponseLR[2][0]},
//            {m16lr, barResponseLR[3][3]},
//            {m2lrB, barResponseLR[0][0]},
//            {m4lrB, barResponseLR[1][0]},
//            {m8lrB, barResponseLR[2][2]},
//            {m16lrB, barResponseLR[3][0]},
//        }
//            , p);

        new Thread(new Runnable() {
            public void run() {
                double t1 = System.currentTimeMillis();
//                s.simulate();
//                s.DEBUG = true;
                for (int i = 0; i < 1; i++) {
//                    System.out.println("L = " + s.getValue());
                }
                double t2 = System.currentTimeMillis();
                System.out.println("Took: " + (t2 - t1) / 1000 + " seconds.");
            }
        }).start();
    }


    public void readParameters(double[] x) {
        int xIndex = 0;

        x0 = x[xIndex++];
        y0 = x[xIndex++];
        dyReversingGratings = x[xIndex++];
        dyDriftingSinusoinds = x[xIndex++];
        dtReversingGratings = x[xIndex++];
        centerA = x[xIndex++];
        surroundSigmaX = x[xIndex++];

        nlfA = x[xIndex++];
        nlfI0 = x[xIndex++];
        nlfB = x[xIndex++];

//        cgcB = x[xIndex++];
//        cgcT = x[xIndex++];

        double alpha = x[xIndex++];
        centerSigmaX = x[xIndex++];

        centerTCt1 = x[xIndex++];
        centerTCt2 = x[xIndex++];
        centerTCa = x[xIndex++];
        centerTCn1 = x[xIndex++];
        centerTCn2 = x[xIndex++];

        surroundTCt1 = x[xIndex++];
        surroundTCt2 = x[xIndex++];
        surroundTCa = x[xIndex++];
        surroundTCn1 = x[xIndex++];
        surroundTCn2 = x[xIndex++];

        // derived
        surroundSigmaY = surroundSigmaX;
        centerSigmaY = centerSigmaX;
        surroundA = - (centerA / alpha) * Math.pow(centerSigmaX / surroundSigmaX, 2);
    }
}
