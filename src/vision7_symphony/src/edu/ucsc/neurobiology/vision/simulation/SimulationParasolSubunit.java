package edu.ucsc.neurobiology.vision.simulation;

import java.io.*;

import java.awt.*;
import javax.swing.*;

import edu.ucsc.neurobiology.vision.analysis.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.math.fitting.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.stimulus.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class SimulationParasolSubunit
    extends Function implements FunctionMinimization.IterationInfo {

    static boolean DEBUG = false;
    boolean SECOND_FILETRING = false;
    boolean USE_CNTRAST_GAIN = false;

//    boolean CONTRAST_GAIN = true;
    public static final int CENTER_SURROUND = 0;
    public static final int SUBUNIT = 1;
    public static final int COMPLETE = 2;

    JPanel panel;
    AdvancedMovie[] movie;
    double refreshTime, stixelWidth, stixelHeight, stixelArea;
    double pSize;
    int width, height, staDepth;
    int nData;
    Spline[] spline;
    private double[][] currentX;
    double[][] dataON;
//    float[][] subunitSTA;
    int[] dataLength;
    DoubleHistogram[] dataHistON;
    double[] subunitCenterTC, subunitSurroundTC /*, secondFilter*/;
    double[] subunitWeight;
//    float[][][] projectedMovie;

    // 0  1  2   3    4    5
    // 1  7  19  37   61   91
    double[][] currentON, middleCurrent;
    double[][] simulatedDataON;


    // the number of layers in the pool
    static final int subunitPoolRange = 100;
    // this is the sum of an algebraic series
    int nSubunits;
    double[] subunitX, subunitY;

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

//    double dt;

//    double cgcT;
//    double cgcB;
    double weightingS;
    double nlfA, nlfI0, nlfB;
//    double gain;
    double x0;
    double y0, dyRevGratings, dtRevGratings, dyDriftingGratings;
    double subunitCenterA;
    double subunitSurroundA;
    double subunitCenterS;
    double subunitSurroundS;
    double subunitSpacing;

    double subunitCenterTCt1;
    double subunitCenterTCt2;
    double subunitCenterTCa;
    double subunitCenterTCn1;
    double subunitCenterTCn2;
    double subunitSurroundTCt1;
    double subunitSurroundTCt2;
    double subunitSurroundTCa;
    double subunitSurroundTCn1;
    double subunitSurroundTCn2;
    double nonlinearityFraction;

    int centerSubunitIndex;
    STATimeFunction1 subunitCenterTCFunction, subunitSurroundTCFunction;

    JFrame frame = new JFrame();
    PlotPanel[] currentPanel;

    final double staLength = 500; // in ms

    double[] parameters;
    double[] freeParameters;
    boolean[] freeParameter;


    public SimulationParasolSubunit(Object[][] input, double[][] x) throws IOException {
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
        panel = new JPanel(new GridLayout(0, 2));
//        fr.add(new JScrollPane(panel));
        frame.add(panel);

        this.movie = movie;
        this.dataHistON = dON;
//        this.dataHistOFF = dOFF;
        this.nData = movie.length;
        this.dataON = new double[nData][];
//        this.dataOFF = new double[nData][];
        this.currentON = new double[nData][];
//        this.currentOFF = new double[nData][];
        this.middleCurrent = new double[nData][];
//        this.simulatedDataOFF = new double[nData][];
        this.simulatedDataON = new double[nData][];
        this.refreshTime = movie[0].getRefreshTime();
        this.staDepth = (int) Math.ceil(staLength / refreshTime);
        subunitCenterTC = new double[staDepth];
        subunitSurroundTC = new double[staDepth];
//        secondFilter = new double[staDepth];
        width = movie[0].getWidth();
        height = movie[0].getHeight();
        pSize = movie[0].getFrame(0).getStixelWidth();
//        subunitSTA = new float[staDepth][width * 3];

        currentX = new double[nData][];
        spline = new Spline[nData];
        currentPanel = new PlotPanel[nData];
        for (int i = 0; i < nData; i++) {
            currentPanel[i] = new PlotPanel();
        }

        this.stixelWidth = movie[0].getFrame(0).getStixelWidth();
        this.stixelWidth = movie[0].getFrame(0).getStixelHeight();
        this.stixelArea = stixelWidth * stixelHeight;
        dataLength = new int[nData];
//        projectedMovie = new float[nData][][];
        for (int n = 0; n < nData; n++) {
            this.dataON[n] = dON[n].toArray();
            dataLength[n] = dataON[n].length;
            currentX[n] = MathUtil.doubleArray(0, dataLength[n]);
            System.out.println("dataLength[n] " + dataLength[n] + ", " +
                               dON[n].getBinCount());

//            currentOFF[n] = new double[dataLength[n]];
            currentON[n] = new double[dataLength[n]];
            middleCurrent[n] = new double[dataLength[n]];
            simulatedDataON[n] = new double[dataLength[n]];

            // project the movies
//            projectedMovie[n] = new float[dataLength[n]][width * 3];
//            for (int f = 0; f < dataLength[n]; f++) {
//                ImageFrame frame = getFrame(movie[n], f * refreshTime);
//                for (int i = 0; i < width; i++) {
//                    projectedMovie[n][f][3 * i + 0] = frame.getPixel(i, 0, 0) - 0.5f;
//                    projectedMovie[n][f][3 * i + 1] = frame.getPixel(i, 0, 1) - 0.5f;
//                    projectedMovie[n][f][3 * i + 2] = frame.getPixel(i, 0, 2) - 0.5f;
//                }
//            }
        }

        for (int i = 0; i < nData; i++) {
            spline[i] = new Spline(dataLength[i]);
        }
    }


    public static ImageFrame getFrame(Movie movie, final double time) throws IOException {
        return movie.getFrame( (int) Math.floor(time / movie.getRefreshTime()));
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
            1, subunitCenterTCt1, subunitCenterTCa, subunitCenterTCt2,
            subunitCenterTCn1, subunitCenterTCn2);
        subunitSurroundTCFunction = new STATimeFunction1(
            1, subunitSurroundTCt1, subunitSurroundTCa,
            subunitSurroundTCt2, subunitSurroundTCn1, subunitSurroundTCn2);

        /*
//        STATimeFunction1 secondTimeFilter = new STATimeFunction1(
//            TCa1, TCt1, TCa2, TCt2, 12);
                for (int f = 0; f < staDepth; f++) {
                    final double a = - (staDepth - 1 - f + 1) * refreshTime;
                    final double b = - (staDepth - 1 - f + 0) * refreshTime;
                    final double t = - (staDepth - 1 - f + 0.5) * refreshTime;
                        secondFilter[f] = Integrator.qsimp(secondTimeFilter, a, b, 1e-5);
                }
         */
        // layout the subunits
        layoutSubunits();

        // construct the STA of the middle subunit
//        constructSTA(subunitX[centerSubunitIndex], subunitY[centerSubunitIndex],
//                     subunitSTA);

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
                double dx = currentON[n][i] - nlfI0;
                simulatedDataON[n][i] = (dx > 0) ? nlfB + nlfA * Math.pow(dx, 1) : nlfB;
            }
        }

        if (DEBUG) {
            showResult();
        }

        // calculate the likelihood
        double l = 0;

        /*
                 for (int n = 0; n < nData; n++) {
            for (int i = 0; i < dataLength[n]; i++) {
                l += simulatedDataON[n][i] -
                    dataON[n][i] * Math.log(simulatedDataON[n][i]);
                if (dataHistOFF[n] != null) {
                    l += simulatedDataOFF[n][i] -
                        dataOFF[n][i] * Math.log(simulatedDataOFF[n][i]);
                }
            }
                 }
         */

        for (int n = 0; n < nData; n++) {
            for (int i = 0; i < dataLength[n]; i++) {
                l += simulatedDataON[n][i] -
                    dataON[n][i] * Math.log(simulatedDataON[n][i]);
            }
        }

        return l;
    }


    void layoutSubunits() {
        PlotPanel layoutPanel = new PlotPanel();
        layoutPanel.setRange(0, 3.2, 0, 1.6);
        FunctionStyle st = new FunctionStyle("f");

        // layout the subunits
        DoubleList xList = new DoubleList();
        DoubleList yList = new DoubleList();
        for (int nx = -2 * subunitPoolRange; nx <= +2 * subunitPoolRange; nx++) {
            for (int ny = -2 * subunitPoolRange; ny <= +2 * subunitPoolRange; ny++) {
                if (nx * nx - nx * ny + ny * ny <= subunitPoolRange * subunitPoolRange) {
                    double rand = 0.01 * subunitSpacing * (Math.random() - 0.5) * 2;
                    rand = 0;
                    double x = x0 + 0.5 * subunitSpacing * (2 * nx - ny) + rand;
                    double y = y0 + 0.5 * subunitSpacing * ny * Math.sqrt(3) + rand;
                    if (x > 0 && x < 3.2 && y > 0 && y < 1.6) {
                        xList.add(x);
                        yList.add(y);

                        if (x == x0 && y == y0) {
                            centerSubunitIndex = xList.size() - 1;
                        }

                        if (DEBUG) {
                            layoutPanel.addData(new ParametricEllipse(x, y,
                                subunitCenterS, subunitCenterS, 0), st);
                        }
                    }
                }
            }
        }
        nSubunits = xList.size();
        subunitX = xList.toArray();
        subunitY = yList.toArray();
        subunitWeight = new double[nSubunits];

        for (int subunitIndex = 0; subunitIndex < nSubunits; subunitIndex++) {
            double dx = x0 - subunitX[subunitIndex];
            double dy = y0 - subunitY[subunitIndex];
            subunitWeight[subunitIndex] +=
                Math.exp( -0.5 * (dx * dx + dy * dy) / (weightingS * weightingS));
        }

        if (DEBUG) {
            layoutPanel.addData(new ParametricEllipse(x0, y0,
                weightingS, weightingS, 0), new FunctionStyle("", Color.black, 2));
            PlotUtil.showData("", layoutPanel);
            System.out.println("subunits " + nSubunits);
        }

    }


    /*
        void constructSTA(double x, double y, float[][] sta) {
            double a, b, t, xx, yy, subunitCenterG, subunitSurroundG;
            Gaussian2DFunction subunitCenterGaussian = new Gaussian2DFunction(
                subunitCenterA, x, y, subunitCenterS, subunitCenterS, 0, 0);
            Gaussian2DFunction subunitSurroundGaussian = new Gaussian2DFunction(
                subunitSurroundA, x, y, subunitSurroundS, subunitSurroundS, 0, 0);
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
                        sta[f][3 * i + 0] += (float) c;
                        sta[f][3 * i + 1] += (float) c;
                        sta[f][3 * i + 2] += (float) c;
                    }
                }
            }
        }
     */

    double halfRectify(double x) {
        if (x > 0) {
            return x;
        } else {
            return 0;
        }
    }


//    double fullRectify(double x) {
//        return Math.abs(x);
//    }


    void calculateMovingFilterCurrent(int n) throws CannotEvaluateException {
        MovingBarMovie m = (MovingBarMovie) movie[n];
        final double velocity = m.getVelocity() * pSize / refreshTime;
        final double w = m.getThickness() * pSize;
        final double I0 = m.getContrast();
        final double ac = subunitCenterS * subunitCenterS + w * w;
        final double as = subunitSurroundS * subunitSurroundS + w * w;
        final double Ac = 2 * Math.PI * I0 * subunitCenterA *
                          subunitCenterS * subunitCenterS * w / Math.sqrt(ac);
        final double As = 2 * Math.PI * I0 * subunitSurroundA *
                          subunitSurroundS * subunitSurroundS * w / Math.sqrt(as);
        final double y1 = m.getY1() * pSize;
        FunctionDataAdapter f = new FunctionDataAdapter() {
            public double getValueAt(double t) {
                double yb = y1 + velocity * t;
                return
                    Ac * Math.exp( - (y0 - yb) * (y0 - yb) / (2 * ac)) *
                    subunitCenterTCFunction.getValueAt(t - t0)
                    +
                    As * Math.exp( - (y0 - yb) * (y0 - yb) / (2 * as)) *
                    subunitSurroundTCFunction.getValueAt(t - t0);
            }
        };
        for (int i = 0; i < dataLength[n]; i++) {
            t0 = i * refreshTime;
            middleCurrent[n][i] = Integrator.gauss8(f, t0 - 500, t0, 1e-5);
            currentON[n][i] = 0;
        }
        // calculate the current for the middle subunit
        // interpolate the current by a cubic spilne
        spline[n].reSpline(currentX[n], middleCurrent[n]);
        double dt;
        for (int subunitIndex = 0; subunitIndex < nSubunits; subunitIndex++) {
            // resample the splite to time-shift the current
            dt = (y0 - subunitY[subunitIndex]) / (velocity * refreshTime); // in frames
            for (int i = 0; i < dataLength[n]; i++) {
                double c = spline[n].getValueAt(i - dt);

                // the linear part
//                currentON[n][i] +=
//                    (1 - nonlinearityFraction) * c * subunitWeight[subunitIndex];

                // the nonlinear part
                currentON[n][i] +=
                    halfRectify(nonlinearityFraction * c * subunitWeight[subunitIndex]);
            }
        }
        /*
                double[] c = new double[dataLength[n]];
                double[] g = new double[dataLength[n] + 1];
                double[] y = new double[dataLength[n]];
                Arrays.fill(g, 1);
                for (int i = 0; i < dataLength[n]; i++) {
                    t0 = (i + 0) * refreshTime;

//            if (USE_CNTRAST_GAIN) {
//                y[i] = g[i] * currentON[n][i];
//                for (int j = 0; j <= i; j++) {
//                    c[i] += rMath.abs(y[j]) * Math.exp( - (i - j) * refreshTime / cgcT);
//                }
//                c[i] *= (cgcB / cgcT) * refreshTime;
//                g[i + 1] = 1 / (1 + Math.pow(c[i], 4));
//            currentON[n][i] = y[i];
//            }
                }
         */
        if (DEBUG) {
            currentPanel[n].addData(new ScatterPlot(currentON[n]), redThin);
        }

        // add the second level of linear time filtering
        if (SECOND_FILETRING) {
//            currentON[n] = secondFilter(secondFilter, currentON[n], refreshTime);
//            currentOFF[n] = secondFilter(secondFilter, currentOFF[n], refreshTime);
//            currentON[n] = dotProduct(secondFilter, currentON[n]);
//            currentOFF[n] = dotProduct(secondFilter, currentOFF[n]);
        }
    }


    double t0 = 0;
    void calculateDriftingSinusoidCurrent(int n) throws CannotEvaluateException {
        DriftingSinusoidMovie m = ( (DriftingSinusoidMovie) movie[n]);
        final double Ty = m.getSpatialPeriod() * pSize;
        final double T = m.getTemporalPeriod();
        final double wy = 2 * Math.PI / Ty;
        final double w = 2 * Math.PI / T;
        final double I0 = m.getContrast();

        double cs2 = subunitCenterS * subunitCenterS;
        double ss2 = subunitSurroundS * subunitSurroundS;
        final double cA = 2 * subunitCenterA * I0 * Math.PI * cs2 *
                          Math.exp( -cs2 * wy * wy / 2);
        final double sA = 2 * subunitSurroundA * I0 * Math.PI * ss2 *
                          Math.exp( -ss2 * wy * wy / 2);
        FunctionDataAdapter C_COS = new FunctionDataAdapter() {
            public double getValueAt(double t) {
                return Math.cos(w * t) * subunitCenterTCFunction.getValueAt(t - t0);
            }
        };
        FunctionDataAdapter C_SIN = new FunctionDataAdapter() {
            public double getValueAt(double t) {
                return Math.sin(w * t) * subunitCenterTCFunction.getValueAt(t - t0);
            }
        };
        FunctionDataAdapter S_COS = new FunctionDataAdapter() {
            public double getValueAt(double t) {
                return Math.cos(w * t) * subunitSurroundTCFunction.getValueAt(t - t0);
            }
        };
        FunctionDataAdapter S_SIN = new FunctionDataAdapter() {
            public double getValueAt(double t) {
                return Math.sin(w * t) * subunitSurroundTCFunction.getValueAt(t - t0);
            }
        };
        double[] c_sin = new double[dataLength[n]];
        double[] c_cos = new double[dataLength[n]];
        double[] s_sin = new double[dataLength[n]];
        double[] s_cos = new double[dataLength[n]];
        for (int i = 0; i < dataLength[n]; i++) {
            t0 = i * refreshTime;
            double t = Math.min(t0, 30 * 1000);
            c_sin[i] = Integrator.gauss8(C_SIN, t - 500, t, 1e-5);
            c_cos[i] = Integrator.gauss8(C_COS, t - 500, t, 1e-5);
            s_sin[i] = Integrator.gauss8(S_SIN, t - 500, t, 1e-5);
            s_cos[i] = Integrator.gauss8(S_COS, t - 500, t, 1e-5);
        }
        for (int subunitIndex = 0; subunitIndex < nSubunits; subunitIndex++) {
            double sin = Math.sin(wy * (subunitY[subunitIndex] - dyDriftingGratings));
            double cos = Math.cos(wy * (subunitY[subunitIndex] - dyDriftingGratings));
            for (int i = 0; i < dataLength[n]; i++) {
                double c = cA * (sin * c_cos[i] + cos * c_sin[i]) +
                           sA * (sin * s_cos[i] + cos * s_sin[i]);

                // the nonlinear part
                currentON[n][i] +=
                    halfRectify(nonlinearityFraction * c * subunitWeight[subunitIndex]);
            }
        }
        // add the second level of linear time filtering
        if (SECOND_FILETRING) {
//            currentON[n] = dotProduct(secondFilter, currentON[n]);
//            currentON[n] = secondFilter(secondFilter, currentON[n], refreshTime);
        }
        if (DEBUG) {
            currentPanel[n].addData(new ScatterPlot(currentON[n]), blackThin);
        }
    }


    /*
        void calculateDriftingSinusoidCurrent(int n) {
            final double Tx = ( (DriftingSinusoidMovie) movie[n]).getSpatialPeriod() *
                              pSize;
            final double T = ( (DriftingSinusoidMovie) movie[n]).getTemporalPeriod() *
                             refreshTime;
            System.out.println("T = " + T);
            final double wx = 2 * Math.PI / Tx;
            final double w = 2 * Math.PI / T;
            final double I0 = 0.48;
            double cs2 = subunitCenterS * subunitCenterS;
            double ss2 = subunitSurroundS * subunitSurroundS;
            final double cA = 2 * subunitCenterA * I0 * Math.PI * cs2 *
                              Math.exp( -cs2 * wx * wx / 2);
            final double sA = 2 * subunitSurroundA * I0 * Math.PI * ss2 *
                              Math.exp( -ss2 * wx * wx / 2);
            FunctionDataAdapter C_COS = new FunctionDataAdapter() {
                public double getValueAt(double t) {
                    return Math.cos(w * t) * subunitCenterTCFunction.getValueAt(t - t0);
                }
            };
            FunctionDataAdapter C_SIN = new FunctionDataAdapter() {
                public double getValueAt(double t) {
                    return Math.sin(w * t) * subunitCenterTCFunction.getValueAt(t - t0);
                }
            };
            FunctionDataAdapter S_COS = new FunctionDataAdapter() {
                public double getValueAt(double t) {
         return Math.cos(w * t) * subunitSurroundTCFunction.getValueAt(t - t0);
                }
            };
            FunctionDataAdapter S_SIN = new FunctionDataAdapter() {
                public double getValueAt(double t) {
         return Math.sin(w * t) * subunitSurroundTCFunction.getValueAt(t - t0);
                }
            };
            FunctionDataAdapter f = new FunctionDataAdapter() {
                public double getValueAt(double t) {
                    return
                        - (cA * subunitCenterTCFunction.getValueAt(t - t0) +
                           sA * subunitSurroundTCFunction.getValueAt(t - t0)) *
                        Math.sin(w * t - wx * x0);
                }
            };
            double[] y = new double[dataLength[n]];
            double[] g = new double[dataLength[n]];
            g[0] = 1;
            double[] c = new double[dataLength[n]];
            double tau = 5 * 1000;
            double B = 750;
            double sin = Math.sin(wx * x0);
            double cos = Math.cos(wx * x0);
            for (int i = 1; i < dataLength[n] / 2; i++) {
                t0 = i * refreshTime;
                double t = Math.min(t0, 30 * 1000);
                double c_sin = Integrator.gauss8(C_SIN, t - 500, t, 1e-5);
                double c_cos = Integrator.gauss8(C_COS, t - 500, t, 1e-5);
                double s_sin = Integrator.gauss8(S_SIN, t - 500, t, 1e-5);
                double s_cos = Integrator.gauss8(S_COS, t - 500, t, 1e-5);
//            y[i] = g[i - 1] *  Integrator.gauss8(f, t - 500, t, 1e-6);
                y[i] = g[i - 1] *  (cA * (sin * c_cos + cos * c_sin) +
                                       sA * (sin * s_cos + cos * s_sin));
                for (int j = 0; j < i; j++) {
                    c[i] += (B / tau) * Math.abs(y[j]) *
                        Math.exp( - ( (i - j) * refreshTime) / tau) * refreshTime;
                }
                g[i] = 1 / (1 + Math.pow(c[i], 4));
            }
            PlotUtilities.showData("c", new ScatterPlot(c), blackThin).autoscale();
            PlotUtilities.showData("g", new ScatterPlot(g), blackThin).autoscale();
            PlotUtilities.showData("y", new ScatterPlot(y), blackThin).autoscale();
            spline[n].reSpline(currentX[n], y);
            for (int subunitIndex = 0; subunitIndex < nSubunits; subunitIndex++) {
                for (int i = 0; i < dataLength[n] / 2; i++) {
                    double v = spline[n].getValueAt(i + 2 * (Math.random() - 0.5) * 32);
                    currentON[n][i] += ( (v > 0) ? v : nlfK2 * v) *
                        subunitWeight[subunitIndex];
                }
            }
            // add the second level of linear time filtering
            if (SECOND_FILETRING) {
                currentON[n] = dotProduct(secondFilter, currentON[n]);
                currentOFF[n] = dotProduct(secondFilter, currentOFF[n]);
            }
            if (DEBUG) {
                currentPanel[n].addData(new ScatterPlot(currentON[n]), blackThin);
            }
        }
     */
    public static double heaviside(double x) {
        return (x < 0) ? 0 : 1;
    }


    void calculateFullScreenFlashCurrent(int n) throws CannotEvaluateException {
        double Io = 0.48;
        final double C =
            2 * Math.PI * Io * subunitCenterA * subunitCenterS * subunitCenterS;
        final double S =
            2 * Math.PI * Io * subunitSurroundA * subunitSurroundS * subunitSurroundS;
        FunctionDataAdapter f = new FunctionDataAdapter() {
            public double getValueAt(double t) {
                return (heaviside(t) - heaviside(t - 1000)) *
                    (C * subunitCenterTCFunction.getValueAt(t - t0) +
                     S * subunitSurroundTCFunction.getValueAt(t - t0));
            }
        };

//        double[] c = new double[dataLength[n]];
//        double[] g = new double[dataLength[n] + 1];
//        Arrays.fill(g, 1);
        double[] y = new double[dataLength[n]];
        for (int i = 0; i < dataLength[n]; i++) {
            t0 = i * refreshTime;
            y[i] = Integrator.gauss8(f, t0 - staLength, t0, 1e-4);
            currentON[n][i] = 0;
        }

        for (int subunitIndex = 0; subunitIndex < nSubunits; subunitIndex++) {
            // the linear part
            for (int i = 0; i < dataLength[n]; i++) {
                // the linear part
//                currentON[n][i] +=
//                    (1 - nonlinearityFraction) * y[i] * subunitWeight[subunitIndex];

                // the nonlinear part
                currentON[n][i] +=
                    halfRectify(nonlinearityFraction * y[i] * subunitWeight[subunitIndex]);
//                currentON[n][i] +=
//                    nonlinearityFraction * y[i] * subunitWeight[subunitIndex];
            }
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
            final double wy = 2 * Math.PI / Ty;
            final double w = 2 * Math.PI / T;
            final double I0 = m.getContrast();
            double cs2 = subunitCenterS * subunitCenterS;
            double ss2 = subunitSurroundS * subunitSurroundS;
            final double cA = 2 * subunitCenterA * I0 * Math.PI * cs2 *
                              Math.exp( -cs2 * wy * wy / 2);
            final double sA = 2 * subunitSurroundA * I0 * Math.PI * ss2 *
                              Math.exp( -ss2 * wy * wy / 2);

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
                t0 = i * refreshTime - dtRevGratings;
                c_sin[i] = Integrator.gauss8(C_SIN, t0 - 500, t0, 1e-5);
                s_sin[i] = Integrator.gauss8(S_SIN, t0 - 500, t0, 1e-5);
                currentON[n][i] = 0;
            }

            for (int subunitIndex = 0; subunitIndex < nSubunits; subunitIndex++) {
     double sin = Math.sin(wy * (subunitY[subunitIndex] + phase + dyRevGratings));
                for (int i = 0; i < dataLength[n]; i++) {
                    double c = (cA * c_sin[i] + sA * s_sin[i]) * sin;

                    // the linear part
//                currentON[n][i] +=
//                    (1 - nonlinearityFraction) * c * subunitWeight[subunitIndex];

                    // the nonlinear part
                    currentON[n][i] +=
     halfRectify(nonlinearityFraction * c * subunitWeight[subunitIndex]);
                }
            }

            // add the second level of linear time filtering
            if (SECOND_FILETRING) {
//            currentON[n] = secondFilter(secondFilter, currentON[n], refreshTime);
//            currentON[n] = dotProduct(secondFilter, currentON[n]);
            }

            if (DEBUG) {
                currentPanel[n].addData(new ScatterPlot(currentON[n]), blackThin);
//            currentPanel[n].addData(new ScatterPlot(middleCurrent[n]), blackThin);
            }
        }
     */

    void showResult() {
        /*
                 PlotPanel pp = new PlotPanel();
                 pp.addData(new Linear1DFunction(1, 0), new FunctionStyle());
                 pp.addData(new Linear1DFunction(nlfK2, 0), new FunctionStyle());
                 pp.setRange( -2, 2, -2, 2);
                 PlotUtilities.showData("Subunit Rectification", pp);
         */
        System.out.println("C/S = " +
                           (subunitCenterA * subunitCenterS * subunitCenterS) /
                           ( -subunitSurroundA * subunitSurroundS * subunitSurroundS));

        System.out.println(subunitCenterTC.length);

        PlotPanel p = new PlotPanel();
        ScatterPlot st1 = new ScatterPlot();
        ScatterPlot st2 = new ScatterPlot();
        for (int i = 0; i < subunitCenterTC.length; i++) {
            double t = - (subunitCenterTC.length - i - 1) * refreshTime;
            st1.add(t, subunitCenterTC[i]);
            st2.add(t, subunitSurroundTC[i]);
        }
        p.addData(st1, redThick);
        p.addData(st2, redThin);
//        p.addData(new ScatterPlot(secondFilter), style3);
        p.autoscale();
        PlotUtil.showData("Subunit TC", p);

        // show spacial filters
        PlotPanel spacePanel = new PlotPanel();
        spacePanel.addData(
            new DOG1DFunction(subunitCenterA, subunitSurroundA, x0, subunitCenterS,
                              subunitSurroundS),
            new FunctionStyle("", Color.red, 1));
        spacePanel.setRange(0, 3.2, -0.5, 5);
        PlotUtil.showData("Receptive Field", spacePanel);

        int pW = 0, pH = 100;

        // show results
        for (int i = 0; i < nData; i++) {
//            currentPanel[n].addData(new ScatterPlot( currentON[n]), redThick);
//            currentPanel[n].addData(new ScatterPlot( currentOFF[n]), blueThick);
            currentPanel[i].autoscale();
            currentPanel[i].setYRange( -2, 3);
            currentPanel[i].setPreferredSize(new Dimension(pW, pH));
            currentPanel[i].setAxisVisible(false);
            panel.add(currentPanel[i]);

            PlotPanel hPanelON = new PlotPanel();
            hPanelON.addData(dataHistON[i], new HistogramStyle());
            ScatterPlot spON = new ScatterPlot("");
            for (int j = 0; j < dataLength[i]; j++) {
                spON.add(j * refreshTime / 1000, simulatedDataON[i][j]);
            }
            hPanelON.addData(spON, redThin);
            hPanelON.autoscale();
            hPanelON.setYRange(0, 180);
            hPanelON.setPreferredSize(new Dimension(pW, pH));
            hPanelON.setAxisVisible(false);
            panel.add(hPanelON);
        }

        frame.setBounds(0, 0, 1280, 1000);
        frame.setVisible(true);
    }


    /*
        public double[] secondFilter(double[] f, double[] x, double dt) {
            double[] y = new double[x.length];
            double[] c = new double[x.length];
            double[] g = new double[x.length];
            g[0] = 1;
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


    public static void dotProduct1(float[][] staFrames, float[][] movie, double[] current) {
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

//        return current;
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

        double precisionFactor = 2;
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

        DoubleHistogram[][] barResponseLR = SimulationParasolVictor.getMovingBarResponses(
            id, dt, baseName, LR);
        DoubleHistogram[][] barResponseUD = SimulationParasolVictor.getMovingBarResponses(
            id, dt, baseName, UD);

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

        /* 37
          {0, 0.5}
                                      , {0, -0.1147}
                                      , {0, 16.900794766943584}
//                             weightingS = x[xIndex++];
//                             subunitCenterA = x[xIndex++];
//                             subunitSurroundS = x[xIndex++];
//                             subunitSpacing = x[xIndex++];
//                             subunitCenterS = x[xIndex++];

                                      , {1, 0.038}
                                      , {1, 10.1}
                                      , {1, 0.542}
                                      , {1, 0.020}
                                      , {1, 0.010}

                                      , {0, 114.5278971752148}
                                      , {0, 0.25122928532113836}
                                      , {0, 1.2702031604478015}

                                      , {0, 0.0}
                                      , {0, -37.81}
                                      , {0, -103.07}
                                      , {0, -0.313}
                                      , {0, 12.0}
                                      , {0, 10.0}
                                      , {0, -53.47}
                                      , {0, -134.22}
                                      , {0, -0.256}
                                      , {0, 23.0}
                                      , {0, 9.0}
                                      ,
         */
        final double[][] p = { {0, 1.225019116982974}, {0, -0.0852}, {0,
                             16.900794766943584}, {1, 0.1}, {0, 0.038}, {0, 10.1}, {0,
                             0.542}, {0, 0.02}, {0, 114.5278971752148}, {0, 0.225}, {0,
                             1.2702031604478015},

                             {0, 1},

                             {0, 2.73}, {0, -39.5202}, {0, -90.443}, {0,
                             -0.24372261456062377}, {0, 14.0}, {0, 5.0}, {0, -51.5616},
                             {0, -116.8916}, {0, -0.15619865430674443}, {0, 30.0}, {0,
                             39.0},
        };

//        final SimulationParasolSubunit s = new SimulationParasolSubunit(new Object[][] {
        // contrast 0.24
//            {rGratingMovie[0][0][0], revGratingsResponse[0][0 * 4 + 0]},
//            {rGratingMovie[0][0][1], revGratingsResponse[0][0 * 4 + 1]},
//            {rGratingMovie[0][0][2], revGratingsResponse[0][0 * 4 + 2]},
//            {rGratingMovie[0][0][3], revGratingsResponse[0][0 * 4 + 3]},
//           , {rGratingMovie[0][1][0], revGratingsResponse[0][1 * 4 + 0]}
//            , {rGratingMovie[0][1][1], revGratingsResponse[0][1 * 4 + 1]}
//            , {rGratingMovie[0][1][2], revGratingsResponse[0][1 * 4 + 2]}
//            , {rGratingMovie[0][1][3], revGratingsResponse[0][1 * 4 + 3]}
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
//            {rGratingMovie[1][3][0], revGratingsResponse[1][3 * 4 + 0]},
//            {rGratingMovie[1][3][1], revGratingsResponse[1][3 * 4 + 1]},
//            {rGratingMovie[1][3][2], revGratingsResponse[1][3 * 4 + 2]},
//            {rGratingMovie[1][3][3], revGratingsResponse[1][3 * 4 + 3]},
        // DRIFTING CONTRAST SINUSOIDS
//            {driftingMovie[0], driftingResponse[0]},
//            {driftingMovie[1], driftingResponse[1]},
//            {driftingMovie[2], driftingResponse[2]},
//            {driftingMovie[3], driftingResponse[3]},
//            {driftingMovie[4], driftingResponse[4]},
//            // DRIFTING SINUSOIDS
//            {sinusoid1, sinData1},
//            {sinusoid2, sinData2},
//            {sinusoid3, sinData3},
        // U -> D bars
//            {m2udW, barResponseUD[0][0]}, {m4udW, barResponseUD[1][1]}, {m8udW,
//            barResponseUD[2][3]}, {m16udW, barResponseUD[3][2]}, {m2udB,
//            barResponseUD[0][1]}, {m4udB, barResponseUD[1][2]}, {m8udB,
//            barResponseUD[2][0]}, {m16udB,
//            barResponseUD[3][0]}, {fm, hFlash1},

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
//        s.SECOND_FILETRING = true;

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

        y0 = x[xIndex++];
        dyRevGratings = x[xIndex++];
        dtRevGratings = x[xIndex++];
        dyDriftingGratings = x[xIndex++];
        x0 = 1.6;

        weightingS = x[xIndex++];
        subunitCenterA = x[xIndex++];
        subunitSurroundS = x[xIndex++];
        subunitCenterS = x[xIndex++];
        subunitSpacing = 2 * subunitCenterS;

        nlfA = x[xIndex++];
        nlfI0 = x[xIndex++];
        nlfB = x[xIndex++];

        nonlinearityFraction = x[xIndex++];

        // constants
        double alpha = x[xIndex++];

        subunitCenterTCt1 = x[xIndex++];
        subunitCenterTCt2 = x[xIndex++];
        subunitCenterTCa = x[xIndex++];
        subunitCenterTCn1 = x[xIndex++];
        subunitCenterTCn2 = x[xIndex++];

        subunitSurroundTCt1 = x[xIndex++];
        subunitSurroundTCt2 = x[xIndex++];
        subunitSurroundTCa = x[xIndex++];
        subunitSurroundTCn1 = x[xIndex++];
        subunitSurroundTCn2 = x[xIndex++];

        // derived
        subunitSurroundA =
            - (subunitCenterA / alpha) * Math.pow(subunitCenterS / subunitSurroundS, 2);
    }
}
