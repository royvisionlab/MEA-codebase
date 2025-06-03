package edu.ucsc.neurobiology.vision.simulation;

import java.io.*;

import java.awt.*;
import javax.swing.*;

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
public class SimulationLargeONUD
    extends Function implements FunctionMinimization.IterationInfo {

    boolean DEBUG = false;
    boolean SECOND_FILETRING = false;
//    boolean CONTRAST_GAIN = true;
    public static final int CENTER_SURROUND = 0;
    public static final int SUBUNIT = 1;
    public static final int COMPLETE = 2;

    AdvancedMovie[] movie;
    float[][][] projectedMovie;
    double refreshTime, stixelWidth, stixelHeight, stixelArea;
    double pSize;
    int width, height, staDepth;
    int nData;
    Spline[] spline;
    private double[][] currentX;
    double[][] dataON, dataOFF;
    float[][] subunitSTA;
    int[] dataLength;
    DoubleHistogram[] dataHistON, dataHistOFF;
    double[] subunitCenterTC, subunitSurroundTC, secondFilter;
    double[] subunitWeight;

    // 0  1  2   3    4    5
    // 1  7  19  37   61   91
    int cellRange = 10;
    double[][] currentON, currentOFF, middleCurrent;
    double[][] simulatedDataON, simulatedDataOFF;


    // the number of layers in the pool
    final int subunitPoolRange = 100;
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

    double dt;

    double cgcT;
    double cgcB;
    double cellSpacing;
    double amacrineWeightA;
    double amacrineWeightS;
    double weightingS;
    double nlfK2;
    double TCa1;
    double TCt1;
    double TCa2;
    double TCt2;
    double nlfA;
    double nlfB;
    double x0;
    double y0;
    double subunitCenterA;
    double subunitSurroundA;
    double subunitCenterS;
    double subunitSurroundS;
    double subunitSpacing;
    double subunitCenterTCt1;
    double subunitCenterTCt2;
    double subunitSurroundTCt1;
    double subunitSurroundTCt2;
    double subunitCenterTCa1, subunitCenterTCa2, subunitSurroundTCa1,
        subunitSurroundTCa2;

    STATimeFunction1 subunitCenterTCFunction;
    STATimeFunction1 subunitSurroundTCFunction;
    int centerSubunitIndex;

    JFrame fr = new JFrame();
    PlotPanel[] currentPanel;
    PlotPanel layoutPanel = new PlotPanel();
    FunctionStyle st = new FunctionStyle("f");

    final double staLength = 250; // in ms

    double[] parameters;
    double[] freeParameters;
    boolean[] freeParameter;


    public SimulationLargeONUD(Object[][] input, double[][] x) throws IOException {
        movie = new AdvancedMovie[input.length];
        DoubleHistogram[] dON = new DoubleHistogram[input.length];
        DoubleHistogram[] dOFF = new DoubleHistogram[input.length];

        for (int i = 0; i < input.length; i++) {
            movie[i] = (AdvancedMovie) input[i][0];
            dON[i] = (DoubleHistogram) input[i][1];
            dOFF[i] = (DoubleHistogram) input[i][2];
        }

        initialize(movie, dON, dOFF);

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


    private void initialize(
        AdvancedMovie[] movie, DoubleHistogram[] dON, DoubleHistogram[] dOFF) throws
        IOException {

        fr.getContentPane().setLayout(new GridLayout(3, 0));
        layoutPanel.setRange(0, 3.2, 0, 1.6);

        this.movie = movie;
        this.dataHistON = dON;
        this.dataHistOFF = dOFF;
        this.nData = movie.length;
        this.dataON = new double[nData][];
        this.dataOFF = new double[nData][];
        this.currentON = new double[nData][];
        this.currentOFF = new double[nData][];
        this.middleCurrent = new double[nData][];
        this.simulatedDataOFF = new double[nData][];
        this.simulatedDataON = new double[nData][];
        this.refreshTime = movie[0].getRefreshTime();
        this.staDepth = (int) Math.ceil(staLength / refreshTime);
        subunitCenterTC = new double[staDepth];
        subunitSurroundTC = new double[staDepth];
        secondFilter = new double[staDepth];
        width = movie[0].getWidth();
        height = movie[0].getHeight();
        pSize = movie[0].getFrame(0).getStixelWidth();
        subunitSTA = new float[staDepth][width * 3];

        currentX = new double[nData][];
        spline = new Spline[nData];
        currentPanel = new PlotPanel[nData];
        for (int i = 0; i < nData; i++) {
            currentPanel[i] = new PlotPanel();
        }

        this.stixelWidth = movie[0].getFrame(0).getStixelWidth();
        this.stixelHeight = movie[0].getFrame(0).getStixelHeight();
        this.stixelArea = stixelWidth * stixelHeight;
        dataLength = new int[nData];
        projectedMovie = new float[nData][][];
        for (int n = 0; n < nData; n++) {
            this.dataON[n] = dON[n].toArray();
            if (dOFF[n] != null) {
                this.dataOFF[n] = dOFF[n].toArray();
//                if (dataON[n].length != dataOFF[n].length) {
//                    throw new IllegalArgumentException(
//                        "dataON[n].length != dataOFF[n].length");
//                }
            }
            dataLength[n] = dataON[n].length;
            currentX[n] = MathUtil.doubleArray(0, dataLength[n]);
            System.out.println("dataLength[n] " + dataLength[n] + ", " +
                               dON[n].getBinCount());

            currentOFF[n] = new double[dataLength[n]];
            currentON[n] = new double[dataLength[n]];
            middleCurrent[n] = new double[dataLength[n]];
            simulatedDataON[n] = new double[dataLength[n]];
            simulatedDataOFF[n] = new double[dataLength[n]];

            // project the movies
            projectedMovie[n] = new float[dataLength[n]][width * 3];
            for (int f = 0; f < dataLength[n]; f++) {
                ImageFrame frame = getFrame(movie[n], f * refreshTime);
                for (int i = 0; i < width; i++) {
                    projectedMovie[n][f][3 * i + 0] = frame.getPixel(i, 0, 0) - 0.5f;
                    projectedMovie[n][f][3 * i + 1] = frame.getPixel(i, 0, 1) - 0.5f;
                    projectedMovie[n][f][3 * i + 2] = frame.getPixel(i, 0, 2) - 0.5f;
                }
            }
        }

        for (int i = 0; i < nData; i++) {
            spline[i] = new Spline(dataLength[i]);
        }
    }


    public static ImageFrame getFrame(Movie movie, final double time) throws IOException {
        return movie.getFrame( (int) Math.floor(time / movie.getRefreshTime()));
    }


//    subunitSurroundA =
//        -subunitCenterA * Math.pow(subunitCenterS / subunitSurroundS, 2);

    public void readParameters(double[] x) {
        subunitCenterA = 1;
        subunitCenterTCa1 = 1;
        subunitSurroundTCa1 = 1;
        TCa1 = 1;

        int xIndex = 0;

        cgcT = x[xIndex++];
        cgcB = x[xIndex++];
        subunitSpacing = x[xIndex++];
        weightingS = x[xIndex++];
        cellSpacing = x[xIndex++];
        amacrineWeightA = x[xIndex++];
        amacrineWeightS = x[xIndex++];
        nlfK2 = x[xIndex++];
        subunitSurroundA = x[xIndex++];
        subunitCenterS = x[xIndex++];
        subunitSurroundS = x[xIndex++];

        x0 = 1.6;
        y0 = x[xIndex++];
        subunitCenterTCt1 = x[xIndex++];
        subunitCenterTCa2 = x[xIndex++];
        subunitCenterTCt2 = x[xIndex++];
        subunitSurroundTCt1 = x[xIndex++];
        subunitSurroundTCa2 = x[xIndex++];
        subunitSurroundTCt2 = x[xIndex++];
        TCt1 = x[xIndex++];
        TCa2 = x[xIndex++];
        TCt2 = x[xIndex++];

        nlfA = x[xIndex++];
        nlfB = x[xIndex++];

        dt = x[xIndex++];
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
        subunitCenterTCFunction = new STATimeFunction1(subunitCenterTCa1,
            subunitCenterTCt1, subunitCenterTCa2, subunitCenterTCt2, 12, 12);
        subunitSurroundTCFunction = new STATimeFunction1(subunitSurroundTCa1,
            subunitSurroundTCt1, subunitSurroundTCa2, subunitSurroundTCt2, 12, 12);

        STATimeFunction1 secondTimeFilter = new STATimeFunction1(
            TCa1, TCt1, TCa2, TCt2, 12, 12);
        for (int f = 0; f < staDepth; f++) {
            final double a = - (staDepth - 1 - f + 1) * refreshTime;
            final double b = - (staDepth - 1 - f + 0) * refreshTime;
            final double t = - (staDepth - 1 - f + 0.5) * refreshTime;
            secondFilter[f] = Integrator.qsimp(secondTimeFilter, a, b, 1e-5);
        }

        // layout the subunits
        layoutSubunits();

        // construct the STA of the middle subunit
        constructSTA(subunitX[centerSubunitIndex], subunitY[centerSubunitIndex],
                     subunitSTA);

        // calculate the currents
        for (int n = 0; n < nData; n++) {
            // clear the currents
            for (int i = 0; i < dataLength[n]; i++) {
                currentON[n][i] = 0;
                currentOFF[n][i] = 0;
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
                simulatedDataON[n][i] = (currentON[n][i] > 0) ?
                                        nlfA * currentON[n][i] + nlfB : nlfB;
                simulatedDataOFF[n][i] = (currentOFF[n][i] > 0) ?
                                         nlfA * currentOFF[n][i] + nlfB : nlfB;
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
                if (dataHistOFF[n] != null) {
                    l += simulatedDataOFF[n][i] -
                        dataOFF[n][i] * Math.log(simulatedDataOFF[n][i]);
                }
            }
        }

        return l;
    }


    void layoutSubunits() {
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

        xList.clear();
        yList.clear();
        DoubleList cellWeight = new DoubleList();
        for (int nx = -2 * cellRange; nx <= +2 * cellRange; nx++) {
            for (int ny = -2 * cellRange; ny <= +2 * cellRange; ny++) {
                if (nx * nx - nx * ny + ny * ny <= cellRange * cellRange) {
                    double x = x0 + 0.5 * cellSpacing * (2 * nx - ny);
                    double y = y0 + 0.5 * cellSpacing * ny * Math.sqrt(3);
                    if (x > 0 - 0.3 && x < 3.2 + 0.3 && y > 0 - 0.3 && y < 1.6 + 0.3) {
                        xList.add(x);
                        yList.add(y);
                        double weight = amacrineWeightA * Math.exp(
                            -0.5 * ( (x - x0) * (x - x0) + (y - y0) * (y - y0)) /
                            (amacrineWeightS * amacrineWeightS));
                        if (x == x0 && y == y0) {
                            weight = 1;
                        }
                        cellWeight.add(weight);

                        if (DEBUG) {
                            layoutPanel.addData(new ParametricEllipse(x, y,
                                weightingS,
                                weightingS, 0), new FunctionStyle("", Color.black, 2));
                        }
                    }
                }
            }
        }
        int nCells = cellWeight.size();

        for (int subunitIndex = 0; subunitIndex < nSubunits; subunitIndex++) {
            subunitWeight[subunitIndex] = 0;
            for (int cellIndex = 0; cellIndex < nCells; cellIndex++) {
                double dx = xList.get(cellIndex) - subunitX[subunitIndex];
                double dy = yList.get(cellIndex) - subunitY[subunitIndex];
                subunitWeight[subunitIndex] += cellWeight.get(cellIndex) *
                    Math.exp( -0.5 * (dx * dx + dy * dy) / (weightingS * weightingS));
            }
        }

        if (DEBUG) {
            layoutPanel.addData(new ParametricEllipse(x0, y0,
                amacrineWeightS, amacrineWeightS, 0),
                                new FunctionStyle("", Color.black, 2));
            PlotUtil.showData("", layoutPanel);
            System.out.println("subunits " + nSubunits);
        }

    }


    void constructSTA(double x, double y, float[][] sta) throws CannotEvaluateException {
        double a, b, t, xx, yy, subunitCenterG, subunitSurroundG;

        Gaussian2DFunction subunitCenterGaussian = new Gaussian2DFunction(
            subunitCenterA, x, y, subunitCenterS, subunitCenterS, 0, 0);
        Gaussian2DFunction subunitSurroundGaussian = new Gaussian2DFunction(
            subunitSurroundA, x, y, subunitSurroundS, subunitSurroundS, 0, 0);

        for (int f = 0; f < staDepth; f++) {
            a = - (staDepth - 1 - f + 1) * refreshTime;
            b = - (staDepth - 1 - f + 0) * refreshTime;
            t = - (staDepth - 1 - f + 0.5) * refreshTime;
            subunitCenterTC[f] = Integrator.qsimp(subunitCenterTCFunction, a, b, 1e-5);
            subunitSurroundTC[f] = Integrator.qsimp(subunitSurroundTCFunction, a, b, 1e-5);

            for (int i = 0; i < width; i++) {
                sta[f][3 * i + 0] = 0;
                sta[f][3 * i + 1] = 0;
                sta[f][3 * i + 2] = 0;

                for (int j = 0; j < height; j++) {
                    xx = (i + 0.5) * stixelWidth;
                    yy = (j + 0.5) * stixelHeight;
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


    void calculateMovingFilterCurrent(int n) throws CannotEvaluateException {
        MovingBarMovie m = ( (MovingBarMovie) movie[n]);
        final double velocity = m.getVelocity() * pSize / refreshTime;
        final double w = m.getThickness() * pSize;
        final double I0 = 0.48;
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
                return 3 * (
                    Ac * Math.exp( - (y0 - yb) * (y0 - yb) / (2 * ac)) *
                    subunitCenterTCFunction.getValueAt(t - t0)
                    +
                    As * Math.exp( - (y0 - yb) * (y0 - yb) / (2 * as)) *
                    subunitSurroundTCFunction.getValueAt(t - t0));
            }
        };

        for (int i = 0; i < dataLength[n]; i++) {
            t0 = i * refreshTime;
            middleCurrent[n][i] = Integrator.gauss8(f, t0 - 500, t0, 1e-5);
        }

        // calculate the current for the middle subunit
        // interpolate the current by a cubic spilne
        spline[n].reSpline(currentX[n], middleCurrent[n]);

        double dt, s;
        for (int subunitIndex = 0; subunitIndex < nSubunits; subunitIndex++) {
            // resample the splite to time-shift the current
            dt = (y0 - subunitY[subunitIndex]) / (velocity * refreshTime); // in frames
            for (int i = 0; i < dataLength[n]; i++) {
                s = spline[n].getValueAt(i - dt);
                currentON[n][i] += ( (s > 0) ? s : nlfK2 * s) *
                    subunitWeight[subunitIndex];
                currentOFF[n][i] += ( (s < 0) ? -s : -nlfK2 * s) *
                    subunitWeight[subunitIndex];
            }
        }

        if (DEBUG) {
            currentPanel[n].addData(new ScatterPlot(currentON[n]), redThin);
            currentPanel[n].addData(new ScatterPlot(currentOFF[n]), blueThin);
        }

        // add the second level of linear time filtering
        if (SECOND_FILETRING) {
//            currentON[n] = secondFilter(secondFilter, currentON[n], refreshTime);
//            currentOFF[n] = secondFilter(secondFilter, currentOFF[n], refreshTime);
            currentON[n] = dotProduct(secondFilter, currentON[n]);
            currentOFF[n] = dotProduct(secondFilter, currentOFF[n]);
        }
    }


    double t0 = 0;
    void calculateDriftingSinusoidCurrent(int n) throws CannotEvaluateException {
        final double Tx = ( (DriftingSinusoidMovie) movie[n]).getSpatialPeriod() *
                          pSize;
        final double T = ( (DriftingSinusoidMovie) movie[n]).getTemporalPeriod() *
                         refreshTime;
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
            double sin = Math.sin(wx * subunitX[subunitIndex]);
            double cos = Math.cos(wx * subunitX[subunitIndex]);
            for (int i = 0; i < dataLength[n]; i++) {
                double c = 3 * cA * (sin * c_cos[i] + cos * c_sin[i]) +
                           3 * sA * (sin * s_cos[i] + cos * s_sin[i]);
                currentON[n][i] += ( (c > 0) ? c : nlfK2 * c) *
                    subunitWeight[subunitIndex];
            }
        }

        // add the second level of linear time filtering
        if (SECOND_FILETRING) {
//            currentON[n] = dotProduct(secondFilter, currentON[n]);
            currentON[n] = secondFilter(secondFilter, currentON[n], refreshTime);
        }

        if (DEBUG) {
            currentPanel[n].addData(new ScatterPlot(currentON[n]), blackThin);
        }
    }


    /*
        double t0 = 0;
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
//            y[i] = g[i - 1] * 3 * Integrator.gauss8(f, t - 500, t, 1e-6);
                y[i] = g[i - 1] * 3 * (cA * (sin * c_cos + cos * c_sin) +
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

    void calculateFullScreenFlashCurrent(int n) {
        // calculate the current for the middle subunit
        dotProduct1(subunitSTA, projectedMovie[n], middleCurrent[n]);
        if (DEBUG) {
//            currentPanel[n].addData(new ScatterPlot( middleCurrent), blackThin);
        }

        // apply the rectifying nonlinearuty
        for (int subunitIndex = 0; subunitIndex < nSubunits; subunitIndex++) {
            for (int i = 0; i < dataLength[n]; i++) {
                double c = (middleCurrent[n][i] > 0) ?
                           middleCurrent[n][i] : nlfK2 * middleCurrent[n][i];
                currentON[n][i] += c * subunitWeight[subunitIndex];
            }
        }

        if (DEBUG) {
            currentPanel[n].addData(new ScatterPlot(currentON[n]), redThin);
        }

        // add the second level of linear time filtering
        if (SECOND_FILETRING) {
//            currentON[n] = secondFilter(secondFilter, currentON[n], refreshTime);
            currentON[n] = dotProduct(secondFilter, currentON[n]);
        }
    }


    /*
        void calculateReversingGratingsCurrent(int n) throws CannotEvaluateException {
     final double Ty = ( (ReversingGratingMovie) movie[n]).getSpatialPeriod() * pSize;
            final double phase = ( (ReversingGratingMovie) movie[n]).getPhase() * pSize;
            final double T = ( (ReversingGratingMovie) movie[n]).getTemporalPeriod();
            final double wy = 2 * Math.PI / Ty;
            final double w = 2 * Math.PI / T;
            final double I0 = 0.48;
            double cs2 = subunitCenterS * subunitCenterS;
            double ss2 = subunitSurroundS * subunitSurroundS;
            final double cA = 2 * subunitCenterA * I0 * Math.PI * cs2 *
                              Math.exp( -cs2 * wy * wy / 2);
            final double sA = 2 * subunitSurroundA * I0 * Math.PI * ss2 *
                              Math.exp( -ss2 * wy * wy / 2);

            FunctionDataAdapter C_SIN = new FunctionDataAdapter() {
                public double getValueAt(double t) {
                    t -= dt;
//                double f = ( (t % 8000) <= 6000) ? 1 : 0;
                    return Math.sin(w * t) * subunitCenterTCFunction.getValueAt(t - t0);
                }
            };

            FunctionDataAdapter S_SIN = new FunctionDataAdapter() {
                public double getValueAt(double t) {
                    t -= dt;
//                double f = ( (t % 8000) <= 6000) ? 1 : 0;
     return Math.sin(w * t) * subunitSurroundTCFunction.getValueAt(t - t0);
                }
            };

            double[] c_sin = new double[dataLength[n]];
            double[] s_sin = new double[dataLength[n]];
            for (int i = 0; i < dataLength[n]; i++) {
                t0 = i * refreshTime - dt;
                c_sin[i] = Integrator.gauss8(C_SIN, t0 - 500, t0, 1e-5);
                s_sin[i] = Integrator.gauss8(S_SIN, t0 - 500, t0, 1e-5);
            }

            for (int subunitIndex = 0; subunitIndex < nSubunits; subunitIndex++) {
                double sin = Math.sin(wy * (subunitY[subunitIndex] + phase));
                for (int i = 0; i < dataLength[n]; i++) {
                    double c = 3 * (cA * c_sin[i] + sA * s_sin[i]) * sin;
                    currentON[n][i] += ( (c > 0) ? c : nlfK2 * c) *
                        subunitWeight[subunitIndex];
                }
            }

///////////////
//        double sin = Math.sin(wx * (x0 + phase));
//        for (int i = 0; i < dataLength[n]; i++) {
//            double c = 3 * (cA * c_sin[i] + sA * s_sin[i]) * sin;
//            middleCurrent[n][i] = ( (c > 0) ? c : nlfK2 * c);
//        }
////////////////

            // add the second level of linear time filtering
            if (SECOND_FILETRING) {
//            currentON[n] = secondFilter(secondFilter, currentON[n], refreshTime);
                currentON[n] = dotProduct(secondFilter, currentON[n]);
            }

            if (DEBUG) {
                currentPanel[n].addData(new ScatterPlot(currentON[n]), blackThin);
//            currentPanel[n].addData(new ScatterPlot(middleCurrent[n]), blackThin);
            }
        }
     */

    void showResult() {
        PlotPanel pp = new PlotPanel();
        pp.addData(new Linear1DFunction(1, 0), new FunctionStyle("f"));
        pp.addData(new Linear1DFunction(nlfK2, 0), new FunctionStyle("f"));
        pp.setRange( -2, 2, -2, 2);
        PlotUtil.showData("Subunit Rectification", pp);

        System.out.println("C/S = " +
                           (subunitCenterA * subunitCenterS * subunitCenterS) /
                           ( -subunitSurroundA * subunitSurroundS * subunitSurroundS));

        // show timecourses
        PlotPanel p = new PlotPanel();
        p.addData(new ScatterPlot(subunitCenterTC), redThick);
        p.addData(new ScatterPlot(subunitSurroundTC), redThin);
        p.addData(new ScatterPlot(secondFilter), style3);
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

        // show results
        for (int n = 0; n < nData; n++) {
//            currentPanel[n].addData(new ScatterPlot( currentON[n]), redThick);
//            currentPanel[n].addData(new ScatterPlot( currentOFF[n]), blueThick);
            currentPanel[n].autoscale();
            fr.add(currentPanel[n]);
        }

        for (int i = 0; i < nData; i++) {
            PlotPanel hPanelON = new PlotPanel();
            hPanelON.addData(dataHistON[i], new HistogramStyle());
            ScatterPlot spON = new ScatterPlot("");
            for (int j = 0; j < dataLength[i]; j++) {
                spON.add(j * refreshTime / 1000, simulatedDataON[i][j]);
            }
            hPanelON.addData(spON, blackThin);
            hPanelON.autoscale();
            fr.add(hPanelON);
        }

        for (int i = 0; i < nData; i++) {
            if (dataHistOFF[i] != null) {
                PlotPanel hPanelOFF = new PlotPanel();
                hPanelOFF.addData(dataHistOFF[i], new HistogramStyle());
                ScatterPlot spOFF = new ScatterPlot("");
                for (int j = 0; j < dataLength[i]; j++) {
                    spOFF.add(j * refreshTime / 1000, simulatedDataOFF[i][j]);
                }
                hPanelOFF.addData(spOFF, blackThin);
                hPanelOFF.autoscale();
                fr.add(hPanelOFF);
            } else {
                fr.add(new JLabel("Empty"));
            }
        }
        fr.setBounds(0, 0, 1280, 1000);
        fr.setVisible(true);
    }


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

        PlotUtil.showData("c", new ScatterPlot(c), blackThin).autoscale();
        PlotUtil.showData("g", new ScatterPlot(g), blackThin).autoscale();

        return y;
    }


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

//        return current;
    }


    public static void main(String[] args) throws IOException {
        int id = 1625; //on

        String[][] nNames = { {
                            "c:\\data\\2003-09-19-3\\data010\\data010.neurons"}
                            , {"c:\\data\\2003-09-19-3\\data012\\data012.neurons"}
                            , {"c:\\data\\2003-09-19-3\\data014\\data014.neurons"}
                            , {"c:\\data\\2003-09-19-3\\data016\\data016.neurons"}
        };

        String[][] sNames = { {"c:\\data\\2003-09-19-3\\stimuli\\s10.txt"}
                            , {"c:\\data\\2003-09-19-3\\stimuli\\s12.txt"}
                            , {"c:\\data\\2003-09-19-3\\stimuli\\s14.txt"}
                            , {"c:\\data\\2003-09-19-3\\stimuli\\s16.txt"}
        };

        DoubleHistogram[][] d = new DoubleHistogram[nNames.length][4];
        double precisionFactor = 1;
        final double dt = 8.34 / precisionFactor;
        final double f = 8;
        final int w = (int) (640 / f), h = (int) (320 / f);
        final double pSize = 5e-3 * f;

        for (int vIndex = 0; vIndex < nNames.length; vIndex++) {
//            int[][] stimulus = VisionUtil.loadMovingFilterStimuli(
//                sNames[vIndex], null);

            int[] totalRepeats = new int[4];

            for (int n = 0; n < 1; n++) {
                NeuronFile nf = new NeuronFile(nNames[vIndex][n]);
//                DoubleHistogram[] _h = MovingBar.getMovingFilterResponse(
//                    nf.getSpikeTimes(id), nf.getTTLTimes(), stimulus[n], dt / 1000);

                for (int runIndex = 0; runIndex < 4; runIndex++) {
//                    int nRepeats = MathUtil.countValues(runIndex, stimulus[n]);
//                    totalRepeats[runIndex] += nRepeats;
//                    _h[runIndex].scale(nRepeats);
//                    if (d[vIndex][runIndex] == null) {
//                        d[vIndex][runIndex] = _h[runIndex];
//                    } else {
//                        d[vIndex][runIndex].fill(_h[runIndex]);
//                    }
                }
            }

            for (int runIndex = 0; runIndex < 4; runIndex++) {
                d[vIndex][runIndex].scale(1.0 / totalRepeats[runIndex]);
            }
        }
        /*
                double x1 = w / 2, x2 = w / 2;
                double y1 = - 3 * 16 / f, y2 = 1000;
                MovingBarMovie m2ud = new MovingBarMovie(
                    500, w, h, pSize, dt, 1000, 16 / f, x1, y1, x2, y2,
                    2 / f / precisionFactor, 0.5, 0.48);
                MovingBarMovie m4ud = new MovingBarMovie(
                    500, w, h, pSize, dt, 1000, 16 / f, x1, y1, x2, y2,
                    4 / f / precisionFactor, 0.5, 0.48);
                MovingBarMovie m8ud = new MovingBarMovie(
                    500, w, h, pSize, dt, 1000, 16 / f, x1, y1, x2, y2,
                    8 / f / precisionFactor, 0.5, 0.48);
                MovingBarMovie m16ud = new MovingBarMovie(
                    500, w, h, pSize, dt, 1000, 16 / f, x1, y1, x2, y2,
                    16 / f / precisionFactor, 0.5, 0.48);
         */
        /*
                FullScreenFlashMovie fm = new FullScreenFlashMovie(
                    w, h, pSize, dt, 1000, 0.48f, 1000, 0f);
                NeuronFile nfFlash1 = new NeuronFile(
                    "c:\\data\\2003-09-19-3\\data005\\data005.neurons");
                DoubleHistogram hFlash1 = VisionUtilities.getFlashResponse(
                    nfFlash1.getSpikeTimes(id), nfFlash1.getTTLTimes(), dt / 1000.0);
         */
        /*
                // DRIFTING SINUSOID 64 spatial, 32 temporal
                NeuronFile dsnf = new NeuronFile(
             "e:\\data\\2003-09-19-3\\data028\\data028-mapped-data001.neurons");
                double[] dir = VisionUtilities.loadDriftingSinusoidStimulus(
                    "e:\\data\\2003-09-19-3\\stimuli\\s28.txt");
                DriftingSinusoidMovie sinusoid = new DriftingSinusoidMovie(
                    7194, w, h, pSize, dt, 64 / f, 32, 3597);
         DoubleHistogram sinData1 = VisionUtilities.getAverageDriftingSinusoidsResponse(
                    dsnf.getSpikeTimes(id), dsnf.getTTLTimes(), dir, dt / 1000);
                DoubleHistogram sinData = new DoubleHistogram("", 0, 30, dt / 1000);
                for (int i = 0; i < sinData.getBinCount(); i++) {
                    sinData.setBin(i, sinData1.getBin(i));
                }
         */
        /*
                // REVERSING GRATINGS
                double[] periods = {96, 48, 24, 12, 9, 6};
                ReversingGratingMovie[][] rGratingMovie =
                    new ReversingGratingMovie[periods.length][4];
                for (int period = 0; period < periods.length; period++) {
                    for (int i = 0; i < rGratingMovie[period].length; i++) {
                        double phase = i * periods[period] / 8;
                        System.out.println("phase = " + phase);
             rGratingMovie[period][i] = new ReversingGratingMovie(w, h, pSize, dt,
                            0.48, periods[period] / f, phase / f, 30 * 8.34);
                    }
                }
                ReversingGratings m = new ReversingGratings(
                    "c:\\data\\2003-09-19-3\\stimuli\\s31.txt");
                NeuronFile nf = new NeuronFile(
                    "c:\\data\\2003-09-19-3\\data031\\data031.neurons");
                DoubleHistogram[] hh = m.getResponse(
                    nf.getSpikeTimes(id), nf.getTTLTimes(), dt / 1000, 4, 1);
         */
// bars
//        {0, 3479.7484620781515},
//        {0, 1.1352764379251896},
//        {1, 0.036671687712853764},
//        {1, 0.11159714258872196},
//        {1, 0.31811566667529106},
//        {1, -0.040142814357977466},
//        {1, 0.482399454726},
//        {1, 0.3399018158693463},
//        {1, -0.07316234229856555},
//        {1, 0.0206756929079698},
//        {1, 0.08754312631221173},
//        {1, 0.762206808711252 - 0.24},
//        {1, -53.82362900570574},
//        {1, -0.5247323618410539},
//        {1, -96.06415930795963},
//        {1, -67.74483112022595},
//        {1, -0.7517586866933716},
//        {1, -85.71694150031941},
//        {1, -13.212785615200337},
//        {1, -0.03621633905323732},
//        {1, -37.10663199238221},
//        {1, 9.687263201779237},
//        {1, 0.10910936611024602},
//        {0, 0.0},

// gratings
//        {0, 3479.7484620781515},
//        {0, 1.1352764379251896},
//        {0, 0.036671687712853764},
//        {0, 0.11159714258872196},
//        {0, 0.31811566667529106},
//        {0, -0.040142814357977466},
//        {0, 0.482399454726},
//        {1, -0.07013220020890618},
//        {0, -0.07316234229856555},
//        {0, 0.0206756929079698},
//        {0, 0.08754312631221173},
//        {1, 0.4852310888468043},
//        {1, -61.232857438464855},
//        {1, -0.5108049413991363},
//        {1, -85.66819415118891},
//        {1, -70.23337341023786},
//        {1, -0.6880940838704132},
//        {1, -89.71502507638093},
//        {1, -4.691218092807134},
//        {1, -0.22144875761191068},
//        {1, -20.16912340700043},
//        {1, 69.35889444006888},
//        {0, 0.10910936611024602},
//        {1, 16.834854341689063},

        final double[][] p = { {0, 3479.7484620781515}
                             , {0, 1.1352764379251896}
                             , {1, 0.036671687712853764}
                             , {1, 0.11159714258872196}
                             , {1, 0.31811566667529106}
                             , {1, -0.040142814357977466}
                             , {1, 0.482399454726}
                             , {1, 0.3399018158693463}
                             , {1, -0.07316234229856555}
                             , {1, 0.0206756929079698}
                             , {1, 0.08754312631221173}
                             , {1, 0.762206808711252 - 0.24}
                             , {1, -53.82362900570574}
                             , {1, -0.5247323618410539}
                             , {1, -96.06415930795963}
                             , {1, -67.74483112022595}
                             , {1, -0.7517586866933716}
                             , {1, -85.71694150031941}
                             , {1, -13.212785615200337}
                             , {1, -0.03621633905323732}
                             , {1, -37.10663199238221}
                             , {1, 9.687263201779237}
                             , {1, 0.10910936611024602}
                             , {0, 0.0}
                             ,
        };

        final SimulationLargeONUD s = new SimulationLargeONUD(new Object[][] {
//           {rGratingMovie[0][0], hh[0 * 4 + 0], null},
//           {rGratingMovie[0][1], hh[0 * 4 + 1], null},
//           {rGratingMovie[0][2], hh[0 * 4 + 2], null},
//           {rGratingMovie[0][3], hh[0 * 4 + 3], null},
//             {rGratingMovie[1][0], hh[1 * 4 + 0], null},
//             {rGratingMovie[1][1], hh[1 * 4 + 1], null},
//             {rGratingMovie[1][2], hh[1 * 4 + 2], null},
//             {rGratingMovie[1][3], hh[1 * 4 + 3], null},
//            {rGratingMovie[2][0], hh[2 * 4 + 0], null},
//            {rGratingMovie[2][1], hh[2 * 4 + 1], null},
//            {rGratingMovie[2][2], hh[2 * 4 + 2], null},
//            {rGratingMovie[2][3], hh[2 * 4 + 3], null},

//            {m2ud, d[0][0], d[0][1]}
//            , {m4ud, d[1][1], d[1][2]}
//            , {m8ud, d[2][3], d[2][0]}
//            , {m16ud, d[3][2], d[3][0]}
//            , {fm, hFlash1, null},
//            {sinusoid, sinData, null}
        }
            , p);
        s.SECOND_FILETRING = true;

        new Thread(new Runnable() {
            public void run() {
                double t1 = System.currentTimeMillis();
//                s.simulate();
                s.DEBUG = true;
                for (int i = 0; i < 1; i++) {
//                    System.out.println("L = " + s.getValue());
                }
                double t2 = System.currentTimeMillis();
                System.out.println("Took: " + (t2 - t1) / 1000 + " seconds.");
            }

            /*
                         double[] pp = new double[p.length];
                         double dx = 1e-6;
                        public double f(int dim1, int dim2, int dx1, int dx2) {
                            for (int i = 0; i < p.length; i++) {
                                pp[i] = p[i];
                            }
                            pp[dim1] += dx1 * dx;
                            pp[dim2] += dx2 * dx;
                            double value = s.getValue(pp);
                            return value;
                        }
             */
        }).start();
    }

    /*
         DenseDoubleMatrix2D inverseVariance = new DenseDoubleMatrix2D(
        p.length, p.length);
                     for (int i = 0; i < p.length; i++) {
//                    int j = i;
        for (int j = i; j < p.length; j++) {
            double derivative = (f(i, j, 1, 1) + f(i, j, -1, -1)
                                 - f(i, j, -1, 1) - f(i, j, 1, -1)) /
                                (4 * dx * dx);
//                    System.out.println(
//                        i + ": " + VisionUtilities.format(p[i], 5) + " - " +
//                        VisionUtilities.format(Math.sqrt(1 / derivative), 5) + " : " +
//                        derivative);
            inverseVariance.set(i, j, derivative);
            inverseVariance.set(j, i, derivative);
        }
                     }
                     System.out.println(inverseVariance.toString());
         DoubleMatrix2D variance = new Algebra().inverse(inverseVariance);
                     System.out.println(variance.toString());
                     if (true) {
        return;
                     }
     */
}
