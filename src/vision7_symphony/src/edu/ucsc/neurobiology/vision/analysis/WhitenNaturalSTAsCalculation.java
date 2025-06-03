package edu.ucsc.neurobiology.vision.analysis;


import java.io.*;
import java.util.*;

import java.awt.*;
import cern.colt.matrix.impl.*;
import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.calculations.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.stimulus.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 *
 * This code takes the natural power stasand calculates what the white noise
 * stas would be.
 *
 * To understand the math, see Theunissen, FE et al., 2001.
 *
 * The sta whitening functions well, but is not optimized for speed.
 *
 * The stv whitening is not mathematicaly correct, and should not be trusted.
 *
 * @author Matthew Grivich, University of California, Santa Cruz
 *
 */
public class WhitenNaturalSTAsCalculation
    extends AbstractCalculation {

    String filePath;
    String moviePath;

    //These parameters are determined at run time.
    //Good sample values are listed here.

    //First two guesses for whitening.  Used to determine third value.
    double firstPowerCutoff = .07, secondPowerCutoff = .03;

    //A cutoff of .03 works well for small, fast neurons.
    //For large, slow neurons with few spikes, a higher cutoff is required.
    //The code now searches for an appropriate value using the following parameters.
    //This values are set by the user, not here.
    //The values are recorded here in case the config.xml gets corrupted.
    double minimumPowerCutoff = .01, maximumPowerCutoff = .3;
    double desiredSignalToNoise = 12.0;
    double desiredAccuracy = 1.0;
    int maxIterations = 10;
    int pixelsToAverage = 10;

    //It is crude and slow, but it works.  On all neurons with an STA, the final
    //signal to noise of the max pixelsToAverage will be desiredSignalToNoise plus or minus
    //desiredAccuracy


    double stimulusVariance;
    boolean calculateSTVs;
    STAFile naturalSTAFile, naturalSTVFile;


    RawMovie movie;
    int width, height, staDepth, staOffset;
    int preTrimFrames = 0, postTrimFrames = 0; //These can be used to reduce
    //the size of the output STA.
    boolean DEBUG = true;
    double refreshTime, stixelWidth, stixelHeight;
    ComplexAlgebra algebra;
    DenseDoubleMatrix2D[] qX, qY, qT, qXDagger, qYDagger, qTDagger;
    Vision app;


    public void startCalculation() throws Exception {
        app = Vision.getInstance();
        Date start = new Date();
        long startTime = System.currentTimeMillis();

        algebra = new ComplexAlgebra();
        String datasetName = new File(filePath).getName();
        String fileName = filePath + File.separator + datasetName;
        System.out.println("STAs loaded from: " + fileName);

        String[] dividedPath = StringUtil.decomposeString(filePath, File.separator);
        String whitePath = File.separator;
        for (int i = 0; i < dividedPath.length - 1; i++) {
            whitePath = whitePath + dividedPath[i] + File.separator;
        }
        whitePath = whitePath + dividedPath[dividedPath.length - 1] + "W";
        System.out.println("Whitened STAs written to: " + whitePath);

        File staF = new File(fileName + ".sta");
        naturalSTAFile = new STAFile(staF.getAbsolutePath());

        if (calculateSTVs) {
            naturalSTVFile = new STAFile(new File(fileName + ".stv").getAbsolutePath());
        }
        PowerLawMovieFile movieFile = null;

        movie = new RawMovie(moviePath, fileName + ".globals");
        new File(whitePath).mkdirs();

        STA sta0 = naturalSTAFile.getSTA(naturalSTAFile.getIDList()[0]);
        width = movie.getWidth(); //sta0.getWidth();
        height = movie.getHeight(); //sta0.getHeight();
        staDepth = naturalSTAFile.getSTADepth();
        staOffset = naturalSTAFile.getSTAOffset();
        refreshTime = movie.getRefreshTime(); //sta0.getRefreshTime();
        stixelWidth = naturalSTAFile.getStixelWidth();

        STAFile whiteSTAFile = new STAFile(
            whitePath + File.separator + dividedPath[dividedPath.length - 1] + "W" +
            ".sta", 10000,
            width, height, staDepth - preTrimFrames - postTrimFrames, staOffset,
            stixelWidth, stixelHeight, refreshTime);

        STAFile whiteSTVFile = null;
        if (calculateSTVs) {
            whiteSTVFile = new STAFile(
                whitePath + File.separator + dividedPath[dividedPath.length - 1] +
                "W" +
                ".stv", 10000,
                width, height, staDepth - preTrimFrames - postTrimFrames, staOffset,
                stixelWidth, stixelHeight, refreshTime);
        }
        generateEigenvectorMatrices();

        double[][][][] power = powerLawCalculateFullMoviePower(new PowerLawMovieFile(
            moviePath));

//        whiten.calculateFullMoviePower(fileNameRoot);
// double[][][][] covariance = whiten.read4DFile(fileNameRoot + ".lam2");

//
//        whiten.displayMovieSpectrum(covariance);
//    whiten.ftSTAs(fileNameRoot, new STAFile(fileNameRoot + ".sta001"));
// whiten.writeEigenFile(covariance, fileNameRoot + ".ft");
//        double[][][][] covariance = whiten.read4DFile(fileNameRoot + ".ft");

        //Find maximum power
        double max = 0.0;
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                for (int t = 0; t < 2 * staDepth - 2; t++) {
                    for (int p = 0; p < 2; p++) {
                        if (Math.abs(power[x][y][t][p]) > max) {
                            max = Math.abs(power[x][y][t][p]);
                        }
                    }
                }
            }
        }

        app.startProgressBar();
        minimumPowerCutoff *= max;
        maximumPowerCutoff *= max;
        firstPowerCutoff *= max;
        secondPowerCutoff *= max;
        whitenSTAs(whiteSTAFile, whiteSTVFile, power, firstPowerCutoff,
                   secondPowerCutoff, minimumPowerCutoff, maximumPowerCutoff,
                   desiredSignalToNoise, desiredAccuracy, pixelsToAverage, maxIterations,
                   stimulusVariance);
        app.endProgressBar();

        app.sendMessage(
            "Done in: " + (System.currentTimeMillis() - startTime) / 1000. + " s.");

        Vision.getInstance().getCalculationManager().calculationDone();
    }


// These matrixes do not rescale: q (qDagger V) = V
    public void generateEigenvectorMatrices() {
        qT = algebra.generateFTMatrix(2 * staDepth - 2, 1, Math.sqrt(2 * staDepth - 2));
        qTDagger = algebra.generateFTMatrix(2 * staDepth - 2, 1,
                                            Math.sqrt(2 * staDepth - 2));

        qX = algebra.generateFTMatrix(width, 1, Math.sqrt(width));
        qXDagger = algebra.generateFTMatrix(width, -1, Math.sqrt(width));

        qY = algebra.generateFTMatrix(height, 1, Math.sqrt(height));
        qYDagger = algebra.generateFTMatrix(height, -1, Math.sqrt(height));

    }


    double[] calculateTemporalMovieCovariance(
        int dataPoints, int localDepth, int x0, int y0, int x1, int y1) throws
        IOException {

        double[] covariance = new double[localDepth];
        double[] sumA = new double[localDepth];
        double[] sumB = new double[localDepth];

        ArrayList<ImageFrame> frames = new ArrayList<ImageFrame>();

        for (int i = 0; i < localDepth; i++) {
            frames.add( (ImageFrame) movie.getFrame(i).clone());

        }

        for (int i = localDepth; i < dataPoints + localDepth; i++) {
            for (int j = 0; j < localDepth; j++) {
                sumA[j] += ( (ImageFrame) frames.get(localDepth - 1)).getPixel(x0, y0, 0);
                sumB[j] += ( (ImageFrame) frames.get(localDepth - 1)).getPixel(x1, y1, 0);
                covariance[j] +=
                    ( (ImageFrame) frames.get(localDepth - 1)).getPixel(x0, y0, 0) *
                    ( (ImageFrame) frames.get(localDepth - 1 - j)).getPixel(x1, y1, 0);
            }

            frames.remove(0);
            frames.add((ImageFrame) movie.getFrame(i).clone());
        }
        for (int i = 0; i < localDepth; i++) {

            covariance[i] = ( (covariance[i] -
                               sumA[i] * sumB[i] / dataPoints) / (dataPoints - 1));
        }

        /*      double[] expandedCovariance = new double[staDepth * 2 - 2];
            for (int i = 0; i < staDepth * 2 - 2; i++) {
                if (i < staDepth) {
                    expandedCovariance[i] = covariance[i];
                } else {
                  expandedCovariance[i] = covariance[staDepth * 2 - 2 - i];
                }
            }*/
        return covariance;
    }


    public void write4DFile(double[][][][] matrix, String fileName) throws IOException {
        int nA = matrix.length;
        int nB = matrix[0].length;
        int nC = matrix[0][0].length;
        int nD = matrix[0][0][0].length;
        File file = new File(fileName);
        DataOutputStream dos = new DataOutputStream(
            new BufferedOutputStream(new FileOutputStream(file), 1024 * 1024));
        dos.writeInt(nA);
        dos.writeInt(nB);
        dos.writeInt(nC);
        dos.writeInt(nD);

        for (int a = 0; a < nA; a++) {
            for (int b = 0; b < nB; b++) {
                for (int c = 0; c < nC; c++) {
                    for (int d = 0; d < nD; d++) {
                        dos.writeDouble(matrix[a][b][c][d]);
                    }
                }

            }
        }

        dos.flush();
        dos.close();
    }


    public double[][][][] read4DFile(String fileName) throws IOException {
        System.out.println("Reading " + fileName + ".");
        File file = new File(fileName);
        double[][][][] matrix = null;
        DataInputStream dis = new DataInputStream(
            new BufferedInputStream(new FileInputStream(file), 1024 * 1024));
        int nA = dis.readInt();
        int nB = dis.readInt();
        int nC = dis.readInt();
        int nD = dis.readInt();

        matrix = new double[nA][nB][nC][nD];

        for (int a = 0; a < nA; a++) {
            for (int b = 0; b < nB; b++) {
                for (int c = 0; c < nC; c++) {
                    for (int d = 0; d < nD; d++) {
                        matrix[a][b][c][d] = dis.readDouble();
                    }
                }

            }
        }

        dis.close();

        return matrix;
    }


//Calculate the power spectrum for a power law movie.
    public double[][][][] powerLawCalculateFullMoviePower(PowerLawMovieFile movieFile) throws
        IOException {

        int spatialPeriodCutoff = movieFile.spatialPeriodCutoff;
        int temporalPeriodCutoff = movieFile.temporalPeriodCutoff;
        double spatialPower = movieFile.spatialPower;
        double temporalPower = movieFile.temporalPower;
        double standardDeviationSq = movieFile.sigma *
                                     movieFile.sigma;
        int omegaSize = 262144;
        int kXSize = width;
        int kYSize = height;

        double aspectRatio = (double) kXSize / kYSize;
        double maxSpatialNormalization = Math.pow(2.0 * kXSize / spatialPeriodCutoff *
                                                  kYSize * aspectRatio /
                                                  spatialPeriodCutoff,
                                                  spatialPower);
        double maxTemporalNormalization = Math.pow( (double) omegaSize /
            temporalPeriodCutoff,
            2 * temporalPower);
        double[][][][] powerSpect = new double[width][height][staDepth][2];
        double[] temporalPowerSpectRe = new double[omegaSize];
        double[] temporalPowerSpectIm = new double[omegaSize];
        //Calculate the power in the temporal direction
        for (int t = 0; t < omegaSize; t++) {
            if (t > 0 /* omegaSize / 2 - 1*/) {
                double temporalNormalization = Math.pow(t <
                    omegaSize / 2 + 1 ? t :
                    (omegaSize - t),
                    2 * temporalPower);
                temporalPowerSpectRe[t] = Math.min(temporalNormalization,
                    maxTemporalNormalization);
                temporalPowerSpectIm[t] = 0.0;
            } else {
                temporalPowerSpectRe[t] = 0.0;
                temporalPowerSpectIm[t] = 0.0;
            }

        }

        ScatterPlot scatter = new ScatterPlot("");
        ScatterPlot ref = new ScatterPlot("");
        for (int i = 0; i < temporalPowerSpectRe.length; i++) {
            //If frequency = nyquist or above, display it as negative
            double freq = 60 * (i <
                                omegaSize / 2 - 1 ? i :
                                - (omegaSize - i)) / (omegaSize / 2.0);

            scatter.add(i /*freq*/,
                        temporalPowerSpectRe[i] / maxTemporalNormalization);
        }

        /*
                PlotPanel panel = new PlotPanel();
                panel.addData(scatter, new ScatterPlotStyle(
                    SymbolType.NONE, 0, Color.black, true, Color.black, 1));
                panel.setLabels("frequency(Hz)", "power/maxPower");
                panel.autoscale();
         panel.setRange(panel.getRange()[0], panel.getRange()[1], 0, panel.getRange()[3]);
                PlotUtilities.showData("Temporal Power Spectrum", panel);
         */
        //Fourier transform the temporal power spectrum, which gives the temporal
        //covariance
        FFT.fft(temporalPowerSpectRe, temporalPowerSpectIm, 1);

        //caculate the spatial power spectrum, and mulitply it by the temporal
        //covariance
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {

                for (int t = 0; t < staDepth; t++) {

                    double spatialNormalization = Math.pow( (x <
                        kXSize / 2 + 1 ?
                        x * x :
                        (kXSize - x) * (kXSize - x)) +
                        (y < kYSize / 2 + 1 ? y * y * aspectRatio * aspectRatio :
                         (kYSize - y) * (kYSize - y) * aspectRatio * aspectRatio),
                        spatialPower);
                    double normalization =
                        Math.min(spatialNormalization,
                                 maxSpatialNormalization) * temporalPowerSpectRe[t];
                    powerSpect[x][y][t][0] = normalization;
                    powerSpect[x][y][t][1] = 0.0;
                }

            }
        }

        scatter = new ScatterPlot("");

        for (int i = 0; i < powerSpect.length; i++) {
            //If frequency = nyquist or above, display it as negative
            double freq = .02 * (i <
                                 powerSpect.length / 2 - 1 ? i :
                                 - (powerSpect.length - i)) / (powerSpect.length / 2.0);

            scatter.add(i /*freq*/,
                        powerSpect[i][0][0][0] / powerSpect[1][0][0][0]);

        }
        /*
                panel = new PlotPanel();
                panel.addData(scatter, new ScatterPlotStyle(
                    SymbolType.NONE, 0, Color.black, true, Color.black, 1));
                panel.setLabels("frequency(micron^(-1))", "power/maxPower");
                panel.autoscale();
         panel.setRange(panel.getRange()[0], panel.getRange()[1], 0, panel.getRange()[3]);
                PlotUtilities.showData("X Power Spectrum", panel);
         */
//        double[][] XY = new double[width][height];
//
//        for (int x = 0; x < width; x++) {
//            for (int y = 0; y < height; y++) {
//                XY[x][y] = powerSpect[x][y][0][0];powerSpect[ (x + width / 2) % width
//                           ][ (y + height / 2) % height
//                           ][0][0];

//            }
//        }

// DoubleHistogram2D hist = new DoubleHistogram2D("", -.02, .02, -.02, .02, XY);
// PlotUtilities.showData("Spatial Power Spectrum",
//                        hist,
//                        new HistogramStyle());

        //Expand the matrix in the temporal direction, per instructions of
        //whitening paper.
        //[x][y][t][Re or Im]
        double[][][][] fullPower = new double[width][height][2 *
                                   staDepth - 2][2];
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                for (int t = 0; t < 2 * staDepth - 2; t++) {
                    fullPower[x][y][t][0] = powerSpect[x]
                                            [y]
                                            [t < staDepth ? t :
                                            2 * staDepth - 2 - t]
                                            [0];
                    fullPower[x][y][t][1] = 0.0;
                }
            }
        }
        //Transform temporal direction back from covariance to power.
        DenseDoubleMatrix2D[] ftT = algebra.generateFTMatrix(2 * staDepth - 2, -1,
            Math.sqrt(2 * staDepth - 2));
//                  System.out.println(totalVariance);
        //amplitude scales like 1/multiplier, exactly.
        transformT(ftT, fullPower, 1);
//        showMatrix(fullPower, "Theoretical");

        return fullPower;
//        this.write4DFile(fullCovariance, fileNameRoot + ".lam2");
    }


    public void calculateFullMoviePower(String fileNameRoot) throws IOException {
        int localWidth = width;
        int localHeight = height;
        int localDepth = staDepth;
        int statisticsLength = 60000;

        //[x][y][OmegaT][Re or Im]
        double[][][][] fullCovariance = new double[localWidth][localHeight][
                                        localDepth][2];
        for (int x = 0; x < localWidth; x++) {
            for (int y = 0; y < localHeight; y++) {
                System.out.println("X:  " + x + " Y:  " + y);
                double[][] temporalCovariance = new double[2][];
                temporalCovariance[0] = calculateTemporalMovieCovariance(statisticsLength,
                    localDepth, 0, 0, x, y);
                temporalCovariance[1] = new double[temporalCovariance[0].length];
                Arrays.fill(temporalCovariance[1], 0.0);
                for (int t = 0; t < localDepth; t++) {
                    fullCovariance[x][y][t][0] = temporalCovariance[0][t];
                    fullCovariance[x][y][t][1] = 0.0;
                }

            }
        }

        algebra.transformT(algebra.generateFTMatrix(localDepth, 1, localDepth),
                           fullCovariance, 1);
        algebra.transformX(algebra.generateFTMatrix(localWidth, 1, localWidth),
                           fullCovariance, 1);
        algebra.transformY(algebra.generateFTMatrix(localHeight, 1, localHeight),
                           fullCovariance, 1);
//        showMatrix(fullCovariance, "Measured");

        ScatterPlot scatter = new ScatterPlot("");

        for (int i = 0; i < localDepth; i++) {

            scatter.add(i,
                        fullCovariance[0][0][i][0]);

        }

        PlotPanel panel = new PlotPanel();
        panel.addData(scatter, new ScatterPlotStyle(
            SymbolType.NONE, 0, Color.black, true, Color.black, 1));

        panel.setLabels("frequency(Hz)", "power/maxPower");

        panel.autoscale();

        PlotUtil.showData("Temporal Power Spectrum", panel);

        scatter = new ScatterPlot("");

        for (int i = 0; i < localWidth; i++) {

            scatter.add(i,
                        fullCovariance[i][0][0][0]);

        }

        panel = new PlotPanel();
        panel.addData(scatter, new ScatterPlotStyle(
            SymbolType.NONE, 0, Color.black, true, Color.black, 1));

        panel.setLabels("frequency(Hz)", "power/maxPower");

        panel.autoscale();

        PlotUtil.showData("X Power Spectrum", panel);

        scatter = new ScatterPlot("");

        for (int i = 0; i < localHeight; i++) {

            scatter.add(i,
                        fullCovariance[0][i][0][0]);

        }

        panel = new PlotPanel();
        panel.addData(scatter, new ScatterPlotStyle(
            SymbolType.NONE, 0, Color.black, true, Color.black, 1));

        panel.setLabels("frequency(Hz)", "power/maxPower");

        panel.autoscale();

        PlotUtil.showData("Y Power Spectrum", panel);

        this.write4DFile(fullCovariance, fileNameRoot + ".power");
    }


    public void transformT(DenseDoubleMatrix2D[] matrix, double[][][][][] sta,
                           double multiplier) {
        for (int c = 0; c < 3; c++) {
            transformT(matrix, sta[c], multiplier);
        }
    }


    public void transformT(DenseDoubleMatrix2D[] matrix, double[][][][] sta,
                           double multiplier) {
        DenseDoubleMatrix1D[] staSubset = new DenseDoubleMatrix1D[2];

        staSubset[0] = new DenseDoubleMatrix1D(sta[0][0].length);
        staSubset[1] = new DenseDoubleMatrix1D(sta[0][0].length);

        for (int x = 0; x < sta.length; x++) {
            for (int y = 0; y < sta[0].length; y++) {

                for (int t = 0; t < sta[0][0].length; t++) {
                    staSubset[0].set(t, sta[x][y][t][0]);
                    staSubset[1].set(t, sta[x][y][t][1]);
                }

                staSubset = algebra.matrixMultiply(matrix, staSubset);

                for (int t = 0; t < sta[0][0].length; t++) {
                    sta[x][y][t][0] = staSubset[0].get(t) * multiplier;
                    sta[x][y][t][1] = staSubset[1].get(t) * multiplier;
                }

            }
        }

    }


    public void transformX(DenseDoubleMatrix2D[] matrix, double[][][][][] sta,
                           double multiplier) {
        for (int c = 0; c < 3; c++) {
            transformX(matrix, sta[c], multiplier);
        }
    }


    public void transformX(DenseDoubleMatrix2D[] matrix, double[][][][] sta,
                           double multiplier) {
        DenseDoubleMatrix1D[] staSubset = new DenseDoubleMatrix1D[2];

        staSubset[0] = new DenseDoubleMatrix1D(sta.length);
        staSubset[1] = new DenseDoubleMatrix1D(sta.length);

        for (int t = 0; t < sta[0][0].length; t++) {
            for (int y = 0; y < sta[0].length; y++) {

                for (int x = 0; x < sta.length; x++) {
                    staSubset[0].set(x, sta[x][y][t][0]);
                    staSubset[1].set(x, sta[x][y][t][1]);
                }

                staSubset = algebra.matrixMultiply(matrix, staSubset);

                for (int x = 0; x < sta.length; x++) {
                    sta[x][y][t][0] = staSubset[0].get(x) * multiplier;
                    sta[x][y][t][1] = staSubset[1].get(x) * multiplier;
                }

            }
        }

    }


    public void transformY(DenseDoubleMatrix2D[] matrix, double[][][][][] sta,
                           double multiplier) {
        for (int c = 0; c < 3; c++) {
            transformY(matrix, sta[c], multiplier);
        }
    }


    public void transformY(DenseDoubleMatrix2D[] matrix, double[][][][] sta,
                           double multiplier) {
        DenseDoubleMatrix1D[] staSubset = new DenseDoubleMatrix1D[2];

        staSubset[0] = new DenseDoubleMatrix1D(sta[0].length);
        staSubset[1] = new DenseDoubleMatrix1D(sta[0].length);

        for (int t = 0; t < sta[0][0].length; t++) {
            for (int x = 0; x < sta.length; x++) {

                for (int y = 0; y < sta[0].length; y++) {
                    staSubset[0].set(y, sta[x][y][t][0]);
                    staSubset[1].set(y, sta[x][y][t][1]);
                }

                staSubset = algebra.matrixMultiply(matrix, staSubset);

                for (int y = 0; y < sta[0].length; y++) {
                    sta[x][y][t][0] = staSubset[0].get(y) * multiplier;
                    sta[x][y][t][1] = staSubset[1].get(y) * multiplier;
                }

            }
        }

    }


//[x][y][t][Re or Im]
    public void whitenSTAs(STAFile staFile, STAFile stvFile,
                           double[][][][] covarianceFT, double initialCutoff,
                           double secondCutoff, double minCutoff, double maxCutoff,
                           double desiredSignalToNoise, double desiredAccuracy,
                           int pixelsToAverage,
                           int maxIterations, double stimulusVariance) throws
        IOException {

        int preTrimFrames = 0; //40;
        int postTrimFrames = 0; //120;

        int[] idList = naturalSTAFile.getIDList();
        for (int i = 0; i < idList.length; i++) {

            app.sendMessage("Whitening:  " + idList[i] + " (" + (i + 1) + ")");
            STA sta = naturalSTAFile.getSTA(idList[i]);
            double[][][][][] staArray = new double[3][width][height][2 * staDepth -
                                        2][2];

            double[][][][][] tempSTAArray = new double[3][width][height][2 * staDepth -
                                            2][2];

            for (int t = 0; t < staDepth; t++) {
                STAFrame frame = (STAFrame) sta.getFrame(t);
                float[] colors = frame.getBuffer();
                for (int x = 0; x < width; x++) {
                    for (int y = 0; y < height; y++) {
                        for (int c = 0; c < 3; c++) {
                            staArray[c][x][y][t][0] = colors[3 * (y * width + x) +
                                c];

                        }
                    }
                }
            }
            double amplitudeSquared = 0.0;
            for (int x = 0; x < width; x++) {
                for (int y = 0; y < height; y++) {
                    for (int t = 0; t < staDepth * 2 - 2; t++) {
                        for (int c = 0; c < 3; c++) {
                            staArray[c][x][y][t][0] = staArray[c][x]
                                [y]
                                [t < staDepth ? t : staDepth * 2 - 2 - t]
                                [0];
                            amplitudeSquared += staArray[c][x][y][t][0] *
                                staArray[c][x][y][t][0];
                        }

                    }
                }
            }

            transformT(qTDagger, staArray, 1.0);
            transformX(qXDagger, staArray, 1.0);
            transformY(qYDagger, staArray, 1.0);

            STA newSTA = null;
            int j = 0;
            double currentSignalToNoise = -1.0;
            double previousSignalToNoise = -1.0;
            double currentCutoff = initialCutoff;
            double previousCutoff = -1;
            double rms = -1;

            while (Math.abs(currentSignalToNoise - desiredSignalToNoise) >
                   desiredAccuracy && j < maxIterations) {
                if (DEBUG) {
                    System.out.println("j: " + j);
                    System.out.println("currentCutoff: " + currentCutoff);
                }

                //make copy of array so that we don't keep transforming the same matrix
                for (int x = 0; x < width; x++) {
                    for (int y = 0; y < height; y++) {
                        for (int t = 0; t < staDepth * 2 - 2; t++) {
                            for (int c = 0; c < 3; c++) {
                                tempSTAArray[c][x][y][t][0] = staArray[c][x][y][t][0];
                                tempSTAArray[c][x][y][t][1] = staArray[c][x][y][t][1];
                            }
                        }
                    }
                }

                for (int c = 0; c < tempSTAArray.length; c++) {
                    for (int x = 0; x < tempSTAArray[0].length; x++) {
                        for (int y = 0; y < tempSTAArray[0][0].length; y++) {
                            for (int t = 0; t < tempSTAArray[0][0][0].length; t++) {
                                if (covarianceFT[x][y][t][0] > currentCutoff ||
                                    covarianceFT[x][y][t][0] < -currentCutoff) {
                                    tempSTAArray[c][x][y][t][0] /= covarianceFT[x][y][t][
                                        0];
                                    //imaginary part of covariance FT is zero.
                                    tempSTAArray[c][x][y][t][1] /= covarianceFT[x][y][t][
                                        0];
                                } else {
                                    tempSTAArray[c][x][y][t][0] = 0.0;
                                    tempSTAArray[c][x][y][t][1] = 0.0;

                                }

                            }

                        }
                    }
                }

                transformT(qT, tempSTAArray, 1.0);
                transformX(qX, tempSTAArray, 1.0);
                transformY(qY, tempSTAArray, 1.0);

                newSTA = new STA(staDepth - preTrimFrames - postTrimFrames,
                                 sta.getWidth(), sta.getHeight(), sta.getRefreshTime(),
                                 sta.getStixelWidth(), sta.getStixelHeight());

                double amplitudeSquaredAfter = 0.0;
                for (int c = 0; c < tempSTAArray.length; c++) {
                    for (int x = 0; x < tempSTAArray[0].length; x++) {
                        for (int y = 0; y < tempSTAArray[0][0].length; y++) {
                            for (int t = 0; t < tempSTAArray[0][0][0].length; t++) {
                                amplitudeSquaredAfter += tempSTAArray[c][x][y][t][0] *
                                    tempSTAArray[c][x][y][t][0];
                            }
                        }
                    }
                }
                double scaler = Math.sqrt(amplitudeSquared / amplitudeSquaredAfter);
                MeanVarianceCalculator mvc = new MeanVarianceCalculator();
                for (int t = preTrimFrames; t < staDepth - postTrimFrames; t++) {
                    STAFrame newFrame = (STAFrame) newSTA.getFrame(t - preTrimFrames);
                    float[] colors = newFrame.getBuffer();

                    for (int x = 0; x < width; x++) {
                        for (int y = 0; y < height; y++) {
                            for (int c = 0; c < 3; c++) {
                                colors[3 * (y * width + x) +
                                    c] = (float) tempSTAArray[c][x][y][t][0] *
                                         (float) scaler;
                                if (t == preTrimFrames && j == 0) {
                                    mvc.add( (float) tempSTAArray[c][x][y][t][0] *
                                            (float) scaler);
                                }

                            }
                        }
                    }

                    newFrame.setBuffer(colors);
                    newFrame.setErrorBuffer( ( (STAFrame) sta.getFrame(t)).getErrorBuffer());
                    newSTA.setFrame(t - preTrimFrames, newFrame);
                }

                int[] params = newSTA.getMainFrameParams();
//                this.showSTA(idList[i] + " " + currentCutoff, newSTA, params[0], false);

                if (j == 0) {
                    rms = mvc.getStandardDeviation();
                }

                previousSignalToNoise = currentSignalToNoise;
                currentSignalToNoise = Math.sqrt(newSTA.getFrame(params[0]).
                                                 getPeakMagnitudeSquared(pixelsToAverage)) /
                                       rms;

                if (j == 0) {
                    previousCutoff = currentCutoff;
                    currentCutoff = secondCutoff;

                } else {

                    double stepSize = (currentCutoff - previousCutoff) /
                                      (currentSignalToNoise - previousSignalToNoise);

                    previousCutoff = currentCutoff;

                    //This removes bad behaviour when the signal to noise function
                    //is not monotonic and cutoff is too low.
                    //If cutoff is too high, it will not work and the neuron
                    //will be very blurry in the final result, because the cutoff
                    //will hit the maximum value.
                    if (stepSize <= 0) {
                        currentCutoff = currentCutoff +
                                        ( (desiredSignalToNoise - currentSignalToNoise >
                                           0)
                                         ? currentCutoff * .2 : -currentCutoff * .2);
                        currentCutoff = Math.min(Math.max(currentCutoff, minCutoff),
                                                 maxCutoff);

                    } else {
                        currentCutoff = Math.min(Math.max( (desiredSignalToNoise -
                            currentSignalToNoise) *
                            stepSize + currentCutoff, minCutoff), maxCutoff);
                    }

                }

                if (DEBUG) {
                    System.out.println("currentSignalToNoise: " + currentSignalToNoise);
                    System.out.println();
                }
                // This will only occur if the min or max cutoff is used twice in a row.
                // That is, the search procedure failed.  The min or max cutoff sta
                // will be saved.
                if (currentCutoff == previousCutoff) {
                    break;
                }

                j++;
            }
            staFile.addSTA(idList[i], newSTA);
            if (calculateSTVs) {
                stvFile.addSTA(idList[i], whitenSTA(naturalSTVFile.getSTA(idList[i])
                    , covarianceFT, currentCutoff, stimulusVariance));
            }
            app.setProgress( (100 * i) / idList.length);
        }
        staFile.close();
    }


    //[x][y][t][Re or Im]
    public STA whitenSTA(STA sta, double[][][][] covarianceFT, double cutoff,
                         double offset) throws
        IOException {

        int preTrimFrames = 0; //40;
        int postTrimFrames = 0; //120;

        double[][][][][] staArray = new double[3][width][height][2 * staDepth -
                                    2][2];

        for (int t = 0; t < staDepth; t++) {
            STAFrame frame = (STAFrame) sta.getFrame(t);
            float[] colors = frame.getBuffer();
            for (int x = 0; x < width; x++) {
                for (int y = 0; y < height; y++) {
                    for (int c = 0; c < 3; c++) {
                        staArray[c][x][y][t][0] = colors[3 * (y * width + x) +
                                                  c] - offset;

                    }
                }
            }
        }
        double amplitudeSquared = 0.0;
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                for (int t = 0; t < staDepth * 2 - 2; t++) {
                    for (int c = 0; c < 3; c++) {
                        staArray[c][x][y][t][0] = staArray[c][x]
                                                  [y]
                                                  [t < staDepth ? t :
                                                  staDepth * 2 - 2 - t]
                                                  [0];
                        amplitudeSquared += staArray[c][x][y][t][0] *
                            staArray[c][x][y][t][0];
                    }

                }
            }
        }

        transformT(qTDagger, staArray, 1.0);
        transformX(qXDagger, staArray, 1.0);
        transformY(qYDagger, staArray, 1.0);

        STA newSTA = null;

        for (int c = 0; c < staArray.length; c++) {
            for (int x = 0; x < staArray[0].length; x++) {
                for (int y = 0; y < staArray[0][0].length; y++) {
                    for (int t = 0; t < staArray[0][0][0].length; t++) {
                        if (covarianceFT[x][y][t][0] > cutoff ||
                            covarianceFT[x][y][t][0] < -cutoff) {
                            staArray[c][x][y][t][0] /= covarianceFT[x][y][t][
                                0];
                            //imaginary part of covariance FT is zero.
                            staArray[c][x][y][t][1] /= covarianceFT[x][y][t][
                                0];
                        } else {
                            staArray[c][x][y][t][0] = 0.0;
                            staArray[c][x][y][t][1] = 0.0;

                        }

                    }

                }
            }
        }

        transformT(qT, staArray, 1.0);
        transformX(qX, staArray, 1.0);
        transformY(qY, staArray, 1.0);

        newSTA = new STA(staDepth - preTrimFrames - postTrimFrames,
                         sta.getWidth(), sta.getHeight(), sta.getRefreshTime(),
                         sta.getStixelWidth(), sta.getStixelHeight());

        double amplitudeSquaredAfter = 0.0;
        for (int c = 0; c < staArray.length; c++) {
            for (int x = 0; x < staArray[0].length; x++) {
                for (int y = 0; y < staArray[0][0].length; y++) {
                    for (int t = 0; t < staArray[0][0][0].length; t++) {
                        amplitudeSquaredAfter += staArray[c][x][y][t][0] *
                            staArray[c][x][y][t][0];
                    }
                }
            }
        }
        double scaler = Math.sqrt(amplitudeSquared / amplitudeSquaredAfter);

        for (int t = preTrimFrames; t < staDepth - postTrimFrames; t++) {
            STAFrame newFrame = (STAFrame) newSTA.getFrame(t - preTrimFrames);
            float[] colors = newFrame.getBuffer();

            for (int x = 0; x < width; x++) {
                for (int y = 0; y < height; y++) {
                    for (int c = 0; c < 3; c++) {
                        colors[3 * (y * width + x) +
                            c] = (float) (staArray[c][x][y][t][0] *
                                          scaler + offset);

                    }
                }
            }

            newFrame.setBuffer(colors);
            newFrame.setErrorBuffer( ( (STAFrame) sta.getFrame(t)).getErrorBuffer());
            newSTA.setFrame(t - preTrimFrames, newFrame);
        }

//                this.showSTA("STV", newSTA, sta.getMainFrame(), false);

        return newSTA;

    }


    public void displayMovieSpectrum(double[][][][] powerSpectrum) {

        double[][] XYre = new double[2 * width - 2][2 * width - 2];
        double[][] XYim = new double[2 * width - 2][2 * width - 2];
        for (int x = 0; x < 2 * width - 2; x++) {
            for (int y = 0; y < 2 * height - 2; y++) {
                XYre[x][y] = powerSpectrum[x][y][2 * staDepth - 2 - 1][0];
                XYim[x][y] = powerSpectrum[x][y][2 * staDepth - 2 - 1][1];
            }
        }

        PlotUtil.showData("Spatial Power Spectrum",
                          new DoubleHistogram2D("", -1, 1, -1, 1, XYre),
                          new HistogramStyle());
        PlotUtil.showData("Spatial Power Spectrum",
                          new DoubleHistogram2D("", -1, 1, -1, 1, XYim),
                          new HistogramStyle());

        double average = 0.0;
        double[] xFT = new double[2 * width - 2];
        for (int x = 0; x < 2 * width - 2; x++) {
            xFT[x] = Math.sqrt(powerSpectrum[x][0][0][0] * powerSpectrum[x][0][0][0] +
                               powerSpectrum[x][0][0][1] * powerSpectrum[x][0][0][1]); // imaginary part is zero
            average += powerSpectrum[x][0][0][0];
        }
        average /= 2 * width - 2;
        for (int x = 0; x < 2 * width - 2; x++) {
            xFT[x] /= average;
        }
        ScatterPlot xScatter = new ScatterPlot("");
        ScatterPlot ref = new ScatterPlot("");
        for (int i = 1; i < xFT.length; i++) {
            double omega = Math.PI * 2 * (double) i / (double) powerSpectrum.length;
            xScatter.add(Math.log(omega) / Math.log(10.0),
                         Math.log(xFT[i]) / Math.log(10.0));
            ref.add(Math.log(omega) / Math.log(10.0),
                    Math.log(1.0 / (double) omega) / Math.log(10.0));
        }

        PlotPanel panel = new PlotPanel();
        panel.addData(xScatter, new ScatterPlotStyle(
            SymbolType.NONE, 0, Color.black, true, Color.black, 1));

        panel.addData(ref, new ScatterPlotStyle(
            SymbolType.NONE, 0, Color.red, true, Color.red, 1));
        panel.setLabels("log(frequency(rad))", "log(power/averagePower)");

        panel.autoscale();
//        panel.setRange(0, 1.5, panel.getRange()[2], panel.getRange()[3]);
        PlotUtil.showData("X Power Spectrum", panel);

        average = 0.0;
        double[] tFT = new double[2 * staDepth - 2];
        for (int t = 0; t < 2 * staDepth - 2; t++) {
            tFT[t] = Math.sqrt(powerSpectrum[0][0][t][0] * powerSpectrum[0][0][t][0] +
                               powerSpectrum[0][0][t][1] * powerSpectrum[0][0][t][1]); // imaginary part is zero
            average += powerSpectrum[0][0][t][0];
        }
        average /= 2 * staDepth - 2;
        for (int t = 0; t < 2 * staDepth - 2; t++) {
            tFT[t] /= average;
        }
        ScatterPlot tScatter = new ScatterPlot("");
        ref = new ScatterPlot("");
        for (int i = 1; i < tFT.length; i++) {
            double omega = Math.PI * 2 * (double) i / (double) powerSpectrum.length;
            tScatter.add(Math.log(omega) / Math.log(10.0),
                         Math.log(tFT[i]) / Math.log(10.0));
            ref.add(Math.log(omega) / Math.log(10.0),
                    Math.log(1.0 / (double) omega) / Math.log(10.0));
        }

        panel = new PlotPanel();
        panel.addData(tScatter, new ScatterPlotStyle(
            SymbolType.NONE, 0, Color.black, true, Color.black, 1));

        panel.addData(ref, new ScatterPlotStyle(
            SymbolType.NONE, 0, Color.red, true, Color.red, 1));
        panel.setLabels("log(frequency(rad))", "log(power/averagePower)");

        panel.autoscale();
//        panel.setRange(0, 1.5, panel.getRange()[2], panel.getRange()[3]);
        PlotUtil.showData("Temporal Power Spectrum", panel);

        double[][] ftXY = new double[2 * height - 2][2 * height - 2];
        for (int x = 0; x < 2 * height - 2; x++) {
            for (int y = 0; y < 2 * height - 2; y++) {
                ftXY[x][y] = Math.log(powerSpectrum[2 * ( (x + 32)) %
                                      (2 * width - 2)][ (y + 32) % (2 * height - 2)][0][0]) /
                             Math.log(10.0);
//                ftXY[x][y] = Math.log(covariance[x][y][0][0])/Math.log(10.0);
//                 ftXY[x][y] = covariance[x][y][0][0];

            }
        }

        PlotUtil.showData("Spatial Power Spectrum",
                          new DoubleHistogram2D("", -1, 1, -1, 1, ftXY),
                          new HistogramStyle());
        xScatter = new ScatterPlot("");
        ref = new ScatterPlot("");
        for (int i = 1; i < ftXY.length; i++) {
            double omega = Math.PI * 2 * (double) i / (double) powerSpectrum.length;
            xScatter.add(Math.log(omega) / Math.log(10.0), xFT[i]);
            ref.add(Math.log(omega) / Math.log(10.0),
                    Math.log(1.0 / (double) omega) / Math.log(10.0));

        }

    }


//    public void showMatrix(double[][][][] matrix, String label) throws IOException {
//        STA sta = new STA(matrix[0][0].length, matrix.length, matrix[0].length, 1, 1, 50);
//        for (int t = 0; t < matrix[0][0].length; t++) {
//            STAFrame frame = (STAFrame) sta.getFrame(t);
//            float[] colors = frame.getBuffer();
//            for (int x = 0; x < matrix.length; x++) {
//                for (int y = 0; y < matrix[0].length; y++) {
//                    for (int c = 0; c < 3; c++) {
//
//                        colors[3 * (y * matrix.length + x) +
//                            c] = (float) matrix[x][y][t][0];
//
//                    }
//                }
//            }
//
//            frame.setBuffer(colors);
//            sta.setFrame(t, frame);
//        }
////                 No longer works: showSTA functianlity changed
//        showSTA(label + ": Real", sta, 0, false);
//
//        sta = new STA(matrix[0][0].length, matrix.length, matrix[0].length, 1, 1, 50);
//        for (int t = 0; t < matrix[0][0].length; t++) {
//            STAFrame frame = (STAFrame) sta.getFrame(t);
//            float[] colors = frame.getBuffer();
//            for (int x = 0; x < matrix.length; x++) {
//                for (int y = 0; y < matrix[0].length; y++) {
//                    for (int c = 0; c < 3; c++) {
//                        colors[3 * (y * matrix.length + x) +
//                            c] = (float) matrix[x][y][t][1];
//
//                    }
//                }
//            }
//
//            frame.setBuffer(colors);
//            sta.setFrame(t, frame);
//        }
//
//        showSTA(label + ": Imaginary", sta, 0, false);
//
//    }


    //pasted the showSTA code here, because it got deleted from VisionUtilities.
    // Not called from anywhere
//    public static PlotPanel showSTA(
//        String name, final STA sta, int frameIndex, boolean includeContour) throws
//        IOException {
//
//        ArrayList<JFrame> frameList = new ArrayList<JFrame>();
//        final double significance = 5;
//
//        JFrame f = new JFrame(name);
//        frameList.add(f);
//
//
//        final PlotPanel p = new PlotPanel();
//        p.setXAxisVisible(false);
//        p.setYAxisVisible(false);
//        p.setRange(0, sta.getStixelWidth() * sta.getWidth(), 0, sta.getStixelHeight() * sta.getHeight());
//
//        // add the STA
//        p.addDataPlotter(new ColorPlotPlotter());
//        final SimpleColorPlot scp = new SimpleColorPlot(SimpleColorPlot.NORMALIZE);
//        scp.setFrame(sta, frameIndex);
//        final ColorPlotStyle style = new ColorPlotStyle("STA");
//        p.addData(scp, style);
//
//        // add the control panel
//        JPanel control = new JPanel(new GridLayout(1, 0));
//
//        final JSpinner frameBox = new JSpinner(new SpinnerNumberModel(frameIndex, 0,
//            sta.size() - 1, 1));
//        frameBox.addChangeListener(new ChangeListener() {
//            public void stateChanged(ChangeEvent event) {
//                try {
//                    scp.setFrame(sta, ( (Integer) frameBox.getValue()).intValue());
//                } catch (IOException e) {
//                    Vision.reportException(e);
//                }
//                p.replotAllData();
//            }
//        });
//        control.add(frameBox);
//
//        ParametersTable table = new ParametersTable();
//        control.add(new JScrollPane(table));
//        f.add(control, BorderLayout.NORTH);
//
//        // add the countour
//        if (includeContour) {
//            p.addDataPlotter(new PolygonPlotter());
//            Polygon2DAdapter p2d = new Polygon2DAdapter();
//            p2d.setPolygon(sta.getContour(5.0, 1, true));
//            p.addData(p2d, new FunctionStyle("Contour", Color.red, 2));
//
//            // add the significant pixels markers
//            ScatterPlot sp = new ScatterPlot();
//            ImageFrame ff = sta.getMainFrame();
//            float[] intensity = new float[3];
//            float[] error = new float[3];
//            for (int y = 0; y < ff.getHeight(); y++) {
//                for (int x = 0; x < ff.getWidth(); x++) {
//                    ff.getPixel(x, y, intensity);
//                    ff.getPixelError(x, y, error);
//                    if (Math.abs(intensity[1] / error[1]) > significance) {
//                        sp.add(
//                            (x + 0.5) * ff.getStixelWidth(),
//                            (ff.getHeight() - y - 0.5) * ff.getStixelHeight());
//                    }
//                }
//            }
//
//            p.addData(sp,
//                      new ScatterPlotStyle(SymbolType.SQUARE, 3,
//                                           Color.black, false, Color.black, 1));
//
//            PlotPanel p1 = new PlotPanel();
//            double[][] timeFilter = sta.getTimeFilters(significance);
//            for (int cIndex = 0; cIndex < 3; cIndex++) {
//                ScatterPlotStyle style1 = new ScatterPlotStyle(
//                    SymbolType.SQUARE, 3, PlotUtil.rgb[cIndex], true,
//                    PlotUtil.rgb[cIndex], 1);
//                p1.addData(new ScatterPlot(timeFilter[cIndex]), style1);
//            }
//            p1.autoscale();
//
//            JSplitPane splitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, p, p1);
//            splitPane.setDividerLocation(350);
//            f.add(splitPane, BorderLayout.CENTER);
//            f.setBounds(100, 100, 400, 400);
//        } else {
//            f.add(p, BorderLayout.CENTER);
//            f.setBounds(100, 100, 400, 400);
//        }
//
//        // show the frame
//        f.setVisible(true);
//
//        return p;
//    }


    public void setParameters(HashMap parameters) {
        filePath = new File( (String) parameters.get("File Path")).getAbsolutePath();
        moviePath = new File( (String) parameters.get("Movie Path")).getAbsolutePath();
        firstPowerCutoff = Double.parseDouble( (String) parameters.get(
            "First Power Cutoff"));
        secondPowerCutoff = Double.parseDouble( (String) parameters.get(
            "Second Power Cutoff"));
        minimumPowerCutoff = Double.parseDouble( (String) parameters.get(
            "Minimum Power Cutoff"));
        maximumPowerCutoff = Double.parseDouble( (String) parameters.get(
            "Maximum Power Cutoff"));
        desiredSignalToNoise = Double.parseDouble( (String) parameters.get(
            "Desired Signal to Noise"));
        desiredAccuracy = Double.parseDouble( (String) parameters.get("Desired Accuray"));
        pixelsToAverage = Integer.parseInt( (String) parameters.get("Pixels to Average"));
        maxIterations = Integer.parseInt( (String) parameters.get("Maximum Iterations"));

        calculateSTVs = Boolean.valueOf( (String) parameters.get("Calculate STVs")).
                        booleanValue();
        stimulusVariance = Double.parseDouble( (String) parameters.get(
            "Stimulus Variance"));

    }

}
