package edu.ucsc.neurobiology.vision.analysis;

import java.io.*;
import java.util.*;

import java.awt.*;

import cern.colt.matrix.impl.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.stimulus.*;


/**
 * This code has not been implemented into the gui.  Must currently be run
 * from JBuilder.
 * Used only for the Voss-McCartney algorithm, which is not currently used in
 * experiments.
 * 
 *
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class WhitenPinkSTAsCalculation {

    public RawMovie movie;
    public int width, height, staDepth, staOffset;
    public double refreshTime, stixelWidth, stixelHeight;
    public ComplexAlgebra algebra;
    public DenseDoubleMatrix2D[] qX, qY, qT, qXDagger, qYDagger, qTDagger;


    public WhitenPinkSTAsCalculation() {
    }


    double[] calculateTemporalMovieCovariance(int dataPoints, int staDepth, int x0,
                                              int y0,
                                              int x1, int y1) throws IOException {

        double[] covariance = new double[staDepth];
        double[] sumA = new double[staDepth];
        double[] sumB = new double[staDepth];

        ArrayList<ImageFrame> frames = new ArrayList<ImageFrame>();

        for (int i = 0; i < staDepth; i++) {
            frames.add( (ImageFrame) movie.getFrame(i).clone());

        }

        for (int i = staDepth; i < dataPoints + staDepth; i++) {
            for (int j = 0; j < staDepth; j++) {
                sumA[j] += ( (ImageFrame) frames.get(staDepth - 1)).getPixel(x0, y0, 0);
                sumB[j] += ( (ImageFrame) frames.get(staDepth - 1)).getPixel(x1, y1, 0);
                covariance[j] +=
                    ( (ImageFrame) frames.get(staDepth - 1)).getPixel(x0, y0, 0) *
                    ( (ImageFrame) frames.get(staDepth - 1 - j)).getPixel(x1, y1, 0);
            }

            frames.remove(0);
            frames.add((ImageFrame) movie.getFrame(i).clone());
        }
        for (int i = 0; i < staDepth; i++) {

            covariance[i] = ( (covariance[i] -
                               sumA[i] * sumB[i] / dataPoints) / (dataPoints - 1));
        }

        double[] expandedCovariance = new double[staDepth * 2 - 2];
        for (int i = 0; i < staDepth * 2 - 2; i++) {
            if (i < staDepth) {
                expandedCovariance[i] = covariance[i];
            } else {
                expandedCovariance[i] = staDepth * 2 - 2 - i;
            }
        }
        return expandedCovariance;
    }


    public void write4DFile(double[][][][] matrix, String fileName) {
        int nA = matrix.length;
        int nB = matrix[0].length;
        int nC = matrix[0][0].length;
        int nD = matrix[0][0][0].length;
        File file = new File(fileName);
        try {
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
        } catch (IOException e) {
            e.printStackTrace();
        }

    }


    public double[][][][] read4DFile(String fileName) {
        System.out.println("Reading " + fileName + ".");
        File file = new File(fileName);
        double[][][][] matrix = null;
        try {
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

        } catch (IOException e) {
            e.printStackTrace();
        }
        return matrix;
    }


    public double[][][][] quickCalculateFullMovieCovariance(String fileNameRoot) {
        PinkMovieFile movieFile = null;
        try {
            movieFile = new PinkMovieFile(fileNameRoot + ".rawMovie");
        } catch (IOException e) {
            e.printStackTrace();
        }
        int spatialFreqCount = movieFile.spatialFreqCount;
        int temporalFreqCount = movieFile.tempFreqCount;
        double standardDeviationSq = movieFile.standardDeviation * 2 *
                                     movieFile.standardDeviation * 2;

        //[x][y][t][Re or Im]
        double[][][][] fullCovariance = new double[2 * width - 2][2 * height - 2][2 *
                                        staDepth - 2][2];

        for (int x = 0; x < width; x++) {

            for (int y = 0; y < height; y++) {
                double sCov = 0.0;
                for (int i = 0; i < spatialFreqCount; i++) {
                    sCov += Math.max( (double) ( (1 << i) - x) / (double) (1 << i), 0.0) *
                        Math.max( (double) ( (1 << i) - y) / (double) (1 << i), 0.0);
                }
                sCov /= spatialFreqCount;

                for (int t = 0; t < staDepth; t++) {
                    double tCov = 0.0;
                    for (int i = 0; i < temporalFreqCount; i++) {
                        tCov +=
                            Math.max( (double) ( (1 << i) - t) / (double) (1 << i), 0.0);
                    }
                    tCov /= temporalFreqCount;

                    double totalCov = sCov * tCov * standardDeviationSq;
                    fullCovariance[x][y][t][0] = totalCov;
                }

            }
        }

        for (int x = 0; x < 2 * width - 2; x++) {
            for (int y = 0; y < 2 * height - 2; y++) {
                for (int t = 0; t < 2 * staDepth - 2; t++) {
                    fullCovariance[x][y][t][0] = fullCovariance[x < width ? x :
                                                 2 * width - 2 - x]
                                                 [y < height ? y : 2 * height - 2 - y]
                                                 [t < staDepth ? t :
                                                 2 * staDepth - 2 - t]
                                                 [0];
                    fullCovariance[x][y][t][1] = 0.0;
                }
            }
        }
        return fullCovariance;
//        this.write4DFile(fullCovariance, fileNameRoot + ".lam2");
    }


    public double[][][][] powerLawCalculateFullMovieCovariance(String fileNameRoot) {
        PowerLawMovieFile movieFile = null;
        try {
            movieFile = new PowerLawMovieFile(fileNameRoot + ".rawMovie");
        } catch (IOException e) {
            e.printStackTrace();
        }
        int spatialPeriodCutoff = movieFile.spatialPeriodCutoff;
        int temporalPeriodCutoff = movieFile.temporalPeriodCutoff;
        double spatialPower = movieFile.spatialPower;
        double temporalPower = movieFile.temporalPower;
        double standardDeviationSq = movieFile.sigma *
                                     movieFile.sigma;
        int omegaSize = staDepth;
        int kXSize = width;
        int kYSize = height;
        double totalVariance = 0.0;
        double aspectRatio = (double) kXSize / kYSize;
        double maxSpatialNormalization = Math.pow(2.0 * kXSize / spatialPeriodCutoff *
                                                  kYSize * aspectRatio /
                                                  spatialPeriodCutoff,
                                                  .5 * spatialPower);
        double maxTemporalNormalization = Math.pow( (double) omegaSize /
            temporalPeriodCutoff,
            temporalPower);
        double[][][][] powerSpect = new double[width][height][staDepth][2];
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                //All waves going backwards in time are zero.
                for (int t = 0; t < omegaSize; t++) {
                    if (t > omegaSize / 2 - 1) {
                        double temporalNormalization = Math.pow(t <
                            omegaSize / 2 + 1 ? t :
                            (omegaSize - t),
                            temporalPower);
                        double spatialNormalization = Math.pow( (x <
                            kXSize / 2 + 1 ?
                            x * x :
                            (kXSize - x) * (kXSize - x)) +
                            (y < kYSize / 2 + 1 ? y * y * aspectRatio * aspectRatio :
                             (kYSize - y) * (kYSize - y) * aspectRatio * aspectRatio),
                            .5 * spatialPower);
                        double normalization = Math.min(temporalNormalization,
                            maxTemporalNormalization) *
                                               Math.min(spatialNormalization,
                            maxSpatialNormalization);
                        powerSpect[x][y][t][0] = normalization;
                        powerSpect[x][y][t][1] = 0.0;
                        totalVariance += .5 * normalization *
                            normalization; // sigma from one sine wave
                    } else {
                        powerSpect[x][y][t][0] = 0.0;
                        powerSpect[x][y][t][1] = 0.0;
                    }
                }

            }
        }

        DenseDoubleMatrix2D[] ftX = algebra.generateFTMatrix(width, 1, width);
        DenseDoubleMatrix2D[] ftY = algebra.generateFTMatrix(height, 1, height);
        DenseDoubleMatrix2D[] ftT = algebra.generateFTMatrix(staDepth, 1, staDepth);
        System.out.println("Tranforming powerSpect to covariance matrix.");
        transformX(ftX, powerSpect, 1);
        transformY(ftY, powerSpect, 1);
        transformT(ftT, powerSpect, 1);

        System.out.println("Expanding Covariance Matrix.");
        //[x][y][t][Re or Im]
        double[][][][] fullCovariance = new double[2 * width - 2][2 * height - 2][2 *
                                        staDepth - 2][2];
        for (int x = 0; x < 2 * width - 2; x++) {
            for (int y = 0; y < 2 * height - 2; y++) {
                for (int t = 0; t < 2 * staDepth - 2; t++) {
                    fullCovariance[x][y][t][0] = powerSpect[x < width ? x :
                                                 2 * width - 2 - x]
                                                 [y < height ? y : 2 * height - 2 - y]
                                                 [t < staDepth ? t :
                                                 2 * staDepth - 2 - t]
                                                 [0];
                    fullCovariance[x][y][t][1] = 0.0;
                }
            }
        }
        return fullCovariance;
//        this.write4DFile(fullCovariance, fileNameRoot + ".lam2");
    }


    public void calculateFullMovieCovariance(String fileNameRoot, STAFile pinkSTAFile) throws
        IOException {

        //[x][y][OmegaT][Re or Im]
        double[][][][] fullCovariance = new double[2 * width - 2][2 * height - 2][2 *
                                        staDepth - 2][2];
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                System.out.println("X:  " + x + " Y:  " + y);
                double[][] temporalCovariance = new double[2][];
                temporalCovariance[0] = calculateTemporalMovieCovariance(20000,
                    staDepth, 0, 0, x, y);
                temporalCovariance[1] = new double[temporalCovariance[0].length];
                Arrays.fill(temporalCovariance[1], 0.0);
                for (int t = 0; t < 2 * staDepth - 2; t++) {
                    fullCovariance[x][y][t][0] =
                        fullCovariance[ (2 * width - 2 - x) % (2 * width - 2)][y][t][0] =
                            fullCovariance[x][ (2 * height - 2 - y) %
                            (2 * height - 2)][t][0] =
                                fullCovariance[ (2 * width - 2 - x) %
                                (2 * width - 2)][ (2 * height - 2 - y) %
                                (2 * height - 2)][t][0] =
                                    temporalCovariance[0][t];

                    //all imaginary covariances are actually zero - code can be cleaned
                    fullCovariance[x][y][t][1] =
                        fullCovariance[ (2 * width - 2 - x) % (2 * width - 2)][y][t][1] =
                            fullCovariance[x][ (2 * height - 2 - y) %
                            (2 * height - 2)][t][1] =
                                fullCovariance[ (2 * width - 2 - x) %
                                (2 * width - 2)][ (2 * height - 2 - y) %
                                (2 * height - 2)][t][1] =
                                    temporalCovariance[1][t];

                }

            }
        }
        this.write4DFile(fullCovariance, fileNameRoot + ".lam");
    }


    public void generateEigenvectorMatrices() {
        System.out.println("Calculating eigenvector matrices.");

        qT = algebra.generateFTMatrix(2 * staDepth - 2, 1, Math.sqrt(2 * staDepth - 2));
        qTDagger = algebra.generateFTMatrix(2 * staDepth - 2, 1,
                                            Math.sqrt(2 * staDepth - 2));

        qX = algebra.generateFTMatrix(2 * width - 2, 1, Math.sqrt(2 * width - 2));
        qXDagger = algebra.generateFTMatrix(2 * width - 2, -1, Math.sqrt(2 * width - 2));

        qY = algebra.generateFTMatrix(2 * height - 2, 1, Math.sqrt(2 * height - 2));
        qYDagger = algebra.generateFTMatrix(2 * height - 2, -1, Math.sqrt(2 * height - 2));

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

                staSubset = matrixMultiply(matrix, staSubset);

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

                staSubset = matrixMultiply(matrix, staSubset);

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

                staSubset = matrixMultiply(matrix, staSubset);

                for (int y = 0; y < sta[0].length; y++) {
                    sta[x][y][t][0] = staSubset[0].get(y) * multiplier;
                    sta[x][y][t][1] = staSubset[1].get(y) * multiplier;
                }

            }
        }

    }


    //[x][y][t][Re or Im]
    public void whitenSTAs(String fileNameRoot, STAFile pinkSTAFile,
                           double[][][][] covarianceFT, double cutoff) {
        double cutoffO = cutoff;
        try {
            STAFile staFile = new STAFile(
                fileNameRoot + ".sta", 10000,
                width, height, staDepth, staOffset,
                stixelWidth, stixelHeight, refreshTime);

            int[] idList = pinkSTAFile.getIDList();
            for (int i = 0; i < idList.length; i++) {
                cutoff = (i + 1) * cutoffO;
                System.out.println("Projecting:  " + idList[i]);
//                STA sta = pinkSTAFile.getSTA(idList[i]);
                STA sta = pinkSTAFile.getSTA(863);
                double[][][][][] staArray = new double[3][2 * width - 2][2 * height -
                                            2][staDepth * 2 - 2][2];

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

                for (int x = 0; x < 2 * width - 2; x++) {
                    for (int y = 0; y < 2 * height - 2; y++) {
                        for (int t = 0; t < staDepth * 2 - 2; t++) {
                            staArray[0][x][y][t][0] = staArray[0][x < width ? x :
                                width * 2 - 2 - x]
                                [y < height ? y : height * 2 - 2 - y]
                                [t < staDepth ? t : staDepth * 2 - 2 - t]
                                [0];
                            staArray[1][x][y][t][0] = staArray[1][x < width ? x :
                                width * 2 - 2 - x]
                                [y < height ? y : height * 2 - 2 - y]
                                [t < staDepth ? t : staDepth * 2 - 2 - t]
                                [0];
                            staArray[2][x][y][t][0] = staArray[2][x < width ? x :
                                width * 2 - 2 - x]
                                [y < height ? y : height * 2 - 2 - y]
                                [t < staDepth ? t : staDepth * 2 - 2 - t]
                                [0];

                        }
                    }
                }
                transformT(qTDagger, staArray, 1.0);
                transformX(qXDagger, staArray, 1.0);
                transformY(qYDagger, staArray, 1.0);

                for (int c = 0; c < staArray.length; c++) {
                    for (int x = 0; x < staArray[0].length; x++) {
                        for (int y = 0; y < staArray[0][0].length; y++) {
                            for (int t = 0; t < staArray[0][0][0].length; t++) {
                                if (covarianceFT[x][y][t][0] > cutoff ||
                                    covarianceFT[x][y][t][0] < -cutoff) {
                                    staArray[c][x][y][t][0] /= covarianceFT[x][y][t][0] *
                                        1000;
                                    //imaginary part of covariance FT is zero.
                                    staArray[c][x][y][t][1] /= covarianceFT[x][y][t][0] *
                                        1000;
                                } else {
                                    staArray[c][x][y][t][0] = 0.0;
                                    staArray[c][x][y][t][1] = 0.0;

                                }
                            }
                        }
                    }
                }

//staArray[0][0][0][0][0] = staArray[1][0][0][0][0] = staArray[2][0][0][0][0] = 0.0;
                transformT(qT, staArray, 1.0);
                transformX(qX, staArray, 1.0);
                transformY(qY, staArray, 1.0);

                for (int t = 0; t < staDepth; t++) {
                    STAFrame frame = (STAFrame) sta.getFrame(t);
                    float[] colors = frame.getBuffer();
                    for (int x = 0; x < width; x++) {
                        for (int y = 0; y < height; y++) {
                            for (int c = 0; c < 3; c++) {
                                colors[3 * (y * width + x) +
                                    c] = (float) staArray[c][x][y][t][0];

                            }
                        }
                    }

                    frame.setBuffer(colors);
                    sta.setFrame(t, frame);
                }

                staFile.addSTA(idList[i], sta);

            }
            staFile.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }


    //[x][y][t][Re or Im]
    public void ftSTAs(String fileNameRoot, STAFile pinkSTAFile) {

        try {
            STAFile staFile = new STAFile(
                fileNameRoot + ".sta", 10000,
                width, height, staDepth, staOffset,
                stixelWidth, stixelHeight, refreshTime);

            int[] idList = pinkSTAFile.getIDList();
            for (int i = 0; i < idList.length; i++) {
                System.out.println("Projecting:  " + idList[i]);
                STA sta = pinkSTAFile.getSTA(idList[i]);
                double[][][][][] staArray = new double[3][width][height][staDepth][2];

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

//            transformT(ftT, staArray, 1.0);
//            transformX(ftX, staArray, 1.0);
//            transformY(ftY, staArray, 1.0);

//                for (int c = 0; c < staArray.length; c++) {
//                    for (int x = 0; x < staArray[0].length; x++) {
//                        for (int y = 0; y < staArray[0][0].length; y++) {
//                            for (int t = 0; t < staArray[0][0][0].length; t++) {
//                                if (covarianceFT[x][y][t][0] > cutoff ||
//                                    covarianceFT[x][y][t][0] < -cutoff) {
//                                    staArray[c][x][y][t][0] /= covarianceFT[x][y][t][0];
//                                    //imaginary part of covariance FT is zero.
//                                    staArray[c][x][y][t][1] /= covarianceFT[x][y][t][0];
//                                } else {
//                                    staArray[c][x][y][t][0] = 0.0;
//                                    staArray[c][x][y][t][1] = 0.0;
//
//                                }
//                            }
//                        }
//                    }
//                }

//staArray[0][0][0][0][0] = staArray[1][0][0][0][0] = staArray[2][0][0][0][0] = 0.0;
//                transformT(qT, staArray, 1.0);
//                transformX(qX, staArray, 1.0);
//                transformY(qY, staArray, 1.0);
//
                for (int t = 0; t < staDepth; t++) {
                    STAFrame frame = (STAFrame) sta.getFrame(t);
                    float[] colors = frame.getBuffer();
                    for (int x = 0; x < width; x++) {
                        for (int y = 0; y < height; y++) {
                            for (int c = 0; c < 3; c++) {

                                Math.log(colors[3 * (y * width + x) +
                                         c] =
                                             (float) Math.sqrt(staArray[c][x][y][t][0] *
                                    staArray[c][x][y][t][0] +
                                    staArray[c][x][y][t][1] *
                                    staArray[c][x][y][t][1]));

                            }
                        }
                    }

                    frame.setBuffer(colors);
                    sta.setFrame(t, frame);
                }

                staFile.addSTA(idList[i], sta);

            }
            staFile.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }


    //Multiplies two complex matrices, which is functionality that is not available in algebra.mult
    public DenseDoubleMatrix2D[] matrixMultiply(DenseDoubleMatrix2D[] A,
                                                DenseDoubleMatrix2D[] B) {
        DenseDoubleMatrix2D tempA;
        DenseDoubleMatrix2D tempB;

        DenseDoubleMatrix2D[] result = new DenseDoubleMatrix2D[2];
        result[0] = new DenseDoubleMatrix2D(A[0].rows(), A[0].rows());
        result[1] = new DenseDoubleMatrix2D(A[0].rows(), A[0].rows());
        tempA = (DenseDoubleMatrix2D) algebra.mult(A[0], B[0]);
        tempB = (DenseDoubleMatrix2D) algebra.mult(A[1], B[1]);
        for (int i = 0; i < A[0].rows(); i++) {
            for (int j = 0; j < A[0].rows(); j++) {
                result[0].set(i, j, tempA.get(i, j) - tempB.get(i, j));
            }
        }
        tempA = (DenseDoubleMatrix2D) algebra.mult(A[0], B[1]);
        tempB = (DenseDoubleMatrix2D) algebra.mult(A[1], B[0]);
        for (int i = 0; i < A[0].rows(); i++) {
            for (int j = 0; j < A[0].rows(); j++) {
                result[1].set(i, j, tempA.get(i, j) + tempB.get(i, j));
            }
        }
        return result;
    }


    //Multiplies a complex matrix by a vector, which is functionality that is not available in algebra.mult
    public DenseDoubleMatrix1D[] matrixMultiply(DenseDoubleMatrix2D[] A,
                                                DenseDoubleMatrix1D[] V) {
        DenseDoubleMatrix1D tempA;
        DenseDoubleMatrix1D tempB;

        DenseDoubleMatrix1D[] result = new DenseDoubleMatrix1D[2];
        result[0] = new DenseDoubleMatrix1D(V[0].size());
        result[1] = new DenseDoubleMatrix1D(V[0].size());
        tempA = (DenseDoubleMatrix1D) algebra.mult(A[0], V[0]);
        tempB = (DenseDoubleMatrix1D) algebra.mult(A[1], V[1]);

        for (int i = 0; i < V[0].size(); i++) {
            result[0].set(i, tempA.get(i) - tempB.get(i));
        }

        tempA = (DenseDoubleMatrix1D) algebra.mult(A[0], V[1]);
        tempB = (DenseDoubleMatrix1D) algebra.mult(A[1], V[0]);
        for (int i = 0; i < V[0].size(); i++) {
            result[1].set(i, tempA.get(i) + tempB.get(i));
        }

        return result;
    }


    public void displayMovieSpectrum(double[][][][] powerSpectrum) {

        double average = 0.0;
        double[] xFT = new double[2 * width - 2];
        for (int x = 0; x < 2 * width - 2; x++) {
            xFT[x] = Math.sqrt(powerSpectrum[x][0][0][0] * powerSpectrum[x][0][0][0] +
                               powerSpectrum[x][0][0][1] * powerSpectrum[x][0][0][1]); // imaginary part is zero
            average += powerSpectrum[x][0][0][0];
        }
        average /= width;
        for (int x = 0; x < width; x++) {
            xFT[x] /= average;
        }
        ScatterPlot xScatter = new ScatterPlot("");
        ScatterPlot ref = new ScatterPlot("");
        for (int i = 1; i < xFT.length / 2; i++) {
//            double omega = Math.PI * 2 * (double) i / (double) powerSpectrum.length;
            double omega = (double) i / (2 * width - 2) / .025;
            xScatter.add(Math.log(omega) / Math.log(10.0),
                         Math.log(xFT[i]) / Math.log(10.0));

            ref.add(Math.log(omega) / Math.log(10.0),
                    Math.log(xFT[1] / (double) omega / (2 * width - 2) / .025) /
                    Math.log(10.0));
        }

        PlotPanel panel = new PlotPanel();
        panel.addData(xScatter, new ScatterPlotStyle(
            SymbolType.NONE, 0, Color.black, true, Color.black, 1));

        panel.addData(ref, new ScatterPlotStyle(
            SymbolType.NONE, 0, Color.red, true, Color.red, 1));
        panel.setLabels("log(frequency(mm^-1))", "log(power/averagePower)");

        panel.autoscale();
//        panel.setRange(0, 1.5, panel.getRange()[2], panel.getRange()[3]);
        PlotUtil.showData("X Power Spectrum", panel);

        average = 0.0;
        double[] tFT = new double[staDepth - 1];

        for (int t = 0; t < staDepth - 1; t++) {
            tFT[t] = Math.sqrt(powerSpectrum[0][0][t][0] * powerSpectrum[0][0][t][0] +
                               powerSpectrum[0][0][t][1] * powerSpectrum[0][0][t][1]); // imaginary part is zero
            average += powerSpectrum[0][0][t][0];
        }
        average /= staDepth - 1;
        for (int t = 0; t < staDepth - 1; t++) {
            tFT[t] /= average;
        }

        ScatterPlot tScatter = new ScatterPlot("");
        ref = new ScatterPlot("");
        for (int i = 1; i < staDepth - 1; i++) {
            double omega = (double) i / (2 * staDepth - 2) / .00833;
            tScatter.add(Math.log(omega) / Math.log(10.0),
                         Math.log(tFT[i]) / Math.log(10.0));

            ref.add(Math.log(omega) / Math.log(10.0),
                    Math.log(tFT[1] / (double) omega / (2 * staDepth - 2) / .0083) /
                    Math.log(10.0));
        }

        panel = new PlotPanel();
        panel.addData(tScatter, new ScatterPlotStyle(
            SymbolType.NONE, 0, Color.black, true, Color.black, 1));

        panel.addData(ref, new ScatterPlotStyle(
            SymbolType.NONE, 0, Color.red, true, Color.red, 1));
        panel.setLabels("log(frequency(Hz))", "log(power/averagePower)");

        panel.autoscale();
//        panel.setRange(0, 1.5, panel.getRange()[2], panel.getRange()[3]);
        PlotUtil.showData("Temporal Power Spectrum", panel);

        double[][] ftXY = new double[2 * height - 2][2 * height - 2];
        for (int x = 0; x < 2 * height - 2; x++) {
            for (int y = 0; y < 2 * height - 2; y++) {
                ftXY[x][y] = Math.log(powerSpectrum[2 * ( (x + 64)) %
                                      (2 * width - 2)][ (y + 64) % (2 * height - 2)][0][0]) /
                             Math.log(10.0);
//                ftXY[x][y] = Math.log(covariance[x][y][0][0])/Math.log(10.0);
//                 ftXY[x][y] = covariance[x][y][0][0];

            }
        }

        PlotUtil.showData("Spatial Power Spectrum",
                          new DoubleHistogram2D("", -20, 20, -20, 20, ftXY),
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
    
    public static void testWhitenPinkSTAsCalculation() throws IOException {
        WhitenPinkSTAsCalculation whiten = new WhitenPinkSTAsCalculation();
        whiten.algebra = new ComplexAlgebra();
//        String fileNameRoot = "d:\\data\\2004-03-23-0\\data006\\data006";
        String fileNameRoot = "g:\\data\\2004-08-06-0\\data003\\data003";
        STAFile pinkSTAFile = new STAFile(fileNameRoot + ".staP");
        whiten.movie = new RawMovie(fileNameRoot + ".rawMovie",
                                    fileNameRoot + ".globals");
        STA sta0 = pinkSTAFile.getSTA(pinkSTAFile.getIDList()[0]);
        whiten.width = sta0.getWidth();
        whiten.height = sta0.getHeight();
        whiten.staDepth = pinkSTAFile.getSTADepth();
        whiten.staOffset = pinkSTAFile.getSTAOffset();
        whiten.refreshTime = sta0.getRefreshTime();
        whiten.stixelWidth = pinkSTAFile.getStixelWidth();
        whiten.stixelHeight = pinkSTAFile.getStixelHeight();
        whiten.generateEigenvectorMatrices();

// whiten.calculateFullMovieCovariance(fileNameRoot, pinkSTAFile);
        double[][][][] covariance = whiten.quickCalculateFullMovieCovariance(fileNameRoot);
//        double[][][][] covariance = whiten.powerLawCalculateFullMovieCovariance(fileNameRoot);
// double[][][][] covariance = whiten.read4DFile(fileNameRoot + ".lam2");

        //convert covariance to power spectrum
        System.out.println("Transforming covariance to power spectrum");
        whiten.transformX(whiten.qXDagger, covariance, 1);
        whiten.transformY(whiten.qYDagger, covariance, 1);
        whiten.transformT(whiten.qTDagger, covariance, 1);

        whiten.displayMovieSpectrum(covariance);
//    whiten.ftSTAs(fileNameRoot, new STAFile(fileNameRoot + ".sta001"));
// whiten.writeEigenFile(covariance, fileNameRoot + ".ft");
//        double[][][][] covariance = whiten.read4DFile(fileNameRoot + ".ft");
        /*      double max = 0.0;
              for (int x = 0; x < 2 * whiten.width - 2; x++) {
                  for (int y = 0; y < 2 * whiten.height - 2; y++) {
                      for (int t = 0; t < 2 * whiten.staDepth - 2; t++) {
                          for (int p = 0; p < 2; p++) {
                              if (Math.abs(covariance[x][y][t][p]) > max)
                                  max = Math.abs(covariance[x][y][t][p]);
                          }
                      }
                  }
              }
              whiten.whitenSTAs(fileNameRoot, pinkSTAFile, covariance, .001  * max);
         */
//        whiten.ftSTAs(fileNameRoot, pinkSTAFile);

    }


}
