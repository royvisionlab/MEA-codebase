package edu.ucsc.neurobiology.vision.math;

import static java.lang.Math.*;

import cern.colt.matrix.impl.*;
import cern.colt.matrix.linalg.*;


/**
 *
 * @author Matthew Grivich, University of California, Santa Cruz
 */
@SuppressWarnings("serial")
public class ComplexAlgebra
    extends Algebra {


    public ComplexAlgebra() {

    }


    //Multiplies two complex matrices, which is functionality that is not available in algebra.mult
    public DenseDoubleMatrix2D[] matrixMultiply(DenseDoubleMatrix2D[] A,
                                                DenseDoubleMatrix2D[] B) {
        DenseDoubleMatrix2D tempA;
        DenseDoubleMatrix2D tempB;

        DenseDoubleMatrix2D[] result = new DenseDoubleMatrix2D[2];
        result[0] = new DenseDoubleMatrix2D(A[0].rows(), A[0].rows());
        result[1] = new DenseDoubleMatrix2D(A[0].rows(), A[0].rows());
        tempA = (DenseDoubleMatrix2D)this.mult(A[0], B[0]);
        tempB = (DenseDoubleMatrix2D)this.mult(A[1], B[1]);
        for (int i = 0; i < A[0].rows(); i++) {
            for (int j = 0; j < A[0].rows(); j++) {
                result[0].set(i, j, tempA.get(i, j) - tempB.get(i, j));
            }
        }
        tempA = (DenseDoubleMatrix2D)this.mult(A[0], B[1]);
        tempB = (DenseDoubleMatrix2D)this.mult(A[1], B[0]);
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
        tempA = (DenseDoubleMatrix1D)this.mult(A[0], V[0]);
        tempB = (DenseDoubleMatrix1D)this.mult(A[1], V[1]);

        for (int i = 0; i < V[0].size(); i++) {
            result[0].set(i, tempA.get(i) - tempB.get(i));
        }

        tempA = (DenseDoubleMatrix1D)this.mult(A[0], V[1]);
        tempB = (DenseDoubleMatrix1D)this.mult(A[1], V[0]);
        for (int i = 0; i < V[0].size(); i++) {
            result[1].set(i, tempA.get(i) + tempB.get(i));
        }

        return result;
    }


//gives the complex matrix that can be used to fourier transform a vector
    public DenseDoubleMatrix2D[] generateFTMatrix(int N, double sign,
                                                  double normalization) {

        double[][][] matrix = new double[2][N][N];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {

                matrix[0][j][i] = Math.cos(sign * (2 * PI * (double) i * (double) j) /
                                           (double) N) /
                                  normalization;
                matrix[1][j][i] = Math.sin(sign * (2 * PI * (double) i * (double) j) /
                                           (double) N) /
                                  normalization;

            }
        }
        DenseDoubleMatrix2D[] denseMatrix = new DenseDoubleMatrix2D[2];
        denseMatrix[0] = new DenseDoubleMatrix2D(matrix[0]);
        denseMatrix[1] = new DenseDoubleMatrix2D(matrix[1]);
        return denseMatrix;
    }


    //returns a complex vector, given a 3 index, complex matrix.
    //index specifies what kind of vector to get, a and b specifies where to get it from.
    public DenseDoubleMatrix1D[] matrixSubset(double[][][][] matrix, int index, int a,
                                              int b) {
        DenseDoubleMatrix1D[] subset = new DenseDoubleMatrix1D[2];
        if (index == 0) {
            subset[0] = new DenseDoubleMatrix1D(matrix.length);
            subset[1] = new DenseDoubleMatrix1D(matrix.length);
            for (int i = 0; i < matrix.length; i++) {
                subset[0].set(i, matrix[i][a][b][0]);
                subset[1].set(i, matrix[i][a][b][1]);
            }
        } else if (index == 1) {
            subset[0] = new DenseDoubleMatrix1D(matrix[0].length);
            subset[1] = new DenseDoubleMatrix1D(matrix[0].length);
            for (int i = 0; i < matrix.length; i++) {
                subset[0].set(i, matrix[a][i][b][0]);
                subset[1].set(i, matrix[a][i][b][1]);
            }

        } else if (index == 2) {
            subset[0] = new DenseDoubleMatrix1D(matrix[0][0].length);
            subset[1] = new DenseDoubleMatrix1D(matrix[0][0].length);
            for (int i = 0; i < matrix.length; i++) {
                subset[0].set(i, matrix[a][b][i][0]);
                subset[1].set(i, matrix[a][b][i][1]);
            }

        }
        return subset;
    }


    //multiplies the time component of sta by matrix and multiplier
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


    //multiplies the X component of sta by matrix and multiplier
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


    //multiplies the Y component of sta by matrix and multiplier
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


}
