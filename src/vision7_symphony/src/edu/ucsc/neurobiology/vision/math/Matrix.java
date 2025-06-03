package edu.ucsc.neurobiology.vision.math;

import static java.lang.Math.*;

import cern.colt.matrix.*;
import cern.colt.matrix.impl.*;
import cern.colt.matrix.linalg.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * A simple matrix implementation based on the cern.colt library.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class Matrix {
    private double[][] m;
    private final int nRows, nCols;


    public Matrix(int nRows, int nCols) {
        this.nRows = nRows;
        this.nCols = nCols;
        m = new double[nRows][nCols];
    }


    public Matrix(double[][] m) {
        this.m = m;
        this.nRows = m.length;
        this.nCols = m[0].length;
    }


    public void reset() {
        for (int i = 0; i < nRows; i++) {
            for (int j = 0; j < nCols; j++) {
                this.m[i][j] = 0;
            }
        }
    }


    public double get(int i, int j) {
        return m[i][j];
    }


    public Matrix add(Matrix newMatrix) {
        if (newMatrix.nCols != this.nCols || newMatrix.nRows != this.nRows) {
            throw new IllegalArgumentException(
                "(" + this.nRows + " x " + this.nCols + ")" + " add " +
                "(" + newMatrix.nRows + " x " + newMatrix.nCols + ")"
                );
        }

        for (int i = 0; i < nRows; i++) {
            for (int j = 0; j < nCols; j++) {
                this.m[i][j] += newMatrix.m[i][j];
            }
        }

        return this;
    }


    public Matrix subtract(Matrix newMatrix) {
        if (newMatrix.nCols != this.nCols || newMatrix.nRows != this.nRows) {
            throw new IllegalArgumentException(
                "(" + this.nRows + " x " + this.nCols + ")" + " subtract " +
                "(" + newMatrix.nRows + " x " + newMatrix.nCols + ")"
                );
        }

        for (int i = 0; i < nRows; i++) {
            for (int j = 0; j < nCols; j++) {
                this.m[i][j] -= newMatrix.m[i][j];
            }
        }

        return this;
    }


    public Matrix scale(double sf) {
        for (int i = 0; i < nRows; i++) {
            for (int j = 0; j < nCols; j++) {
                this.m[i][j] *= sf;
            }
        }

        return this;
    }


    public double modulus() {
        return new LUDecomposition(new DenseDoubleMatrix2D(m)).det();
    }


    public Matrix inverse() {
        DoubleMatrix2D inv = Algebra.DEFAULT.inverse(new DenseDoubleMatrix2D(m));
        return new Matrix(inv.toArray());
    }


    public Matrix getAverageColumns() {
        Matrix average = new Matrix(1, nCols);

        for (int j = 0; j < nCols; j++) {
            average.m[0][j] = 0;
            for (int i = 0; i < nRows; i++) {
                average.m[0][j] += this.m[i][j];
            }
            average.m[0][j] /= nRows;
        }

        return average;
    }


    public Matrix getTranspose() {
        Matrix transpose = new Matrix(nCols, nRows);

        for (int i = 0; i < nRows; i++) {
            for (int j = 0; j < nCols; j++) {
                transpose.m[j][i] = this.m[i][j];
            }
        }

        return transpose;
    }


    public static Matrix multiply(Matrix m1, Matrix m2) {
        if (m1.nCols != m2.nRows) {
            throw new IllegalArgumentException("m1.nCols != m2.nRows");
        }

        Matrix m = new Matrix(m1.nRows, m2.nCols);

        for (int i = 0; i < m.nRows; i++) {
            for (int j = 0; j < m.nCols; j++) {
                m.m[i][j] = 0;
                for (int k = 0; k < m1.nCols; k++) {
                    m.m[i][j] += m1.m[i][k] * m2.m[k][j];
                }
            }
        }

        return m;
    }


    public void print(int precision) {
        System.out.println(nRows + " x " + nCols);

        for (int i = 0; i < nRows; i++) {
            for (int j = 0; j < nCols; j++) {
                System.out.print(StringUtil.format(this.m[i][j], precision) + "\t");
            }
            System.out.println();
        }
    }


    public double[] getRow(int rowIndex) {
        double[] row = new double[nCols];

        for (int i = 0; i < nCols; i++) {
            row[i] = m[rowIndex][i];
        }

        return row;
    }


    public static double err = 1e-5;
    static public boolean correct(ParametricEllipse ei, ParametricEllipse e) {
        boolean correct = true;
        for (int i = 0; i < 5; i++) {
            double theta;
            if (i == 4) {
                theta = PI / 6;
            } else {
                theta = i * PI / 2;
            }
            double[] pi = ei.getPointFor(theta);
            double[] p = e.getPointFor(theta);
            if (Math.abs(pi[0] - p[0]) > err && Math.abs(pi[1] - p[1]) > err) {
                correct = false;
                break;
            }
        }
        return correct;
    }

}
