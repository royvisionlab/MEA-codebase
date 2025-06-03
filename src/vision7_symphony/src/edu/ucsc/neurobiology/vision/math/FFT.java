package edu.ucsc.neurobiology.vision.math;

/**
 * An implemnetation of the Fast Fourier Transform.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public final class FFT {
    private static final double log2 = Math.log(2);


    private FFT() {

    }


    private static final int bitrev(int j1, int nu) {
        int k = 0;
        for (int i = 1; i <= nu; i++) {
            k = (k - (j1 >> 1)) * 2 + j1;
            j1 >>= 1;
        }
        return k;
    }


    public static final void fft2(double[][] xre, double[][] xim, int sign) {
        double swap;
        int nx = xre.length;
        int ny = xre[0].length;
        double[] re = new double[nx];
        double[] im = new double[nx];

        // do the row-wise 1D fft
        for (int i = 0; i < nx; i++) {
            fft(xre[i], xim[i], sign);
        }

        // do the column-wise 1D fft
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                re[i] = xre[i][j];
                im[i] = xim[i][j];
            }
            fft(re, im, sign);
            for (int i = 0; i < nx; i++) {
                xre[i][j] = re[i];
                xim[i][j] = im[i];
            }
        }

        // invert
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny / 2; j++) {
                swap = xre[i][j];
                xre[i][j] = xre[i][j + ny / 2];
                xre[i][j + ny / 2] = swap;

                swap = xim[i][j];
                xim[i][j] = xim[i][j + ny / 2];
                xim[i][j + ny / 2] = swap;
            }
        }

        for (int i = 0; i < nx / 2; i++) {
            for (int j = 0; j < ny; j++) {
                swap = xre[i][j];
                xre[i][j] = xre[i + nx / 2][j];
                xre[i + nx / 2][j] = swap;

                swap = xim[i][j];
                xim[i][j] = xim[i + nx / 2][j];
                xim[i + nx / 2][j] = swap;
            }
        }

    }


    public static final void fft(double[] xre, double[] xim, final int sign) {
        if (xre.length != xim.length) {
            throw new ArrayIndexOutOfBoundsException("xre.length != xim.length");
        }

        // assume n is a power of 2
        final int n = xre.length;
        int nu = (int) (Math.log(n) / log2);
        int nHalf = n >> 1;
        int nu1 = nu - 1;
        int k = 0;

        double tr, ti, p, arg, c, s;
        double pre1 = arg = sign * 2 * Math.PI / n;

        for (int l = 1; l <= nu; l++) {
            while (k < n) {
                for (int i = 1; i <= nHalf; i++) {
                    p = bitrev(k >> nu1, nu);
//                    arg = sign * 2 * Math.PI * p / n;
                    arg = pre1 * p;
                    c = Math.cos(arg);
                    s = Math.sin(arg);

                    final int m = k + nHalf;
                    tr = xre[m] * c + xim[m] * s;
                    ti = xim[m] * c - xre[m] * s;
                    xre[m] = xre[k] - tr;
                    xim[m] = xim[k] - ti;
                    xre[k] += tr;
                    xim[k] += ti;
                    k++;
                }
                k += nHalf;
            }
            k = 0;
            nu1--;
            nHalf >>= 1;
        }

        k = 0;
        while (k < n) {
            int r = bitrev(k, nu);
            if (r > k) {
                tr = xre[k];
                ti = xim[k];
                xre[k] = xre[r];
                xim[k] = xim[r];
                xre[r] = tr;
                xim[r] = ti;
            }
            k++;
        }
    }


    public final double[] spectrum(double[] x) {
        // assume n is a power of 2
        int n = x.length;
        int nu = (int) (Math.log(n) / log2);
        int nHalf = n / 2;
        int nu1 = nu - 1;
        double[] xre = new double[n];
        double[] xim = new double[n];
        double tr, ti, p, arg, c, s;
        for (int i = 0; i < n; i++) {
            xre[i] = x[i];
            xim[i] = 0.0f;
        }
        int k = 0;

        for (int l = 1; l <= nu; l++) {
            while (k < n) {
                for (int i = 1; i <= nHalf; i++) {
                    p = bitrev(k >> nu1, nu);
                    arg = 2 * Math.PI * p / n;
                    c = Math.cos(arg);
                    s = Math.sin(arg);
                    tr = xre[k + nHalf] * c + xim[k + nHalf] * s;
                    ti = xim[k + nHalf] * c - xre[k + nHalf] * s;
                    xre[k + nHalf] = xre[k] - tr;
                    xim[k + nHalf] = xim[k] - ti;
                    xre[k] += tr;
                    xim[k] += ti;
                    k++;
                }
                k += nHalf;
            }
            k = 0;
            nu1--;
            nHalf = nHalf / 2;
        }
        k = 0;
        int r;
        while (k < n) {
            r = bitrev(k, nu);
            if (r > k) {
                tr = xre[k];
                ti = xim[k];
                xre[k] = xre[r];
                xim[k] = xim[r];
                xre[r] = tr;
                xim[r] = ti;
            }
            k++;
        }

        double[] mag = new double[nHalf];
        mag[0] = (Math.sqrt(xre[0] * xre[0] + xim[0] * xim[0])) / n;
        for (int i = 1; i < n / 2; i++) {
            mag[i] = 2 * (Math.sqrt(xre[i] * xre[i] + xim[i] * xim[i])) / n;
        }

        return mag;
    }


    public static final double[] spectrum(double[] xre, double[] xim) {
        int n = xre.length;
        int nHalf = n / 2;

        double[] mag = new double[nHalf];
        mag[0] = (Math.sqrt(xre[0] * xre[0] + xim[0] * xim[0])) / n;
        for (int i = 1; i < n / 2; i++) {
            mag[i] = 2 * (Math.sqrt(xre[i] * xre[i] + xim[i] * xim[i])) / n;
        }

        return mag;
    }


    public final double[][] spectrum(double[][] xre, double[][] xim) {
        int nx = xre.length;
        int ny = xre[0].length;
//        int nHalf = n / 2;ectrum

        double[][] spect = new double[nx][ny];
//        mag[0] =  (Math.sqrt(xre[0] * xre[0] + xim[0] * xim[0])) / n;
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                spect[i][j] =
                    2 * (Math.sqrt(xre[i][j] * xre[i][j] + xim[i][j] * xim[i][j])) /
                    (nx * ny);
            }
        }

        return spect;
    }

    /*
        public final double[] phase(double[][] x) {
            int n = x[0].length;
            int nHalf = n / 2;

            double[] phase = new double[nHalf];
//        mag[0] = (Math.sqrt(x[0][0] * x[0][0] + x[1][0] * x[1][0])) / n;
//        for (int i = 1; i < nHalf; i++) {
//            mag[i] = 2 * (Math.sqrt(x[0][i] * x[0][i] + x[1][i] * x[1][i])) / n;
//        }

            return phase;
        }
     */

}
