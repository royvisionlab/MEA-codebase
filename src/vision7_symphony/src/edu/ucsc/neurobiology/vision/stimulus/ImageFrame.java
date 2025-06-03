package edu.ucsc.neurobiology.vision.stimulus;

import java.awt.Point;

import edu.ucsc.neurobiology.vision.math.FitFailedException;
import edu.ucsc.neurobiology.vision.math.MathUtil;
import edu.ucsc.neurobiology.vision.math.fitting.Fitter;
import edu.ucsc.neurobiology.vision.math.fitting.Gaussian2DFunction;


/**
 * This abstract class specifies a single image frame.
 * It is the parent of any frames used in Vision. It defines only basic methods
 * shared by any kind of frame. The way the frame actually stores its colors
 * and the way the frame updates those colors are not specified here by in implementing
 * classes.
 *
 * @author Charles A. Loomis, University of California, Santa Cruz
 *         Dumitru Petrusca, University of California, Santa Cruz
 */
public abstract class ImageFrame {

    /** The width of the frame in pixels */
    protected int width;

    /** The height of the frame in pixels */
    protected int height;

    /** The size of a stixel in microns*/
    protected double stixelWidth;
    protected double stixelHeight;


    public ImageFrame(int width, int height, double pixelWidth, double pixelHeight) {
        this.width = width;
        this.height = height;
        this.stixelWidth = pixelWidth;
        this.stixelHeight = pixelHeight;
    }
    
    
    /**
     * This method returns the width of the Frame in pixels.
     * @return The width of the frame
     */
    public int getWidth() {
        return width;
    }


    /**
     * This method returns the height of the Frame in pixels.
     * @return The height of the frame
     */
    public int getHeight() {
        return height;
    }


    /**
     * Returns the pixel width in microns.
     */
    public double getStixelWidth() {
        return stixelWidth;
    }
    
    /**
     * Returns the pixel height in microns.
     */
    public double getStixelHeight() {
        return stixelHeight;
    }

    public void setPixelSize(double stixelWidth, double stixelHeight) {
        this.stixelWidth = stixelWidth;
        this.stixelHeight = stixelHeight;
    }
    
    public abstract ImageFrame downBin(int downFactor);

    /**
     * This method returns the color components for all pixels in the Frame.
     * The first and second indexes in the array are the x and y coordinates.
     * The third index is the color component index in the order: red, green, blue.
     *
     * @return The internal buffer of this frame
     */
    public abstract float[] getBuffer();


    public abstract float[] getPixel(int x, int y, float[] pixel);


    public abstract float getPixel(int x, int y, int c);

    
    public abstract float[] getPixelError(int x, int y, float[] error);


    public abstract float getPixelError(int x, int y, int c);


    public float getPixelSignificance(int x, int y, int c) {
        return Math.abs(getPixel(x, y, c) / getPixelError(x, y, c));
    }


    public float getPixelSignificance(int x, int y) {
        float[] c = new float[3];
        float[] e = new float[3];
        getPixel(x, y, c);
        getPixelError(x, y, e);
        float significance = (float) Math.sqrt(
            (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]) /
            (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]));
        return significance;
    }


    public Object clone() {
        try {
            return super.clone();
        } catch (Exception e) {
            return null;
        }
    }


    public float getMaxAbsColor(Point p, int colorComponent) {
        if (p == null) {
            throw new IllegalArgumentException("p cannon be null");
        }

        float c, max = Float.NEGATIVE_INFINITY;

        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                c = Math.abs(getPixel(i, j, colorComponent));
                if (c > max) {
                    max = c;
                    p.x = i;
                    p.y = j;
                }
            }
        }

        return max;
    }


    public float getMaxAbsColor(Point p) {
        if (p == null) {
            throw new IllegalArgumentException("p cannot be null");
        }

        float max = Float.NEGATIVE_INFINITY;
        float[] c = new float[3];

        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                getPixel(i, j, c);
                float color = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
                if (color > max) {
                    max = color;
                    p.x = i;
                    p.y = j;
                }
            }
        }

        return max;
    }


    /**
     * Finds the average amplitude squared from the largest nPixels in a frame.
     * @author Matthew Grivich
     */
    public float getPeakMagnitudeSquared(int nPixels) {
        float[] c = new float[3];
        float minimumHighValue = Float.NEGATIVE_INFINITY;
        int index = 0;
        float[] magnitudesSquared = new float[nPixels];
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                getPixel(i, j, c);
                float color = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
                if (color > minimumHighValue) {
                    magnitudesSquared[index] = color;
                    minimumHighValue = color;
                    for (int k = 0; k < magnitudesSquared.length; k++) {
                        if (magnitudesSquared[k] < minimumHighValue) {
                            minimumHighValue = magnitudesSquared[k];
                            index = k;
                        }
                    }
                }

            }
        }

        return (float) MathUtil.mean(magnitudesSquared);

        /*   float minimumHighValue = Float.NEGATIVE_INFINITY;
         float[] c = new float[3];
         int index = 0;
         float[] magnitudesSquared = new float[nPixels];
         for (int i = 0; i < width; i++) {
             for (int j = 0; j < height; j++) {
                 getPixel(i, j, c);
                 float color = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
                 if (color > minimumHighValue) {
                     magnitudesSquared[index] = color;
                     minimumHighValue = color;


                 }
             }
         }

         return magnitudesSquared[index];*/

    }

    
    /**
     * Convert the frame object into a collection of primitives; ugly.  Also normalizes.
     * @param f
     * @param x
     * @param y
     * @param sigma
     * @param colorIndex
     * @param p
     * @param halfSize
     * @return
     */
    public int primitivize(double[][] x, double[] y, double[] sigma, int colorIndex,
            Point p, int halfSize) {

        int n = 0;

        int x1 = Math.max(p.x - halfSize, 0);
        int x2 = Math.min(p.x + halfSize, width - 1);
        int y1 = Math.max(p.y - halfSize, 0);
        int y2 = Math.min(p.y + halfSize, height - 1);
        double average = 0;
        for (int i = x1; i <= x2; i++) {
            for (int j = y1; j <= y2; j++) {
                // the x and y coordinates
                x[0][n] = i + 0.5;
                x[1][n] = height - j - 1 + 0.5;

                // the value and the sigma
                y[n] = Math.abs(this.getPixel(i, j, colorIndex));
                sigma[n] = this.getPixelError(i, j, colorIndex);

                average += y[n];
                n++;
            }
        }
        average /= y.length;
        for (int i = 0; i < y.length; i++) {
            y[i] -= average;
        }

        return n;
    }


    public Gaussian2DFunction fit(int color) {
        boolean DEBUG = false;
        
        Point maxPoint = new Point();
        this.getMaxAbsColor(maxPoint);
        
        double[][] x = new double[2][width * height];
        double[] y = new double[width * height];
        double[] sigma = new double[width * height];
        Point p = new Point(width / 2, height / 2);
        int nPoints = primitivize(x, y, sigma, color, p, Math.max(width, height));
        
        double[] params = new double[7];
        int maxAbsIndex = MathUtil.maxAbsIndex(y);
        Fitter fitter = new Fitter();
        fitter.setLastChiSquaredVariation(1e-3);
        fitter.setMaxIterations(1000);
        Gaussian2DFunction gBest = null;
        double minChi2 = Double.POSITIVE_INFINITY;
        
        // Fit is highly unstable.  Try a range of sigmas to make getting a good fit likely.
        for (double i = 1; i <= 3; i += 0.5) {
            for (double j = 1; j <= 3; j += 0.5) {
                params[0] = y[maxAbsIndex]; // A
                params[1] = maxPoint.x + 0.5; // x0
                params[2] = height - maxPoint.y - 1 + 0.5; // y0
                params[3] = i; // sigX
                params[4] = j; // sigY
                params[5] = 0; // Theta, radians
                params[6] = 0; //B

                Gaussian2DFunction g2d = new Gaussian2DFunction();
                g2d.setParameters(params);
                try {
                    fitter.fit(g2d, x, y, sigma, nPoints);
                    if (DEBUG) {
                        System.out.println(g2d);
                    }
                    double chi2 = g2d.getChiSquared();
                    if (chi2 < minChi2) {
                        minChi2 = chi2;
                        gBest = g2d;
                    }
                } catch (FitFailedException e) {
                    if (DEBUG) {
                        System.out.println(e.getMessage() + "\n");
                    }
                }
            }
        }
        
        // validity cuts
        if (gBest == null) {
            return null;
        } else {
            gBest.normalizeVariables();
            if (gBest.getX0() > width * 2 ||
                    gBest.getX0() < -width ||
                    gBest.getY0() > height * 2 ||
                    gBest.getY0() < -height ||
                    gBest.getSigmaX() > width * 2 ||
                    gBest.getSigmaX() < -width ||
                    gBest.getSigmaY() > height * 2 ||
                    gBest.getSigmaY() < -height) {

                return null;
            }
        }
        
        return gBest;
    }

    
    /**
     * Should be further folded together with normal fit method if possible.
     * @param color
     * @param x0
     * @param y0
     * @return
     */
    public Gaussian2DFunction improveFit(int color, double x0, double y0) {
        boolean DEBUG = false;

        double[][] x = new double[2][width * height];
        double[] y = new double[width * height];
        double[] sigma = new double[width * height];
        Point p0 = new Point();
        getMaxAbsColor(p0);

        if (DEBUG) {
            System.out.println("coords [" + p0 + "]");
            System.out.println("================================================");
        }

        // prepare the initial parameters
        Point p = new Point(width / 2, height / 2);
        int nPoints = primitivize(x, y, sigma, color, p, Math.max(width, height));
        int maxAbsIndex = MathUtil.maxAbsIndex(y);
        Fitter fitter = new Fitter();
        fitter.setLastChiSquaredVariation(1e-3);
        fitter.setMaxIterations(1000);
        Gaussian2DFunction gBest = null;
        double minChi2 = Double.POSITIVE_INFINITY;

        for (double i = 1; i <= 5; i += 0.5) {
            for (double j = 1; j <= 5; j += 0.5) {
                double sigX = i;
                double sigY = j;

                x0 = (x0 > 0) ? x0 : p0.x + 0.5;
                y0 = (y0 > 0) ? y0 : height - p0.y - 1 + 0.5;

                Gaussian2DFunction g2d = new Gaussian2DFunction(
                        y[maxAbsIndex], x0, y0, sigX, sigY, 0, 0);
                g2d.setParameterState(1, false);
                g2d.setParameterState(2, false);

                try {
                    fitter.fit(g2d, x, y, sigma, nPoints);
                    if (DEBUG) {
                        System.out.println(g2d);
                    }
                    double chi2 = g2d.getChiSquared();
                    if (chi2 < minChi2) {
                        minChi2 = chi2;
                        gBest = g2d;
                    }
                } catch (FitFailedException e) {
                    if (DEBUG) {
                        System.out.println(e.getMessage() + "\n");
                    }
                }
            }
        }

        // validity cuts
        if (gBest == null) {
            return null;
        } else {
            gBest.normalizeVariables();
            if (gBest.getX0() > width * 2 || gBest.getX0() < -width ||
                    gBest.getY0() > height * 2 || gBest.getY0() < -height) {
                return null;
            }
        }
        return gBest;
    }
   
}
