package edu.ucsc.neurobiology.vision.test;

import java.awt.*;

import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.math.fitting.*;
import edu.ucsc.neurobiology.vision.stimulus.*;


/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2002</p>
 *
 * <p>Company: </p>
 *
 * @author not attributable
 * @version 1.0
 */
public class STADOG {
    /**
     * Fits the main frame of the STA to a general two-dimeansional gaussian.
     *
     * @param sta the sta to be fitted
     * @param color the collor component to be used (0-red, 1-green, 2-blue)
     * @return an instance of <i>Gaussian2DFunction</i> containing the parameters of
     * the gaussian
     */
    public static DOG2DFunction fitDOG(STA sta, int color) {
        boolean DEBUG = true;

        final int w = sta.getWidth(), h = sta.getHeight();
        double[][] x = new double[2][w * h];
        double[] y = new double[w * h];
        double[] sigma = new double[w * h];
        int[] maxFrameParams = sta.getMainFrameParams();
        int frameIndex = maxFrameParams[0];
        int x0 = maxFrameParams[1];
        int y0 = maxFrameParams[2];

        if (DEBUG) {
            System.out.println(
                "STA " + ", frame " + frameIndex +
                ", coords [" + x0 + ":" + y0 + "]");
            System.out.println("================================================");
        }

        // prepare the initial parameters
        ImageFrame f = sta.getFrame(frameIndex);
        Point p = new Point();
        f.getMaxAbsColor(p, color);
        float A = f.getPixel(p.x, p.y, color);

        int nPoints = convertFrameNew(
            f, x, y, sigma, color, new Point(w / 2, h / 2), Math.max(w, h));

        Gaussian2DFunction gf = sta.fit(color);

        double[] params = new double[8];
        params[0] = gf.getX0(); // x0
        params[1] = gf.getY0(); // y0
        params[2] = gf.getSigmaX(); // sigX1
        params[3] = gf.getSigmaY(); // sigY1
        params[4] = 4; // r
        params[5] = gf.getAmplitude(); // A1
        params[6] = -0.1 * gf.getAmplitude(); // A2
        params[7] = gf.getTheta(); // theta
        DOG2DFunction g2d = new DOG2DFunction();
        g2d.setParameters(params);
        if (DEBUG) {
            System.out.println(g2d);
        }
        Fitter fitter = new Fitter();
        fitter.setLastChiSquaredVariation(1e-4);
        fitter.setMaxIterations(10000);
        try {
            fitter.fit(g2d, x, y, sigma, nPoints);

            if (DEBUG) {
                System.out.println(g2d);
            }
        } catch (FitFailedException e) {
            if (DEBUG) {
                System.out.println(e.getMessage() + "\n");
            }
            return null;
        }

        /*
                // validity cuts
                g2d.normalizeVariables();
                if (g2d.getX0() > sta.getWidth() * 2 ||
                    g2d.getX0() < -sta.getWidth() ||
                    g2d.getY0() > sta.getHeight() * 2 ||
                    g2d.getY0() < -sta.getHeight()
                    ) {
                    return null;
                }
         */
        return g2d;
    }


    private static int convertFrameNew(
        ImageFrame f, double[][] x, double[] y, double[] sigma, int colorIndex,
        Point p, int halfSize) {

        final int w = f.getWidth(), h = f.getHeight();
        int n = 0;

        int x1 = Math.max(p.x - halfSize, 0);
        int x2 = Math.min(p.x + halfSize, w - 1);
        int y1 = Math.max(p.y - halfSize, 0);
        int y2 = Math.min(p.y + halfSize, h - 1);
        for (int i = x1; i <= x2; i++) {
            for (int j = y1; j <= y2; j++) {
                // the x and y coordinates
                x[0][n] = i + 0.5;
                x[1][n] = h - j - 1 + 0.5;
                // the value and the the sigma
                y[n] = f.getPixel(i, j, colorIndex);
                sigma[n] = f.getPixelError(i, j, colorIndex);

                n++;
            }
        }
        return n;
    }


}
