package edu.ucsc.neurobiology.vision.util;

import edu.ucsc.neurobiology.vision.math.*;


/**
 * Uses cubic spline interpolation to find the exact spike time and then resamples the
 * spike so that the peak of the spike exactly corresponds with a sample point.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class SpikeAligner {
    private final int nPoints, nlPoints;
    private int[][] adjacent;
    private double[] spikeData, /*uniformX,*/ alignedSpike;
    private UniformSpline spline;
    private final double minimizationError;


    public SpikeAligner(int nPoints, int nlPoints, int[][] adjacent, double minimizationError) {
        this.nPoints = nPoints;
        this.nlPoints = nlPoints;
        this.adjacent = adjacent;
        this.minimizationError = minimizationError;
        spikeData = new double[nPoints];
        spline = new UniformSpline(nPoints);

        int maxAdjacents = -1;
        for (int electrode = 1; electrode < adjacent.length; electrode++) {
            if (adjacent[electrode].length > maxAdjacents) {
                maxAdjacents = adjacent[electrode].length;
            }
        }
        alignedSpike = new double[nPoints * maxAdjacents];
    }

    // ToDo: 2011-01: If doAlignment is false, why do we do all this work?
    public double[] align(int electrode, short[] spike, boolean doAlignment) throws CannotEvaluateException {
        double x0 = 0;

        for (int el = 0; el < adjacent[electrode].length; el++) {
            for (int i = 0, index1 = el * nPoints; i < nPoints; i++) {
                spikeData[i] = spike[index1 + i];
            }

            spline.reSpline( /*uniformX, */spikeData);

            if (el == 0) {
                if (doAlignment) {
                x0 = FunctionMinimization.brentParabolic(
                    spline, nlPoints - 1, nlPoints + 1, minimizationError);
                } else {
                    x0 = nlPoints;
                }
            }

//            double n = x0 + 2 - nlPoints;
            double n = x0 - nlPoints;
            for (int i = 0, index2 = el * (nPoints - 2); i < nPoints - 2; i++) {
                alignedSpike[index2 + i] = (float) spline.getValueAt(i + n);
            }
        }

        return alignedSpike;
    }
    
    public double[] align(int electrode, int[] spike, boolean doAlignment) throws CannotEvaluateException {
        double x0 = 0;

        for (int el = 0; el < adjacent[electrode].length; el++) {
            for (int i = 0, index1 = el * nPoints; i < nPoints; i++) {
                spikeData[i] = spike[index1 + i];
            }

            spline.reSpline( /*uniformX, */spikeData);

            if (el == 0) {
                if(doAlignment) {
                x0 = FunctionMinimization.brentParabolic(
                    spline, nlPoints - 1, nlPoints + 1, minimizationError);
                } else {
                    x0 = nlPoints;
                }
            }

//            double n = x0 + 2 - nlPoints;
            double n = x0 - nlPoints;
            for (int i = 0, index2 = el * (nPoints - 2); i < nPoints - 2; i++) {
                alignedSpike[index2 + i] = (float) spline.getValueAt(i + n);
            }
        }

        return alignedSpike;
    }
    
}