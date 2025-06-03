package edu.ucsc.neurobiology.vision.math;

import java.util.*;

import java.awt.*;

import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.util.VisionParams;



/**
 * This class implements the Expectation Maximization (EM) algorithm for mixtures of
 * multivariate gaussians.
 * 
 * Clusters are not allowed to rotate.  Dumitru Petrusca and Jon determined that this
 * reduces the tendency of clusters to collapse.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 * @author Matthew Grivich, The Salk Institute
 */
public abstract class ExpectationMaximization {

    //   public static final int ONE_GAUSSIAN_PER_CLUSTER = 0;
    //   public static final int TWO_GAUSSIANS_PER_CLUSTER = 1;

    protected boolean updateMeans = true;
    protected boolean updateVariances = true;
    protected boolean updateWeights = true;

    protected final int nDimensions;
    protected final int maxGaussians;
    protected final float expThreshold = 100;  //used for building e^x lookup table
    protected final int nExpSteps = 100000;  //used for building e^x lookup table
    protected final float stepSize;  //used for building e^x lookup table
    protected final float[] expLookup; //the e^x lookup table

    protected int nGaussians;
//	protected int fullNPoints;
//	protected int truncatedNPoints = -1;
    protected double factor2pi;
    protected float[][] x; // x[dimension][spike], the data points
    protected double[][] mean; // mean[gaussian][dimension], the cluster means
    protected double[][] cjj; // Variances, cjj[gaussian][dimension], the cluster variances
    protected double[] pJ; // pJ[gaussian] Mixing Parameters
    protected float[][] pJX; // pjX[gaussian][spike] probability of each point belonging to each gaussian
    protected double[] variance; //variance[dimension]
    double[][] cjjInverse; //cjjInverse[gaussian][dimension]
    double[] normalizationFactor; //normalizationFactor[gaussian]
    double[] pXJtemp; //pXJtemp[gaussian]
    protected boolean[] variable;  //variable[gaussian] Whether or not a gaussian has collapsed.
    protected int collapsedCount; // check if all gaussians have collapsed

    protected boolean generateReport = false;
    protected boolean valid = false;

    private String identifier = "";  //Identifier for what process is going on.  Used to improve error messages.

    private PointChooser fullPC;  //point chooser with all the spikes
    private PointChooser reducedPC;  //point chooser with a reduced set of spikes


    protected ExpectationMaximization(
            int nDimensions, int maxGaussians, int maxNPoints, float[][] _pJX) {

        this.nDimensions = nDimensions;
        this.maxGaussians = maxGaussians;
        collapsedCount = 0;

        expLookup = new float[nExpSteps];
        stepSize = expThreshold / nExpSteps;
        for (int i = 0; i < nExpSteps; i++) {
            expLookup[i] = (float) Math.exp( -i * stepSize);
        }

        mean = new double[maxGaussians][nDimensions];
        cjj = new double[maxGaussians][nDimensions];
        pJ = new double[maxGaussians];
        variance = new double[nDimensions];
        cjjInverse = new double[maxGaussians][nDimensions];
        normalizationFactor = new double[maxGaussians];
        pXJtemp = new double[maxGaussians];
        variable = new boolean[maxGaussians];


        // memory demanding fields
        if (_pJX == null) {
            pJX = new float[maxGaussians][maxNPoints];
        } else {
            pJX = _pJX;
        }

    }




    public final void setState(
            boolean updateMeans, boolean updateVariances, boolean updateWeights) {

        this.updateMeans = updateMeans;
        this.updateVariances = updateVariances;
        this.updateWeights = updateWeights;
    }


    protected final float exp(float x) {
        x = -x;

        if (x > expThreshold) {
            return 0;
        } else {
            return expLookup[ (int) (x / stepSize)];
        }
    }


    protected final float exp(double x) {
        x = -x;

        if (x > expThreshold) {
            return 0;
        } else {
            return expLookup[ (int) (x / stepSize)];
        }
    }


    public final boolean isValid() {
        return valid;
    }


    protected final void updateProbabilities(PointChooser pc) {
        double pX, temp, chiSquared;
        int i, j, d;


        // update the pXJ's
        for (j = 0; j < nGaussians; j++) {
            if (variable[j]) {
                for (d = 0; d < nDimensions; d++) {
                    cjjInverse[j][d] = 1 / cjj[j][d];
                }
                double covarianceDeterminant = 1;
                for (d = 0; d < nDimensions; d++) {
                    covarianceDeterminant *= cjj[j][d];
                }
                normalizationFactor[j] = 1 / (factor2pi * Math.sqrt(covarianceDeterminant));
            }
        }

        for (i = 0; i < pc.getNChosen(); i++) {
            // calculate pX (the evidence, that is, the normalization)
            // calculate pXJtemp (the likelihood)
            // pJ is the prior probability, and is calculated in iteration()
            pX = 0;
            for (j = 0; j < nGaussians; j++) {
                if (variable[j]) {
                    chiSquared = 0;
                    for (d = 0; d < nDimensions; d++) {
                        temp = x[d][pc.get(i)] - mean[j][d];
                        chiSquared += temp * temp * cjjInverse[j][d];
                    }
                    pXJtemp[j] = (float) (normalizationFactor[j] * exp( -0.5 * chiSquared));
                    pX += pXJtemp[j] * pJ[j];
                }
            }

            // calculate PJX (the posterior probability)
            for (j = 0; j < nGaussians; j++) {
                if (variable[j]) {
                    if (pX != 0) {
                        pJX[j][pc.get(i)] = (float) (pXJtemp[j] * pJ[j] / pX);
                    }
                }
            }
        }

        
        ///print spikes per cluster start
        //SPIKES WITH ZERO PROBABILITY FOR ALL CLUSTERS GET PUT IN CLUSTER O
/*				int[] nSpikes = new int[pJX.length];
        for(int spike=0; spike<pc.getNChosen(); spike++) {
            double maxProb = -1;
            int maxClust = -1;
            for(int cluster=0; cluster<pJX.length; cluster++) {
                if(pJX[cluster][pc.get(spike)] > maxProb) {
                    maxClust = cluster;
                    maxProb = pJX[cluster][pc.get(spike)];
                }
            }
            nSpikes[maxClust]++;
        }

        System.out.println();
        for(int cluster=0; cluster<nSpikes.length; cluster++) {
            System.out.println("Cluster: " + cluster + " Spike count: " + nSpikes[cluster]);
        }*/
        ///print spikes per cluster end
    }


    public int fit(double likelihoodError, int minIterations, int maxIterations) throws
    FitFailedException {
    

        double s1 = System.currentTimeMillis();
        int iter = 0;

        checkVariances();
        if(collapsedCount == nGaussians) return iter;
        iteration(reducedPC);
        double old_l = likelihood(reducedPC);
        checkVariances();
        if(collapsedCount == nGaussians) return iter;

        for (iter = 1; iter <= maxIterations; iter++) {

            iteration(reducedPC);
            checkVariances();
            if(collapsedCount == nGaussians) {
                System.err.println(identifier + "All gaussians collapsed on electrode!");
                return iter;
            }

            double l = likelihood(reducedPC);
            if (Double.isNaN(l)) {
                reset(x, fullPC.getNTotal(), fullPC.getNTotal());
                throw new FitFailedException("The likelihood is NaN");
            }
            double relativeLikelihoodChange = Math.abs( (l - old_l) / old_l);
            if (generateReport) {
                System.out.println(iter + ": " + l + " : " + relativeLikelihoodChange);
            }
            if (iter > minIterations) {
                if (relativeLikelihoodChange < likelihoodError) {
                    break;
                } else {
                    old_l = l;
                }
            }
        }

        double s2 = System.currentTimeMillis();

//		System.out.println(
//				"EM took: " +  (s2 - s1) / 1000.0 + "s for " +
//				iter + " iterations");

        if (reducedPC.getNChosen() != reducedPC.getNTotal()) {
            updateProbabilities(fullPC);
        }

        valid = true;
        return iter;
    }

    public void setIdentifier(String identifier) {
        this.identifier = identifier;
    }

    // see if a cluster collapsed (or if we're mapping, a cluster might be collapsed before we begin)
    private void checkVariances() {
        for (int j = 0; j < nGaussians; j++) {
            for (int d = 0; d < nDimensions; d++) {
                if ((variable[j] == true) && ((Math.sqrt(cjj[j][d]) / variance[d] < 1e-3)  || (Double.isNaN(cjj[j][d])))) {
                    System.err.println(identifier + "Gaussian " + (j + 1) +
                            " collapsed to 0 variance on dimension " + (d + 1) + ".");
                    variable[j] = false;
                    collapsedCount++;
                }
            }
        }
    }


    public final void reset(float[][] x, int reducedNPoints, int fullNPoints) {
        this.x = x;
        
        fullPC = new PointChooser(fullNPoints, fullNPoints);
        reducedPC = new PointChooser(reducedNPoints, fullNPoints);	
    
        
        this.nGaussians = 0;
        Arrays.fill(variance, 0);
        Arrays.fill(normalizationFactor, 0);
        Arrays.fill(pXJtemp, 0);
        Arrays.fill(pJ, 0);
        Arrays.fill(variable, true);

        for (int j = 0; j < maxGaussians; j++) {
            Arrays.fill(mean[j], 0);
            Arrays.fill(cjj[j], 0);
            Arrays.fill(cjjInverse[j], 0);
            Arrays.fill(pJX[j], 0);
        }

        factor2pi = Math.pow(2 * Math.PI, nDimensions / 2.0);
        valid = false;

        calculateVariances();
    }


    public final void reset(float[][] x, int reducedNPoints, int fullNPoints,
            ClusteringModelFile.EMModel model) {
        this.x = x;

        fullPC = new PointChooser(fullNPoints, fullNPoints);
        reducedPC = new PointChooser(reducedNPoints, fullNPoints);

        this.nGaussians = model.nGaussians;
        System.arraycopy(model.probability, 0, pJ, 0, model.probability.length);
        Arrays.fill(variable, true);

        
        for (int j = 0; j < nGaussians; j++) {
            System.arraycopy(model.covariances[j], 0, cjj[j], 0, nDimensions);
            System.arraycopy(model.means[j], 0, mean[j], 0, nDimensions);
            Arrays.fill(pJX[j], 0);
        }

        factor2pi = Math.pow(2 * Math.PI, nDimensions / 2.0);
        updateProbabilities(fullPC);

        calculateVariances();
        valid = true;
    }


    private void calculateVariances() {
        MeanVarianceCalculator mvc = new MeanVarianceCalculator();
        for (int d = 0; d < nDimensions; d++) {
            mvc.reset();
            for (int i = 0; i < fullPC.getNChosen(); i++) {
                mvc.add(x[d][fullPC.get(i)]);
            }
            variance[d] = mvc.getStandardDeviation();
        }
    }


    public abstract void addGaussian(double[] m, double[] c);


    public abstract void removeGaussian(int gIndex);


    public void removeClosestGaussian(double x, double y, int dim1, int dim2) {
        if (getClustersCount() == 0) {
            return;
        }

        double dMin = Double.POSITIVE_INFINITY;
        int jMin = -1;

        for (int j = 0; j < getGaussiansCount(); j++) {
            double[] mean = getMeans(j);
            double[] cov = getCovariances(j);
            double rho =
                Math.pow(x - mean[dim1], 2) / cov[dim1] +
                Math.pow(y - mean[dim2], 2) / cov[dim2];
            double d = Math.abs(1 - rho);
            if (d < dMin) {
                dMin = d;
                jMin = j;
            }
        }

        removeGaussian(jMin);
    }


    public abstract int getClustersCount();


    public final int getGaussiansCount() {
        return nGaussians;
    }


    public final int getPointsCount() {
        return fullPC.getNChosen();
    }


    public final int getDimensionsCount() {
        return nDimensions;
    }


    /**
     * 
     * 
     * @param i  spike, assumed to be an index valid for all spikes (That is,
     * not from a reducedPC)
     * @param pc  spike chooser
     * @return cluster that the spike has the greatest probability of coming from.
     */
    public abstract int getCluster(int i);

    /**
     * 
     * This function is not called locally.
     * 
     * @param xx Points
     * @param m  Means
     * @param s2 Variance Squared
     * @param f ??  Some sort of normalization factor
     * @return
     */
    public static final float _pXJ(double[] xx, double[] m, double[] s2, double f) {
        float chi2 = 0;
        float detCov = 1; //determinant of covariance matrix

        for (int d = 0; d < xx.length; d++) {
            chi2 += (xx[d] - m[d]) * (xx[d] - m[d]) / s2[d];
            detCov *= s2[d];
        }

        return (float) (StrictMath.exp( -0.5 * chi2) / (f * Math.sqrt(detCov)));
    }

    // Return the log  instead of the likelihood.
    // --tamachado
    //Not called from Java - mgrivich
    public static final float _lpXJ(double[] xx, double[] m, double[] s2, double f) {
        float chi2 = 0;
        float detCov = 1;

        for (int d = 0; d < xx.length; d++) {
            chi2 += (xx[d] - m[d]) * (xx[d] - m[d]) / s2[d];
            detCov *= s2[d];
        }

        return (float) ( -0.5 * chi2 - StrictMath.log(f * Math.sqrt(detCov)));
    }


    protected abstract void iteration(PointChooser pc) throws FitFailedException;


    public final float likelihood(PointChooser pc) {
        float l = 0;
        double chiSquared, temp, pXJ;
    
        for (int i = 0; i < pc.getNChosen(); i+= VisionParams.LIKELIHOOD_CALC_STRIDE) { 
            float pi = 0;
            for (int j = 0; j < nGaussians; j++) {
                // do not process collapsed gaussians
                if(!variable[j])
                    continue;
                // compute the likelihood
                chiSquared = 0;
                for (int d = 0; d < nDimensions; d++) {
                    temp = x[d][pc.get(i)] - mean[j][d];
                    chiSquared += temp * temp * cjjInverse[j][d];
                }
                pXJ = (float) (normalizationFactor[j] * exp( -0.5 * chiSquared));

                pi += pJ[j] * pXJ;
            }
            if (pi != 0) {
                l += Math.log(pi);
            }
        }
        return l;
    }

    /**
     * Get means of data for all dimensions
     * 
     * @param j
     * @return
     */
    public final double[] getMeans(int j) {
        double[] m = new double[nDimensions];
        System.arraycopy(mean[j], 0, m, 0, nDimensions);
        return m;
    }

    /**
     * Get covariances of data for all dimensions
     * 
     * @param j
     * @return
     */
    public final double[] getCovariances(int j) {
        double[] c = new double[nDimensions];
        System.arraycopy(cjj[j], 0, c, 0, nDimensions);
        return c;
    }
    
    /**
     * Get mean of particular gaussian and dimension
     * 
     * @param gaussian
     * @param dimension
     * @return
     */
    public double getMean(int gaussian, int dimension) {
        return mean[gaussian][dimension];
    }
    
    /**
     * Set mean of particular gaussian and dimension
     * 
     * @param gaussian
     * @param dimension
     * @param newMean
     */
    
    public void setMean(int gaussian, int dimension, double newMean) {
        mean[gaussian][dimension] = newMean;
    }


    public final double getGaussianProbability(int j) {
        return pJ[j];
    }


    public final ParametricEllipse[] getEllipses(int dim1, int dim2) {
        ParametricEllipse[] e = new ParametricEllipse[nGaussians];
        for (int j = 0; j < nGaussians; j++) {
            double sx = Math.sqrt(cjj[j][dim1]);
            double sy = Math.sqrt(cjj[j][dim2]);
            float angle = 0;

            e[j] = new ParametricEllipse(
                    mean[j][dim1], mean[j][dim2], sx, sy, angle);
        }
        return e;
    }


    public final void setGenerateReport(boolean generateReport) {
        this.generateReport = generateReport;
    }


    


    public final PlotPanel getInitialConditions() {
        PlotPanel p = new PlotPanel();
        p.setAxisVisible(false);

        ScatterPlot sp = new ScatterPlot();
        for (int i = 0; i < fullPC.getNChosen(); i++) {
            sp.add(x[0][fullPC.get(i)], x[1][fullPC.get(i)]);
        }

        p.addData(sp, new ScatterPlotStyle(
                SymbolType.SQUARE, 1, Color.black, false, Color.black, 1));

        for (int j = 0; j < nGaussians; j++) {
            //            double sx = Math.sqrt(cjj[j][0]);
            //            double sy = Math.sqrt(cjj[j][1]);
            double sx = 30;
            double sy = 30;
            p.addData(new ParametricEllipse(mean[j][0], mean[j][1], sx, sy, 0),
                    new FunctionStyle("", Color.red, 4));
        }

        p.setSize(350, 350);
        p.autoscale();

        return p;
    }


    public final PlotPanel getResults(int symbolSize) {
        PlotPanel plot = new PlotPanel();
        plot.setAxisVisible(false);

        ScatterPlot[] sp = new ScatterPlot[nGaussians];
        for (int i = 0; i < sp.length; i++) {
            sp[i] = new ScatterPlot();
            plot.addData(sp[i], new ScatterPlotStyle(
                    SymbolType.SQUARE, symbolSize,
                    PlotUtil.getColor(i), false, Color.black, 1));
        }

        for (int i = 0; i < fullPC.getNChosen(); i++) {
            int cluster = getCluster(i);
            if (cluster != -1) {
                sp[cluster].add(x[0][fullPC.get(i)], x[1][fullPC.get(i)]);
            }
        }

        ParametricEllipse[] ellipse = getEllipses(0, 1);
        for (int j = 0; j < nGaussians; j++) {
            plot.addData(ellipse[j], new FunctionStyle("", Color.red, 3));
        }

        plot.setSize(350, 350);
        plot.autoscale();

        //        PlotUtilities.showData("" + likelihood(), plot);
        return plot;
    }

}
