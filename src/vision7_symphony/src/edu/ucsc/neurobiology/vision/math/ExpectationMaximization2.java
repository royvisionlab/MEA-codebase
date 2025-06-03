package edu.ucsc.neurobiology.vision.math;


/**
 * This class implements the Expectation Maximization (EM) algorithm for mixtures of
 * multivariate gaussians.  There are two gaussians per cluster.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ExpectationMaximization2
    extends ExpectationMaximization {


    public ExpectationMaximization2(
        int nDimensions, int _maxGaussians, int maxNPoints) {

        super(nDimensions, _maxGaussians, maxNPoints, null);
    }


    public ExpectationMaximization2(
        int nDimensions, int _maxGaussians, int maxNPoints, float[][] _pJX) {

        super(nDimensions, _maxGaussians, maxNPoints, _pJX);
    }


    public void addGaussian(double[] m, double[] c) {
        _addGaussian(m, c);
        MathUtil.multiply(c, 4);
        _addGaussian(m, c);

        valid = false;
    }


    private void _addGaussian(double[] m, double[] c) {
        if (nGaussians == maxGaussians) {
//            System.out.println("EM is full.");
            return;
        }

        // set the means and covariances
        for (int d = 0; d < nDimensions; d++) {
            this.mean[nGaussians][d] = m[d];
            this.cjj[nGaussians][d] = c[d];
        }

        // increment the number of components
        nGaussians++;

        // redistribute the probabilities
        for (int j = 0; j < nGaussians; j++) {
            this.pJ[j] = 1.0 / nGaussians;
        }
    }


    public void removeGaussian(int gIndex) {
        boolean[] toRemove = new boolean[nGaussians];

        System.out.println("Remove Cluster: " + (gIndex % 2));

        if (gIndex % 2 == 0) {
            toRemove[gIndex + 0] = true;
            toRemove[gIndex + 1] = true;
        } else {
            toRemove[gIndex - 1] = true;
            toRemove[gIndex + 0] = true;
        }

        int index = 0;
        for (int j = 0; j < nGaussians; j++) {
            if (!toRemove[j]) {
                for (int d = 0; d < nDimensions; d++) {
                    mean[index][d] = mean[j][d];
                    cjj[index][d] = cjj[j][d];
                }
                index++;
            }
        }

        nGaussians = index;

        // redistribute the probabilities
        for (int j = 0; j < nGaussians; j++) {
            pJ[j] = 1.0 / nGaussians;
        }

        valid = false;
    }


    public int getClustersCount() {
        return nGaussians / 2;
    }


    public final int getCluster(int i) {
        float maxProbability = pJX[0][i];
        int maxJ = 0;

        for (int j = 1; j < nGaussians; j++) {
            if (pJX[j][i] > maxProbability) {
                maxProbability = pJX[j][i];
                maxJ = j;
            }
        }

        return maxJ / 2;
    }


    protected void iteration(PointChooser pc) {
        updateProbabilities(pc);

        double temp, temp2;

        // calculate new mean & variance
        for (int j = 0; j < nGaussians; j++) {
            float denom = 0;
            for (int i = 0; i < pc.getNChosen(); i++) {
                denom += pJX[j][pc.get(i)];
            }

            if (updateMeans && variable[j]) {
                // mean
                for (int d = 0; d < nDimensions; d++) {
                    temp = 0;
                    for (int i = 0; i < pc.getNChosen(); i++) {
                        temp += pJX[j][pc.get(i)] * x[d][pc.get(i)];
                    }
                    mean[j][d] = temp / denom;

                    if (j % 2 == 1) {
                        mean[j][d] = mean[j - 1][d];
                    }
                }
            }

            if (updateVariances && variable[j]) {
                // variance
                for (int d = 0; d < nDimensions; d++) {
                    temp = 0;
                    for (int i = 0; i < pc.getNChosen(); i++) {
                        temp2 = x[d][pc.get(i)] - mean[j][d];
                        temp += pJX[j][pc.get(i)] * temp2 * temp2;
                    }
                    cjj[j][d] = temp / denom;

                    if (j % 2 == 1) {
                        if (cjj[j][d] > 9 * cjj[j - 1][d]) {
                            cjj[j][d] = 9 * cjj[j - 1][d];
                        }
                    }
                }
            }

            if (updateWeights && variable[j]) {
                pJ[j] = denom / pc.getNChosen();
            }
        }
    }

}
