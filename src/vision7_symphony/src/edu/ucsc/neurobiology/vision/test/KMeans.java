package edu.ucsc.neurobiology.vision.test;

import java.util.*;

import cern.jet.random.*;
import cern.jet.random.engine.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.plot.*;


/**
 * This class implements the k-means clustering algorithm
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class KMeans {
    private int nClusters;
    private final int nDimensions;
    private int nPoints;
    private double[] min, max, mean;
    private double[][] x; // x[i][d]
    private double[][] means; // mean[j][d]
    private Random random = new Random();
    private int[] clusterID;


    public KMeans(double[][] x, int nPoints) {
        this.nDimensions = x[0].length;
        this.nPoints = nPoints;
        clusterID = new int[nPoints];

        this.x = new double[nPoints][];
        for (int i = 0; i < nPoints; i++) {
            this.x[i] = x[i];
        }

        min = new double[nDimensions];
        max = new double[nDimensions];
        mean = new double[nDimensions];
        for (int d = 0; d < nDimensions; d++) {
            min[d] = Double.POSITIVE_INFINITY;
            max[d] = Double.NEGATIVE_INFINITY;
            mean[d] = 0;
            for (int i = 0; i < nPoints; i++) {
                if (this.x[i][d] < min[d]) {
                    min[d] = this.x[i][d];
                }
                if (this.x[i][d] > max[d]) {
                    max[d] = this.x[i][d];
                }
                mean[d] += this.x[i][d];
            }
            mean[d] /= nPoints;
        }
//        mean = new double[maxGaussians][nDimensions];
    }


    private double distance(int i, int j) {
        double dist = 0;
        for (int d = 0; d < nDimensions; d++) {
            dist += (x[i][d] - means[j][d]) * (x[i][d] - means[j][d]);
        }
        return Math.sqrt(dist);
    }


    public void cluster1(int nClusters) {
        for (int i = 0; i < 10; i++) {
            cluster(nClusters);
        }
    }


    public void cluster(int nClusters) {
        this.nClusters = nClusters;

        // initialize the cluster centers
        means = new double[nClusters][nDimensions];
        for (int j = 0; j < nClusters; j++) {
            for (int d = 0; d < nDimensions; d++) {
                //means[j][d] = min[d] + (max[d] - min[d]) * random.nextDouble();
                means[j][d] = mean[d] +
                              0.05 * (max[d] - min[d]) * (random.nextDouble() - 1);
            }
        }

        int iteration = 0;
        double oldD = -1;
        while (true) {
            // Assign each point to its nearest cluster as represented by the means
            for (int i = 0; i < nPoints; i++) {
                double minDist = Double.POSITIVE_INFINITY;
                clusterID[i] = -1;
                for (int j = 0; j < nClusters; j++) {
                    double d = distance(i, j);
                    if (d < minDist) {
                        minDist = d;
                        clusterID[i] = j;
                    }
                }
            }

            // Update the cluster centers as sample means of the points
            for (int j = 0; j < nClusters; j++) {
                for (int d = 0; d < nDimensions; d++) {
                    means[j][d] = 0;
                    int n = 0;
                    for (int i = 0; i < nPoints; i++) {
                        if (clusterID[i] == j) {
                            means[j][d] += x[i][d];
                            n++;
                        }
                    }
                    means[j][d] /= n;
                }
            }

            // calculate the distance function
            double d = 0;
            for (int j = 0; j < nClusters; j++) {
                for (int i = 0; i < nPoints; i++) {
                    if (clusterID[i] == j) {
                        d += distance(i, j);
                    }
                }
            }
            System.out.println(iteration + " : " + d);

            // chack the termination condition
            if (oldD != -1 && d >= oldD) {
                break;
            }

            oldD = d;
            iteration++;
        }

        showState(oldD);
    }


    public void showState(double iteration) {
        PlotPanel p = new PlotPanel();

        for (int j = 0; j < nClusters; j++) {
            ScatterPlot sp = new ScatterPlot("");
            PlotPanel p1 = new PlotPanel();
            ScatterPlotStyle style1 = new ScatterPlotStyle();
            style1.setConnectingPoints(true);
            style1.setSymbolType(SymbolType.NONE);
            int n = 0;
            for (int i = 0; i < nPoints; i++) {
                if (clusterID[i] == j) {
                    sp.add(x[i][0], x[i][1]);
                    p1.addData(new ScatterPlot(x[i]), style1);
                    n++;
                }
            }
            ScatterPlotStyle style = new ScatterPlotStyle();
            style.setSymbolColor(PlotUtil.getColor(j));
            p.addData(sp, style);

            PlotUtil.showData("" + n, p1);
            p1.autoscale();
        }

        for (int j = 0; j < nClusters; j++) {
            double a = 1;
            p.addData(new ParametricEllipse(means[j][0], means[j][1], a, a, 0),
                      new FunctionStyle("Ellipse"));
        }

        PlotUtil.showData("" + iteration, p);
        p.autoscale();
    }


    public static void main(String[] args) {
        double[] a = {0.3, 0.3, 0.4};
        final double fMax = MathUtil.max(a);
        double[][] m = { { -5, -5}
                       , {0, 0}
                       , {5, 5}
        };
        double b = 2;
        double[][] s = { {b, b}
                       , {b, b}
                       , {b, b}
        };
        Normal[] xNormal = new Normal[a.length];
        Normal[] yNormal = new Normal[a.length];
        for (int i = 0; i < a.length; i++) {
            xNormal[i] = new Normal(m[i][0], Math.sqrt(s[i][0]), new DRand(i));
            yNormal[i] = new Normal(m[i][1], Math.sqrt(s[i][1]), new DRand(100 - i));
        }

        int nPoints = 500;
        double[][] data = new double[nPoints][2];
        for (int i = 0; i < nPoints; i++) {
            double x = Math.random();
            double y = Math.random();
            int n = (int) (x * a.length);
            if (y * fMax < a[n]) {
                data[i][0] = xNormal[n].nextDouble() + 2 * (Math.random() - 1);
                data[i][1] = yNormal[n].nextDouble() + 2 * (Math.random() - 1);
            }
        }

        KMeans km = new KMeans(data, data.length);
        km.cluster1(3);
    }
}
