package edu.ucsc.neurobiology.vision.test;

import static edu.ucsc.neurobiology.vision.plot.SymbolType.NONE;

import java.io.IOException;

import edu.ucsc.neurobiology.vision.analysis.TimeCourseCalculator;
import edu.ucsc.neurobiology.vision.io.NeuronFile;
import edu.ucsc.neurobiology.vision.io.ParametersFile;
import edu.ucsc.neurobiology.vision.io.STAFile;
import edu.ucsc.neurobiology.vision.math.CovarianceMatrix;
import edu.ucsc.neurobiology.vision.math.MathUtil;
import edu.ucsc.neurobiology.vision.math.PCA;
import edu.ucsc.neurobiology.vision.math.TooManyIterationsException;
import edu.ucsc.neurobiology.vision.math.fitting.STATimeFunction1;
import edu.ucsc.neurobiology.vision.plot.DoubleHistogram2D;
import edu.ucsc.neurobiology.vision.plot.FunctionStyle;
import edu.ucsc.neurobiology.vision.plot.HistogramStyle;
import edu.ucsc.neurobiology.vision.plot.PlotPanel;
import edu.ucsc.neurobiology.vision.plot.PlotUtil;
import edu.ucsc.neurobiology.vision.plot.ScatterPlot;
import edu.ucsc.neurobiology.vision.plot.ScatterPlotStyle;
import edu.ucsc.neurobiology.vision.stimulus.STA;
import edu.ucsc.neurobiology.vision.stimulus.STAFrame;
import edu.ucsc.neurobiology.vision.util.VisionParams;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class NeuronsTest {
    NeuronFile nf;
    ParametersFile pf;
    public STAFile sf;

    ScatterPlotStyle style;
    PlotPanel p;


    public NeuronsTest(String fileName) throws IOException {
        nf = new NeuronFile(fileName + VisionParams.NEURON_FILE_EXTENSION);
//        pf = new ParametersFile(fileName + ".params");
        sf = new STAFile(fileName + VisionParams.STA_FILE_EXTENSION);
    }


    //Dumitru, is this code still used?  -- Matthew
    public double staPCA(int id) throws IOException {
        final int nEigenvectors = 5;

        System.out.println("Neuron: " + id);
        STA sta = sf.getSTA(id);
        final int w = sta.getWidth();
        final int h = sta.getHeight();
        final int staDepth = sta.size();
        double[][] tCourses = new double[w * h][3 * staDepth];

        // extract the time-courses
        ScatterPlotStyle style = new ScatterPlotStyle();
        style.setSymbolType(NONE);
        style.setConnectingPoints(true);
        for (int c = 0; c < 3; c++) {
            for (int i = 0; i < w; i++) {
                for (int j = 0; j < h; j++) {
                    for (int f = 0; f < staDepth; f++) {
                        tCourses[i * h + j][staDepth * c + f] =
                            sta.getFrame(f).getPixel(i, j, c);
                    }
                }
            }
        }

        // do the PCA
        CovarianceMatrix c = new CovarianceMatrix(3 * staDepth);
        for (int i = 0; i < tCourses.length; i++) {
            c.addData(tCourses[i]);
        }
        float[] cov = c.getCovariance();
        PCA pca = new PCA(cov);
        try {
            pca.doPCA();
        } catch (TooManyIterationsException e) {
            System.out.println("---------------- Impossible to Analyze ---------------");
            return Double.NaN;
        }

        // show the eigenvectors
//        pca.printPercentageEigenvalues(nEigenvectors);
        pca.showEigenVectors(nEigenvectors);

//        double[] pScoreArray = new double[nEigenvectors];
//        double[] nScoreArray = new double[nEigenvectors];
//        double[] maxScoreArray = new double[nEigenvectors];
//        double[] complexityIndex = new double[nEigenvectors];

        for (int eigen = 0; eigen < nEigenvectors; eigen++) {
            STAFrame frame = new STAFrame(w, h, 60.0, 60.0);
            for (int i = 0; i < w; i++) {
                for (int j = 0; j < h; j++) {
                    double x = pca.project(tCourses[i * h + j], eigen);
                    frame.setPixel(i, j, new float[] { (float) x, (float) x, (float) x});
                }
            }

            // obtain the scores from this eigenvector
            DoubleHistogram2D h2d = new DoubleHistogram2D("", 0, w, 0, h, 1, 1);
//            DoubleHistogram scoreHist = new DoubleHistogram("", -.1, .1, .00025);
//            double[][] scores = new double[1][w * h];
            double maxScore = Double.NEGATIVE_INFINITY;
            double minScore = Double.POSITIVE_INFINITY;
            for (int i = 0, k = 0; i < w; i++) {
                for (int j = 0; j < h; j++) {
                    double score = pca.project(tCourses[i * h + j], eigen);
//                    scoreHist.fill(score, 1);
//                    scores[0][k++] = score;
                    h2d.fill(i, h - j - 1, score);

                    if (score > maxScore) {
                        maxScore = score;
                    }
                    if (score < minScore) {
                        minScore = score;
                    }
                }
            }

            if (maxScore < 0 || minScore > 0) {
                System.out.println("---> problem");
            } else if (maxScore > -minScore) {
                System.out.println("---> ok");
            } else {
                System.out.println("---> has to be inverted");
            }

            PlotUtil.showData("" + (eigen + 1), h2d, new HistogramStyle());
        }
        /*
                double[] sigma = new double[] {
                                 Math.sqrt(em.getCovariances(0)[0]),
                                 Math.sqrt(em.getCovariances(1)[0])};
                final int surroundCluster = VisionUtilities.maxIndex(sigma);
                final int centerCluster = VisionUtilities.minIndex(sigma);
                double centerMean = em.getMeans(centerCluster)[0];
                // calculate the scores above n*Sigma
                double score = 0, pScore = 0, nScore = 0;
                for (int i = 0; i < w * h; i++) {
             double d = Math.abs( (scores[0][i] - centerMean) / sigma[centerCluster]);
                    if (d > 4) {
                        // if (em.getCluster1(i, 0.99) == surroundCluster) {
                        score += Math.abs(scores[0][i]);
                        if (scores[0][i] > 0) {
                            pScore += scores[0][i];
                        } else {
                            nScore += -scores[0][i];
                        }
                    }
                }
                pScoreArray[eigen] = pScore;
                nScoreArray[eigen] = nScore;
                maxScoreArray[eigen] = pScore + nScore;
                complexityIndex[eigen] = Math.abs(pScore - nScore) / (pScore + nScore);
//            ArrayList frames = new ArrayList();
//            frames.add(frame);
//            STA sta1 = new STA(frames, 1);
//            PlotPanel staPanel = new PlotPanel();
//            staPanel.setAxisVisible(false);
//            staPanel.setRange(0, 60 * sta1.getWidth(), 0, 60 * sta1.getHeight());
//            staPanel.addDataPlotter(new ColorPlotPlotter());
//            SimpleColorPlot scp = new SimpleColorPlot(SimpleColorPlot.NORMALIZE);
//            scp.setFrame(sta1, sta1.maxFrame);
//            BasicColorPlotStyle style = new BasicColorPlotStyle();
//            staPanel.addData(scp, style);
//            PlotUtilities.showData("pc " + (eigen + 1), staPanel);
            }
            //        VisionUtilities.printArray(maxScore);
            boolean good = false;
            for (int i = 0; i < nEigenvectors; i++) {
                if (maxScoreArray[i] != 0) {
                    good = true;
                }
            }
            if (!good) {
                System.out.println("Bad STA");
                return Double.NaN;
            }
            final int maxEigen = VisionUtilities.maxIndex(maxScoreArray);
            System.out.println("Eiegenvector " + (maxEigen + 1) +
                               ", c = " + complexityIndex[maxEigen]);
            if (complexityIndex[maxEigen] >= 0.99) {
                int sign = (pScoreArray[maxEigen] > nScoreArray[maxEigen]) ? +1 : -1;
                System.out.println("sign = " + sign);
                double[] eigenvector = VisionUtilities.multiply(
                    pca.getEigenVector(maxEigen), sign);
                p.addData(new ScatterPlot("" + id, eigenvector), style);
                p.autoscale();
            }
            return complexityIndex[maxEigen];
         */
        return 0;
    }


    //Dumitru, is this code still used?   -- Matthew
    public STATimeFunction1 fitSTA(int id, STATimeFunction1 fun) throws
        Exception {

        STA sta = sf.getSTA(id);
        int[] p = sta.getMainFrameParams();
        final int staDepth = sta.size();
        double[] tCourse = new double[staDepth];
        double[] error = new double[staDepth];

        // extract the time-courses
        ScatterPlot sp = new ScatterPlot("");
        double[][] x = new double[1][staDepth];

        tCourse = sta.getTimeFilters(5.0)[1];
        error = sta.getTimeFiltersError()[1];
        MathUtil.multiply(error, 4);
        for (int f = 0; f < staDepth; f++) {
//            tCourse[f] = sta.getFrame(f).getPixel(p[1], p[2], 1);
//            error[f] = sta.getFrame(f).getPixelError(p[1], p[2], 1);
            sp.add( - (staDepth - f - 1) * 8.34, tCourse[f], error[f]);
            x[0][f] = - (staDepth - f - 1) * 8.34;
        }

        int max = MathUtil.maxAbsIndex(tCourse);

        PlotPanel pp = new PlotPanel();
        pp.addData(sp, new ScatterPlotStyle());
        pp.addData(fun, new FunctionStyle("f"));
        pp.autoscale();
        pp.padY();
        PlotUtil.showData("", pp, 600, 300);

        return fun;
    }


    public STATimeFunction1 showTC(double[] timeCourse, double[] timeCourseError,
                                   double refreshTime, STATimeFunction1 fun) throws
        Exception {

        // extract the time-courses
        ScatterPlot sp = new ScatterPlot("");
        double[][] x = new double[1][timeCourse.length];

        //code for STV
//         tCourse = sta.getComparisonTimeFilters(5.0, stv)[1];
//
//         MathUtil.sub(tCourse, MathUtil.mean(tCourse));
//        MathUtil.multiply(error, .4);

        for (int f = 0; f < timeCourse.length; f++) {
//            tCourse[f] = sta.getFrame(f).getPixel(p[1], p[2], 1);
//            error[f] = sta.getFrame(f).getPixelError(p[1], p[2], 1);
            sp.add( - (timeCourse.length - f - 1) * refreshTime, timeCourse[f],
                   timeCourseError[f]);
            x[0][f] = - (timeCourse.length - f - 1) * refreshTime;
        }

        PlotPanel pp = new PlotPanel();
        pp.addData(sp, new ScatterPlotStyle());
        pp.addData(fun, new FunctionStyle("f"));
        pp.autoscale();
        pp.padY();
        PlotUtil.showData("", pp, 600, 300);

        return fun;
    }


    public static void main(String[] args) throws Exception {
        final NeuronsTest ntest = new NeuronsTest(
            "g:\\data\\2005-05-02-0\\data002\\data002");
//        final NeuronsTest ntest = new NeuronsTest("g:\\data\\2005-08-25-0\\data002\\data002");


        int id = 3;
        int[] n1n2n3 = new int[3];
        STA sta = ntest.sf.getSTA(id);
        double[] timeCourse = sta.getTimeFilters(5.0)[1];
        double[] timeCourseError = sta.getTimeFiltersError()[1];

        STATimeFunction1 f = TimeCourseCalculator.fitSTATimeFilter(timeCourse,
            timeCourseError, false, 0,
            30, 0, 16, 8, 16,
            n1n2n3, sta.getRefreshTime(), 0.0);

//        STATimeFunction2 fu = ntest.fitSTATimeFilter2(timeCourse, timeCourseError, false,
//            3, sta.getRefreshTime());

//        System.err.println(fu);
        System.err.println(n1n2n3[0] + ", " + n1n2n3[1] + "," + n1n2n3[2]);

        ntest.showTC(timeCourse, timeCourseError, sta.getRefreshTime(), f);

    }

}
