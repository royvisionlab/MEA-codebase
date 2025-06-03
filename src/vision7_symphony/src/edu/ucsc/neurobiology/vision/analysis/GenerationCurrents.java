package edu.ucsc.neurobiology.vision.analysis;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.calculations.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.math.fitting.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.stimulus.*;
import edu.ucsc.neurobiology.vision.stimulus.AdvancedMovie.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * This calculation generates a series of plots of the average spike rate versus
 * the generation current. The generation current is calculated as the dot product
 * of the STA with the image projected on the screen before the occurrence of a spike.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class GenerationCurrents
    extends AbstractCalculation {

//    private static double minCurrent = -.02, maxCurrent = +.02,
//        currentBinning = 0.00005;

    private static double minCurrent = -10, maxCurrent = +10;
    private static double currentBinning = 0.05;
    private static int nBins = (int) Math.ceil( (maxCurrent - minCurrent) /
                                               currentBinning);


    private NeuronFile neuronFile;
    private STAFile staFile;
    private Vision app;
    private BufferedMovie movie;
    private String fileNameNoExt;
    private Date start;
    private String filePath, neuronFileName, staFileName, movieFileName, globalsFileName;
    private int neuronID;


    public void startCalculation() throws IOException {
        app = Vision.getInstance();
        app.sendMessage("Sr(g)...");

        String datasetName = new File(filePath).getName();

        neuronFileName = filePath + File.separator + datasetName +
                         VisionParams.NEURON_FILE_EXTENSION;
        staFileName = filePath + File.separator + datasetName +
                      VisionParams.STA_FILE_EXTENSION;
        movieFileName = filePath + File.separator + datasetName +
                        VisionParams.MOVIE_FILE_EXTENSION;   
        globalsFileName = filePath + File.separator + datasetName +
        ".globals";

        // Get the spike, neurons and sta's.
        neuronFile = new NeuronFile(neuronFileName);
        movie = new BufferedMovie(new WhiteNoiseMovie(movieFileName, globalsFileName));
        movie.setFrameEncoding(FrameEncoding.FLOAT_ARRAY_ENCODING);
        movie.reset();
        staFile = new STAFile(staFileName);

        start = new Date();

        // create the histogram file
//        String nlfFileName = filePath + File.separator + datasetName + ".nlf";
//        HistogramFile f = null;
//        f = new HistogramFile(nlfFileName, minCurrent, maxCurrent, nBins);

        // calculate the non-linear functions
        if (neuronID < 0) {
            int[] idList = neuronFile.getIDList();
            for (int i = 0; i < idList.length; i++) {
                app.sendMessage("Calculating Sr(I), id " + idList[i] +
                                "(" + i + "/" + idList.length + ")");
                DoubleErrorHistogram[] srg = calculate(new int[] {idList[i]}
                    , staFile, neuronFile, movie);
//                f.addHistogram(idList[i], srg[0]);
            }
        } else {
            app.sendMessage("Calculating Sr(I), id " + neuronID);
            DoubleErrorHistogram[] srg = calculate(new int[] {neuronID}
                , staFile, neuronFile, movie);
//            f.addHistogram(neuronID, srg[0]);
        }

        // close the histogram file
//        f.close();

        // end the calculation
        Date end = new Date();
        double deltaTime = (end.getTime() - start.getTime()) / 1000.;
        app.sendMessage("Sr(g) Finished in " + deltaTime + "s");
        Vision.getInstance().getCalculationManager().calculationDone();
    }


    private static DoubleErrorHistogram[] calculate(
        int[] idList, STAFile staFile, NeuronFile neuronFile, BufferedMovie movie) throws
        IOException {

        movie.reset();
        double refreshTime = movie.getRefreshTime();
        double samples2ms = 1.0 / 20.0;

        Vision app = Vision.getInstance();
        double firstImageChangeTime = neuronFile.getTTLTimes()[0] * samples2ms +
                                      refreshTime;

        float[][][] staFrames = new float[idList.length][][];
        DoubleHistogram[] spikesCount = new DoubleHistogram[idList.length];
        DoubleHistogram[] spikeCountVsCurrentHist = new DoubleHistogram[idList.length];
        DoubleHistogram[] currentHist = new DoubleHistogram[idList.length];

        ScatterPlot sp = new ScatterPlot();

        for (int i = 0; i < idList.length; i++) {
            spikeCountVsCurrentHist[i] = new DoubleHistogram(
                "Spike Rate Vs. Current " + idList[i], minCurrent, maxCurrent,
                currentBinning);
            currentHist[i] = new DoubleHistogram(
                "Current " + idList[i], minCurrent, maxCurrent, currentBinning);

            // load the STA
            STA sta = staFile.getSTA(idList[i]);
            staFrames[i] = new float[sta.size()][];
            for (int j = 0; j < sta.size(); j++) {
                staFrames[i][j] = sta.getFrame(j).getBuffer();
            }

            // calculate the spike rate
            int[] spikeTimes = neuronFile.getSpikeTimes(idList[i]);
            spikesCount[i] = new DoubleHistogram(
                "", firstImageChangeTime, neuronFile.getNumberOfSamples() * samples2ms,
                refreshTime);
            for (int j = 0; j < spikeTimes.length; j++) {
                spikesCount[i].fill(spikeTimes[j], 1);
            }
        }

        // go through the data frame by frame and calculate the Sr(g) scatter plots
        int oldProgress = -1;
        final int nFrames = spikesCount[0].getBinCount();

        // frame zero os happening at firstImageChangeTime
        for (int frame = staFrames.length; frame < nFrames; frame++) {
            for (int neuron = 0; neuron < idList.length; neuron++) {
                double current = calculateCurrent(staFrames[neuron], movie, frame);
                if (Double.isNaN(current)) {
                    continue;
                }

                sp.add(current,
                       spikesCount[neuron].getBin(frame) +
                       0.25 * 2 * (Math.random() - 0.5));

                spikeCountVsCurrentHist[neuron].fill(current,
                    spikesCount[neuron].getBin(frame));
                currentHist[neuron].fill(current, 1);

                // update the progress
                int progress = (frame * 100) / nFrames;
                if (progress != oldProgress) {
                    app.setProgress(progress);
                    oldProgress = progress;
                }
            }
        }

        // normalize the histograms (get the average of each bin)
        final double factor = 1000.0 / refreshTime;
        DoubleErrorHistogram[] srt = new DoubleErrorHistogram[idList.length];
        for (int neuron = 0; neuron < idList.length; neuron++) {
            DoubleHistogram valueHist = new DoubleHistogram(
                "Sr(t) value" + idList[neuron], minCurrent, maxCurrent, currentBinning);
            DoubleHistogram errorHist = new DoubleHistogram(
                "Sr(t) error" + idList[neuron], minCurrent, maxCurrent, currentBinning);
            for (int i = 0; i < nBins; i++) {
                double x = spikeCountVsCurrentHist[neuron].getBin(i);
                double y = currentHist[neuron].getBin(i);
                if (y != 0) {
                    // the average of the spike count for a particular current and
                    // the error of the ration (error propagation, Poisson errors)
                    // the value and the error are converted to spikes/second
                    valueHist.setBin(i, factor * x / y);
                    errorHist.setBin(i, factor * Math.sqrt(x * (x + y) / (y * y * y)));
                } else {
                    valueHist.setBin(i, 0);
                    errorHist.setBin(i, 0);
                }
            }

            srt[neuron] = new DoubleErrorHistogram(valueHist, errorHist);

            PlotUtil.showData("spikeCountVsCurrentHist",
                              spikeCountVsCurrentHist[neuron], new HistogramStyle());
            PlotUtil.showData("currentHist", currentHist[neuron], new HistogramStyle());
            PlotUtil.showData("nonlinear", srt[neuron], new HistogramStyle());
        }

        PlotUtil.showData(sp, new ScatterPlotStyle());

        return srt;
    }


    public static double calculateCurrent(float[][] sta, BufferedMovie movie, int frameIndex) throws IOException {
        final int firstFrame = frameIndex - sta.length + 1;
        if (firstFrame < 0) return Double.NaN;

        double g = 0;
        for (int f = 0; f < sta.length; f++) {
            float[] imageFrame = (float[]) movie.getEncodedFrame(firstFrame + f);
            for (int i = 0; i < sta[f].length; i++)
                g += sta[f][i] * (imageFrame[i] /* - 0.5*/);
        }
        return g;
    }


    /*
        public static double calculateCurrent(
            float[][] sta, Movie movie, int frameIndex) {
            final int firstFrame = frameIndex - sta.length + 1;
            if (firstFrame < 0) {
                return Double.NaN;
            }
            double g = 0;
            for (int f = 0; f < sta.length; f++) {
                float[] imageFrame = movie.getFrame(firstFrame + f).getBuffer();
                for (int i = 0; i < imageFrame.length; i++) {
                    g += sta[f][i] * (imageFrame[i] - 0.5);
                }
            }
            return g;
        }
     */

    /*
        public static double calculateCurrent(
            float[][] sta, WhiteNoiseMovie movie, int frameIndex) {
            final int firstFrame = frameIndex - sta.length + 1;
            if (firstFrame < 0) {
                return Double.NaN;
            }
            float g = 0;
            float[] imageFrame = null;
            for (int f = 0; f < sta.length; f++) {
         imageFrame = (float[]) movie.getEncodedFrame(firstFrame + f, imageFrame);
                for (int i = 0; i < imageFrame.length; i++) {
                    g += sta[f][i] * imageFrame[i];
                }
            }
            return g / (movie.getWidth() * movie.getHeight());
        }
     */

    public void setParameters(HashMap<String, String> p) {
//        minCurrent = ((DoubleParameter)table.getParameter(MIN_CURRENT)).getValue();
//        maxCurrent = ((DoubleParameter)table.getParameter(MAX_CURRENT)).getValue();
//        currentBinning = ((DoubleParameter)table.getParameter(CURRENT_BINNING)).getValue();
        filePath = (String) p.get("File_Path");
        neuronID = Integer.parseInt( (String) p.get("Neuron_ID"));
    }


    public static SigmoidFunction fitNLF(DoubleErrorHistogram h) {
        int nBins = h.getBinCount();
        double min = h.getMin();
        double max = h.getMax();

        DoubleList x = new DoubleList();
        DoubleList y = new DoubleList();
        DoubleList s = new DoubleList();

        double start = min, end = max;
        int n1 = (int) ( (start - h.getMin()) / h.getBinInterval());
        int n2 = (int) ( (end - h.getMin()) / h.getBinInterval());
        for (int i = 0; i < n2 - n1; i++) {
            double error = h.getBinError(n1 + i);
            if (error != 0) {
                x.add(start + (i + 0.5) * (max - min) / nBins);
                y.add(h.getBin(n1 + i));
                s.add(error);
            }
        }

        SigmoidFunction f = new SigmoidFunction(100, 3, 0.2);

        Fitter fitter = new Fitter();
        try {
            fitter.fit(f, new double[][] {x.toArray()}
                       , y.toArray(), s.toArray());
        } catch (FitFailedException e) {
            return null;
        }

        return f;
    }


}
