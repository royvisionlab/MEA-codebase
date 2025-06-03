package edu.ucsc.neurobiology.vision.analysis;

import java.awt.Color;
import java.io.File;
import java.io.IOException;

import edu.ucsc.neurobiology.vision.io.NeuronFile;
import edu.ucsc.neurobiology.vision.io.OMSMovieFile;
import edu.ucsc.neurobiology.vision.io.ParametersFile;
import edu.ucsc.neurobiology.vision.io.RawMovieFile;
import edu.ucsc.neurobiology.vision.io.chunk.ChunkFile;
import edu.ucsc.neurobiology.vision.io.chunk.GlobalsFile;
import edu.ucsc.neurobiology.vision.math.MathUtil;
import edu.ucsc.neurobiology.vision.math.MeanVarianceCalculator;
import edu.ucsc.neurobiology.vision.math.Num;
import edu.ucsc.neurobiology.vision.plot.DoubleErrorHistogram;
import edu.ucsc.neurobiology.vision.plot.DoubleHistogram;
import edu.ucsc.neurobiology.vision.plot.HistogramStyle;
import edu.ucsc.neurobiology.vision.plot.PlotPanel;
import edu.ucsc.neurobiology.vision.plot.ScatterPlot;
import edu.ucsc.neurobiology.vision.plot.ScatterPlotStyle;
import edu.ucsc.neurobiology.vision.plot.SymbolType;
import edu.ucsc.neurobiology.vision.stimulus.OMSFrameGenerator;
import edu.ucsc.neurobiology.vision.util.DoubleList;


/**
 *
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class OMSClassification {
    ParametersFile paramsFile;
    double[][][] omsCenters;
    int nRuns;
    double refreshTime;
    double framesPerRun;
    int width;
    int height;
    int seed;
    int framesNeeded;
    double radius;
    double spacing;
    double maxDistance;
    double stepSize;
    int barWidth;
    int runLength;
    double refMicronsPerStixel;
    double refTotalHeight;
    double omsMicronsPerStixel;
    double xOffset;
    double yOffset;
    int[] droppedFrames = null;
    int bins = 25;
//	boolean showBackground, randomizeMotion;

    private NeuronFile omsFile;


    public OMSClassification(ParametersFile paramsFile) {
        this.paramsFile = paramsFile;
    }


    public OMSClassification(String omsFileName, String runDirectory, String refFileName) throws
    IOException {

        if (omsFileName != null) {
            String filePath =
                runDirectory + File.separator + omsFileName +
                File.separator + omsFileName + ".neurons";
            System.out.println("OMS File: " + filePath);

            omsFile = new NeuronFile(filePath);
            OMSMovieFile omsMovieFile = new OMSMovieFile(runDirectory +
                    File.separator + omsFileName +
                    File.separator + omsFileName + ".rawMovie");
            framesPerRun = omsMovieFile.getRunLength();
            width = omsMovieFile.getWidth();
            height = omsMovieFile.getHeight();
            seed = omsMovieFile.getSeed();
            radius = omsMovieFile.getObjectRadius();
            spacing = omsMovieFile.getObjectSpacing();
            stepSize = omsMovieFile.getStepSize();
            barWidth = omsMovieFile.getBarWidth();
            runLength = omsMovieFile.getRunLength();
//			showBackground = omsMovieFile.getShowBackground();
//			randomizeMotion = omsMovieFile.getRandomizeMotion();

            GlobalsFile omsGlobalsFile = new GlobalsFile(runDirectory +
                    File.separator + omsFileName +
                    File.separator + omsFileName + ".globals", ChunkFile.READ);
            GlobalsFile.RunTimeMovieParams omsRunTime = omsGlobalsFile.getRunTimeMovieParams();
            refreshTime = omsRunTime.refreshPeriod;
            framesNeeded = (int) Math.ceil( (double) omsFile.getNumberOfSamples() /
                    20.0 / refreshTime);
            nRuns = (int) Math.ceil( (double) framesNeeded / framesPerRun);
            omsMicronsPerStixel = (omsRunTime.micronsPerStixelX  + omsRunTime.micronsPerStixelY )/2;

            droppedFrames = omsRunTime.droppedFrames;
//			System.out.println("omsMicronsPerStixel: " + omsMicronsPerStixel);





            maxDistance = spacing * omsMicronsPerStixel / 1.8;

            RawMovieFile refFile = new RawMovieFile(runDirectory + File.separator +
                    refFileName +
                    File.separator + refFileName + ".rawMovie");
            refTotalHeight = refFile.getHeight();
            GlobalsFile refGlobalsFile = new GlobalsFile(runDirectory +
                    File.separator + refFileName +
                    File.separator + refFileName + ".globals", ChunkFile.READ);
            GlobalsFile.RunTimeMovieParams refRunTime = refGlobalsFile.getRunTimeMovieParams();
            refMicronsPerStixel = (refRunTime.micronsPerStixelX  + refRunTime.micronsPerStixelY )/2;

            xOffset = omsRunTime.xOffset - refRunTime.xOffset;
            yOffset = omsRunTime.yOffset - refRunTime.yOffset;



//			System.out.println("refMicronsPerStixel: " + refMicronsPerStixel);
//			System.out.println(refRunTime.micronsPerStixelX);
//			System.out.println(refRunTime.micronsPerStixelY);
//			System.out.println(refRunTime.pixelsPerStixelX);
//			System.out.println(refRunTime.pixelsPerStixelY);
//			System.out.println(xOffset);
//			System.out.println(yOffset);
        }

        float[] contrast = {.48f, .48f, .48f};

        OMSFrameGenerator omsFrameGenerator = new OMSFrameGenerator(width, height,
                contrast,
                seed,
                framesNeeded, radius,
                spacing, stepSize, barWidth, runLength, true, true);

        this.omsCenters = omsFrameGenerator.getRunCenters(nRuns);
//		for(int i=0; i<omsCenters.length;i++) {
//		for(int j=0; j<omsCenters[i].length;j++) {
//		System.out.println(omsCenters[i][j][0] + "   " + omsCenters[i][j][1]);
//		}
//		System.out.println("end of run");
//		}

    }


    synchronized public double[][] calculateOMSParams(int id, double x0, double y0) throws
    IOException {

        double[][] omsParams = new double[5][];
        if (omsFile != null) {
            int[] times = omsFile.getSpikeTimes(id);
            double[] omsDistance = new double[nRuns];
            double[] omsSpikes = new double[nRuns];
//			double[] omsError = new double[nRuns];

            double samplesPerRun = 20 * framesPerRun * refreshTime;

            double[] spikesPerRun = new double[nRuns];

            for (int i = 0; i < times.length; i++) {
                if ( ( (double) times[i] / samplesPerRun % 1) > .1) {
                    spikesPerRun[ (int) Math.floor( (double) times[i] / samplesPerRun)]++;

                }
            }
//			DoubleHistogram spikesPerRunHistogram = new DoubleHistogram(
//			"", spikesPerRun, 0.0, nRuns);
//			PlotUtilities.showData(spikesPerRunHistogram, new HistogramStyle());
            double neuronX = x0 * refMicronsPerStixel;
            double neuronY = (refTotalHeight - y0) * refMicronsPerStixel;

            for (int i = 0; i < nRuns; i++) {
                if (x0 != Double.NaN) {
                    double minDistance = Double.MAX_VALUE;

                    for (int j = 0; j < omsCenters[i].length; j++) {

                        double circleX = omsCenters[i][j][0] * omsMicronsPerStixel +
                        xOffset;
                        double circleY = omsCenters[i][j][1] * omsMicronsPerStixel +
                        yOffset;

                        minDistance = Math.min( (neuronX - circleX) *
                                (neuronX - circleX)
                                +
                                (neuronY - circleY) *
                                (neuronY - circleY), minDistance);
                    }

                    omsDistance[i] = Math.sqrt(minDistance);
                    omsSpikes[i] = spikesPerRun[i];
                } else {
                    omsDistance[i] = 0;
                    omsSpikes[i] = 0; //spikesPerRun[i];
                }
            }

//			PlotUtilities.showData("", this.makeOMSScatter(omsDistance, omsSpikes));

            double distancePerBin = maxDistance / bins;
            int[] count = new int[bins];
            double[] omsHistData = new double[bins];
            double[] omsHistErrorData = new double[bins];
            MeanVarianceCalculator[] omsMeanVar = new MeanVarianceCalculator[bins];
            for (int i = 0; i < bins; i++) {
                omsMeanVar[i] = new MeanVarianceCalculator();
            }
            for (int i = 0; i < omsDistance.length; i++) {
                omsMeanVar[Math.min( (int) (omsDistance[i] / distancePerBin),
                        bins - 1)].add(omsSpikes[i]);
                count[Math.min( (int) (omsDistance[i] / distancePerBin), bins - 1)]++;
            }
            for (int i = 0; i < count.length; i++) {

                omsHistData[i] = omsMeanVar[i].getMean();
                omsHistErrorData[i] = omsMeanVar[i].getMeanVariance();
                
                if(Double.isNaN(omsHistData[i])) {
                    omsHistData[i] = 0;
                    omsHistErrorData[i] = 0;
                }
            }

//			PlotUtilities.showData("", makeOMSHistogram(omsHistData, omsHistErrorData, maxDistance));


//			for(int i=0; i<omsCenters.length;i++) {
//			for(int j=0; j<omsCenters[i].length;j++) {
//			System.out.println(omsCenters[i][j][0] + "   " + omsCenters[i][j][1]);
//			}
//			System.out.println("end of run");
//			}

            omsParams[0] = omsDistance;
            omsParams[1] = omsSpikes;
            omsParams[2] = omsHistData;
            omsParams[3] = omsHistErrorData;
            omsParams[4] = new double[] {maxDistance};

        }

        return omsParams;
    }


    synchronized public PlotPanel makeOMSScatter(int id) {
        return makeOMSScatter(paramsFile.getArrayCell(
                id, "omsDistances"), paramsFile.getArrayCell(
                        id, "omsSpikes"));
    }


    synchronized public PlotPanel makeOMSHistogram(int id) {
        return makeOMSHistogram(paramsFile.getArrayCell(
                id, "omsHistData"), paramsFile.getArrayCell(
                        id, "omsErrorHistData"), paramsFile.getDoubleCell(id, "omsMaxDistance"));
    }


    synchronized public PlotPanel makeOMSScatter(double[] omsDistance, double[] omsSpikes) {
        double[] omsError = new double[omsDistance.length];
        ScatterPlot scatter = new ScatterPlot(omsDistance, omsSpikes, omsError);
        PlotPanel panel = new PlotPanel("omsScatter");
        panel.addData(scatter, new ScatterPlotStyle(
                SymbolType.FILLED_SQUARE, 3, Color.black, false, Color.black, 1));
        panel.setLabels("Distance (microns)", "Spikes/Run");
        panel.autoscale();
        panel.setRange(panel.getRange()[0], panel.getRange()[1], 0, panel.getRange()[3]);
        return panel;
    }


    synchronized public PlotPanel makeOMSHistogram(
            double[] omsHistData, double[] omsHistErrorData, double maxDistance) {
        PlotPanel panel = new PlotPanel("omsHist");

        DoubleHistogram omsHist = new DoubleHistogram("", omsHistData, 0, maxDistance);

        DoubleHistogram omsErrorHist = new DoubleHistogram("", omsHistErrorData, 0,
                maxDistance);
        omsErrorHist = new DoubleHistogram("", omsHistErrorData, 0, maxDistance);
        DoubleErrorHistogram errorHist = new DoubleErrorHistogram(omsHist,
                omsErrorHist);
        panel.addData(errorHist, new HistogramStyle());
        panel.setLabels("Distance (microns)", "Spikes/Run");
        panel.autoscale();
//		PlotUtilities.showData(errorHist, new HistogramStyle());

        return panel;
    }


    synchronized public PlotPanel makeOMSHistgrams(int[] neuronsInClass) {
        PlotPanel panel = new PlotPanel("omsHistograms");
        if (neuronsInClass.length > 0) {
            maxDistance = paramsFile.getDoubleCell(neuronsInClass[0], "omsMaxDistance");
            double[] omsHistData = paramsFile.getArrayCell(
                    neuronsInClass[0], "omsHistData");
            MeanVarianceCalculator[] mvc = new MeanVarianceCalculator[omsHistData.length];
            double[] average = new double[omsHistData.length];
            double[] error = new double[omsHistData.length];
            int[] count = new int[omsHistData.length];

            double[] omsDistances = new double[omsHistData.length];
            for (int j = 0; j < omsDistances.length; j++) {
                mvc[j] = new MeanVarianceCalculator();
                omsDistances[j] = (double) j * maxDistance / (double) omsHistData.length
                + maxDistance / 2.0 / (double) omsHistData.length;
            }
            double[] blankErrorHist = new double[omsDistances.length];
            ScatterPlotStyle style = new ScatterPlotStyle("OMS", SymbolType.NONE, 1,
                    Color.BLACK, true,
                    Color.black, 0.25f);
            for (int i = 0; i < neuronsInClass.length; i++) {
                int id = neuronsInClass[i];
                omsHistData = paramsFile.getArrayCell(
                        id, "omsHistData");
                DoubleList omsList = new DoubleList();
                DoubleList omsDist = new DoubleList();

                if (omsHistData != null) {
                    double sum = MathUtil.sumAbs(omsHistData);

                    for (int j = 0; j < omsHistData.length; j++) {
                        if (sum != 0) {
                            omsHistData[j] /= sum;
                        }
                        if (omsHistData[j] > 0.00001) {
                            count[j]++;
                            omsList.add(new Double(omsHistData[j]));
                            omsDist.add(new Double( (double) j * maxDistance /
                                    (double) omsHistData.length
                                    + maxDistance / 2.0 / (double) omsHistData.length));

                            mvc[j].add(omsHistData[j]);
                        }
                    }

                    ScatterPlot scatter = new ScatterPlot(omsDist.toArray(),
                            omsList.toArray(),
                            blankErrorHist);
                    panel.addData(scatter, style);

                }
            }

            for (int j = 0; j < average.length; j++) {
                average[j] = Double.isNaN(mvc[j].getMean()) ? 0.0 : mvc[j].getMean();
                error[j] = Double.isNaN(mvc[j].getMeanVariance()) ? 0.0 : mvc[j].getMeanVariance();
            }
            ScatterPlot scatter = new ScatterPlot(omsDistances, average, error);

            panel.addData(scatter, new ScatterPlotStyle("Average",
                    SymbolType.NONE, 2, Color.red, true, Color.red, 1));
            panel.setLabels("Distance (microns)", "Spikes/(Run*TotalSpikes)");
            panel.autoscale();
            double[] range = panel.getRange();
            panel.setRange(range[0], range[1], 0, range[3]);
        }
        return panel;

    }


//	This function is designed to be called externally, especially by GrivichThesis.
//	It calculates an average reference value/ an average OMS value.
//	The reference and OMS range is determined by the various Dist values.
    synchronized public static Num calculateOMSParameter(int[] neuronsInClass,
            ParametersFile paramsFile, double omsDistMin, double omsDistMax,
            double refDistMin, double refDistMax) {
        if (neuronsInClass.length > 0) {
            double maxDistance = paramsFile.getDoubleCell(neuronsInClass[0],
            "omsMaxDistance");
            double[] omsHistData = paramsFile.getArrayCell(
                    neuronsInClass[0], "omsHistData");
            MeanVarianceCalculator[] mvc = new MeanVarianceCalculator[omsHistData.length];
            double[] average = new double[omsHistData.length];
            double[] error = new double[omsHistData.length];
            int[] count = new int[omsHistData.length];

            double[] omsDistances = new double[omsHistData.length];
            for (int j = 0; j < omsDistances.length; j++) {
                mvc[j] = new MeanVarianceCalculator();
                omsDistances[j] = (double) j * maxDistance / (double) omsHistData.length
                + maxDistance / 2.0 / (double) omsHistData.length;
            }

            for (int i = 0; i < neuronsInClass.length; i++) {
                int id = neuronsInClass[i];
                omsHistData = paramsFile.getArrayCell(
                        id, "omsHistData");
                DoubleList omsList = new DoubleList();
                DoubleList omsDist = new DoubleList();

                if (omsHistData != null) {
                    double sum = MathUtil.sumAbs(omsHistData);

                    for (int j = 0; j < omsHistData.length; j++) {
                        if (sum != 0) {
                            omsHistData[j] /= sum;
                        }
                        if (omsHistData[j] > 0.00001) {
                            count[j]++;
                            omsList.add(new Double(omsHistData[j]));
                            omsDist.add(new Double( (double) j * maxDistance /
                                    (double) omsHistData.length
                                    + maxDistance / 2.0 / (double) omsHistData.length));
                            mvc[j].add(omsHistData[j]);
                        }
                    }
                }
            }

            for (int j = 0; j < average.length; j++) {
                average[j] = mvc[j].getMean();
                error[j] = mvc[j].getMeanVariance();
            }

            int omsDistCount = 0;
            int omsRefDistCount = 0;
            int earliestOmsDist = Integer.MAX_VALUE;
            int earliestRefDist = Integer.MAX_VALUE;
            for (int j = omsDistances.length - 1; j > -1; j--) {
                if (omsDistances[j] > omsDistMin && omsDistances[j] < omsDistMax) {
                    omsDistCount++;
                    earliestOmsDist = j;
                }
                if (omsDistances[j] > refDistMin && omsDistances[j] < refDistMax) {
                    omsRefDistCount++;
                    earliestRefDist = j;
                }
            }
            if (omsDistCount == 0) {
                throw new IllegalStateException(
                "There are no OMS values available for OMS parameter calculation.");
            }

            if (omsRefDistCount == 0) {
                throw new IllegalStateException(
                "There are no Reference values available for OMS parameter calculation.");
            }

            Num[] nums = new Num[omsDistCount];
            Num[] nums2 = new Num[omsRefDistCount];
            for (int i = earliestOmsDist; i < earliestOmsDist + omsDistCount; i++) {
                nums[i - earliestOmsDist] = new Num(average[i], error[i]);
            }
            for (int i = earliestRefDist; i < earliestRefDist + omsRefDistCount; i++) {
                nums2[i - earliestRefDist] = new Num(average[i], error[i]);
            }
            return Num.average(nums2).div(Num.average(nums));

        } else {
            return null;
        }

    }


}
