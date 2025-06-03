package edu.ucsc.neurobiology.vision.analysis;

import java.io.*;
import java.util.*;

import java.awt.*;

import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.util.*;

/**
 *
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class DriftingGratings {
    private NeuronFile[] neuronFile;
    private double[][][] stimulusInfo;
    public double[][] directionChoices; //[spatial:temporal:direction][values]
    private static double mmPerPixel = .0058;
    private static double frameTime = .0083;

    private static HistogramStyle histogramStyle = new HistogramStyle("Histogram Style",
            HistogramStyle.OutlineType.LINEAR, Color.black, .5f, false,
            Color.black, false, Color.black, .5f);


    public DriftingGratings(String[] sinusoidsFileName, String runDirectory) throws
    IOException {

        // open the sinusoids file
        neuronFile = new NeuronFile[sinusoidsFileName.length];
        stimulusInfo = new double[sinusoidsFileName.length][3][];

        for (int i = 0; i < neuronFile.length; i++) {
            if (sinusoidsFileName[i] != null &&
                    sinusoidsFileName[i].trim().length() != 0) {

                String sinFilePath =
                    runDirectory + File.separator + sinusoidsFileName[i] +
                    File.separator + sinusoidsFileName[i] + ".neurons";
                System.out.println("Sinusoids File: " + sinFilePath);
                neuronFile[i] = new NeuronFile(sinFilePath);
                String stimFilePath =
                    runDirectory + File.separator + "stimuli" + File.separator + "s" +
                    sinusoidsFileName[i].substring(5, 7) + ".txt";
                System.out.println("Sinusoids Stimulus File: " + stimFilePath);
                stimulusInfo[i] = loadDriftingSinusoidStimulus(stimFilePath);
            }
        }

        //Begin block to find and sort all spatial, temporal, direction values
        ArrayList<Double>[] directionInfoList = new ArrayList[3];
        for (int i = 0; i < 3; i++) {
            directionInfoList[i] = new ArrayList<Double>();
        }

        for (int i = 0; i < stimulusInfo.length; i++) {
            for (int j = 0; j < stimulusInfo[i].length; j++) {
                for (int k = 0; k < stimulusInfo[i][j].length; k++) {

                    boolean found = false;
                    for (int ii = 0; ii < directionInfoList[j].size(); ii++) {

                        if (stimulusInfo[i][j][k] ==
                            ( (Double) directionInfoList[j].get(ii)).doubleValue()) {
                            found = true;

                            break;
                        }

                    }
                    if (!found) {
                        directionInfoList[j].add(new Double(stimulusInfo[i][j][k]));
                    }
                }
            }
        }

        directionChoices = new double[3][];
        for (int i = 0; i < 3; i++) {
            directionChoices[i] = new double[directionInfoList[i].size()];
            for (int j = 0; j < directionInfoList[i].size(); j++) {
                directionChoices[i][j] = ( (Double) directionInfoList[i].get(j)).
                doubleValue();
            }
            Arrays.sort(directionChoices[i]);
        }
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < directionChoices[i].length; j++) {
                System.out.print(directionChoices[i][j] + " ");
            }
            System.out.println();
        }
    }


    synchronized public DoubleHistogram calculateSinResponseHist(double[]
                                                                        gratingSpikeRates) {
        int nDirections = directionChoices[2].length;
        double angleIncrement = 360.0 / nDirections;

        DoubleHistogram sinResponseHist =
            new DoubleHistogram("", -angleIncrement / 2.0, 360.0 - angleIncrement / 2.0,
                    angleIncrement);

        for (int i = 0; i < gratingSpikeRates.length; i = i + nDirections) {
            double spatialPeriod = directionChoices[0][i / nDirections /
                                                       directionChoices[1].length];
            double temporalPeriod = directionChoices[1][i / nDirections %
                                                        directionChoices[1].length];

            if (temporalPeriod >= 64.0 && spatialPeriod == 64.0) {
                for (int j = 0; j < nDirections; j++) {
                    sinResponseHist.fill(j * angleIncrement, gratingSpikeRates[i + j]);

                }

            }

        }

//		PlotUtil.showData("", makeSinusoidsPanel(sinResponseHist.toArray()));
        return sinResponseHist;
    }




    synchronized public double[] calculateDSParams(double[] gratingSpikeRates,
            double minSpatial,
            double maxSpatial, double minTemporal, double maxTemporal) {
        int nDirections = directionChoices[2].length;
        double angleIncrement = 360.0 / nDirections;

        DoubleHistogram sinResponseHist =
            new DoubleHistogram("", -angleIncrement / 2.0, 360.0 - angleIncrement / 2.0,
                    angleIncrement);

        for (int i = 0; i < gratingSpikeRates.length; i = i + nDirections) {
            double spatialPeriod = directionChoices[0][i / nDirections /
                                                       directionChoices[1].length];
            double temporalPeriod = directionChoices[1][i / nDirections %
                                                        directionChoices[1].length];

            if (temporalPeriod >= minTemporal && temporalPeriod <= maxTemporal &&
                    spatialPeriod >= minSpatial && spatialPeriod <= maxSpatial) {
                for (int j = 0; j < nDirections; j++) {
                    sinResponseHist.fill(j * angleIncrement, gratingSpikeRates[i + j]);
                }
            }

        }

        double[] sinResponse = sinResponseHist.toArray();
        double[] sinResponseIm = new double[sinResponse.length];

        FFT.fft(sinResponse, sinResponseIm, -1);

        if (sinResponse[0] * sinResponse[0] +
                sinResponseIm[0] * sinResponseIm[0] > 0) {
            //Find magnitude of f1 / "spike rate", or f0.
            double magDS = Math.sqrt( (sinResponse[1] * sinResponse[1] +
                    sinResponseIm[1] * sinResponseIm[1])
                    /
                    (sinResponse[0] * sinResponse[0] +
                            sinResponseIm[0] * sinResponseIm[0]));
            //Find phase, in radians, of f1, normalized by "spike rate", or f0.
            //0 and 2 Pi correspond to drifting grating, direction 0.
            double angDS = Math.atan(sinResponseIm[1] / sinResponse[1]);
            if (sinResponse[1] < 0) {
                angDS += Math.PI;
            }
            if (angDS >= Math.PI) {
                angDS -= 2 * Math.PI;

                //Find magnitude of f2 / "spike rate", or f0.
            }
            double magOS = Math.sqrt( (sinResponse[2] * sinResponse[2] +
                    sinResponseIm[2] * sinResponseIm[2])
                    /
                    (sinResponse[0] * sinResponse[0] +
                            sinResponseIm[0] * sinResponseIm[0]));
            //Find phase, in radians, of f2, normalized by "spike rate", or f0.
            //0 and 4*Pi correspond to drifting grating, direction 0.
            //2*Pi is drifting grating, direction Pi.
            double angOS = Math.atan(sinResponseIm[2] / sinResponse[2]);
            if (sinResponse[2] < 0) {
                angOS += Math.PI;
            }
            if (angOS >= Math.PI) {
                angOS -= 2 * Math.PI;

            }

            double[] dsParams = {magDS * Math.cos(angDS), magDS * Math.sin(angDS), magDS,
                    angDS,
                    magOS * Math.cos(angOS), magOS * Math.sin(angOS), magOS,
                    angOS};
            return dsParams;
        } else {
            double[] dsParams = {0, 0, 0, 0, 0, 0, 0, 0};
            return dsParams;
        }

    }


    public static double[] getDriftingSinusoidsResponse(
            int[] times, int[] ttl, double[] directions) {
        // remove the intermediary ttl's
//		System.err.println("TTLs " + ttl.length);
        IntegerList ttlList = new IntegerList();
        final int ttlsPerrun = ttl.length / directions.length;
//		System.out.println("ttlsPerRun: " + ttlsPerrun);
//		System.out.println("ttl 0: " + ttl[0]);
        for (int run = 0; run < directions.length; run++) {
            /// Chop first TTL
            ttlList.add(ttl[run * ttlsPerrun + 1]);
//			System.out.println("Run: " + run);
//			System.out.println("ttl[run*ttlsPerrun] " + ttl[run*ttlsPerrun]);
            /// Chop last TTL
            ttlList.add(ttl[ (run + 1) * ttlsPerrun - 2]);
//			System.out.println("ttl[ (run + 1) * ttlsPerrun - 1] " + ttl[ (run + 1) * ttlsPerrun - 1]);
        }
        ttl = ttlList.toArray();

        SpikeStream stream = new SpikeStream(times, ttl);

        // process the spikes
        double[] nSpikes = new double[directions.length];
        boolean count = false;
        int currentRun = -1;

        int t;
        while ( (t = stream.getNext()) != Integer.MAX_VALUE) {
            if (t < 0) {
                count = !count;
                if (count == true) {
                    currentRun++;
                }
                if (currentRun >= directions.length) {
                    break;
                }
            } else {
                if (count) {
                    nSpikes[currentRun]++;
                }
            }
        }

        return nSpikes;
    }
    
    synchronized public double[] calculateMeansAndRms(double[] gratingSpikeRates) {



        double lowSpatialFrequency = Math.log(1 /
                (mmPerPixel *
                        (directionChoices[0][directionChoices[0].length - 1]))) / Math.log(10);
        double lowTemporalFrequency = Math.log(1 /
                (frameTime * (directionChoices[1][directionChoices[1].length - 1]))) /
                Math.log(10);
        double spatialBinSize = Math.log(1 /
                (mmPerPixel *
                        directionChoices[0][directionChoices[0].
                                            length - 2])) / Math.log(10) -
                                            lowSpatialFrequency;
        double temporalBinSize = Math.log(1 /
                (frameTime *
                        directionChoices[1][directionChoices[1].
                                            length - 2])) / Math.log(10) -
                                            lowTemporalFrequency;
//		double highSpatialFrequency = lowSpatialFrequency +
//		directionChoices[0].length * spatialBinSize;
//		double highTemporalFrequency = lowTemporalFrequency +
//		directionChoices[1].length * temporalBinSize;
        double[][] bins = new double[directionChoices[0].length] //[spatial][temporal]
                                     [directionChoices[1].length];




        for (int j = 0; j < gratingSpikeRates.length; j++) {
            bins[directionChoices[0].length - 1 -
                 j / directionChoices[1].length][directionChoices[1].length - 1 -
                                                 j % directionChoices[1].length]
                                                 += gratingSpikeRates[j];
        }

        double spatialMean = 0.0;
        double temporalMean = 0.0;
        double spatialRMS = 0.0;
        double temporalRMS = 0.0;
        double norm = 0.0;
        double temp = 0.0;

        for(int i=0; i<bins.length; i++) {
            for(int j=0; j<bins[0].length; j++) {
                spatialMean += (spatialBinSize*i+lowSpatialFrequency)*bins[i][j];
                temporalMean += (temporalBinSize*j+lowTemporalFrequency)*bins[i][j];
                norm += bins[i][j];
            }
        }
        spatialMean/=norm;
        temporalMean/=norm;
        
        for(int i=0; i<bins.length; i++) {
            for(int j=0; j<bins[0].length; j++) {
                temp = spatialBinSize*i+lowSpatialFrequency-spatialMean;
                spatialRMS += temp*temp*bins[i][j];
                
                temp = temporalBinSize*j+lowTemporalFrequency-temporalMean;
                temporalRMS += temp*temp*bins[i][j];
            }
        }
        
        spatialRMS = Math.sqrt(spatialRMS/norm);
        temporalRMS = Math.sqrt(temporalRMS/norm);
        


        return new double[] {spatialMean, temporalMean, spatialRMS, temporalRMS};

    }


    public static double[][] loadDriftingSinusoidStimulus(String fileName) throws
    IOException {

        InputStreamReader r = new InputStreamReader(new FileInputStream(fileName));
        StreamTokenizer st = new StreamTokenizer(r);

        st.whitespaceChars('(', '(');
        st.whitespaceChars(')', ')');
        st.whitespaceChars(':', ':');
        st.whitespaceChars('#', '#');
        st.eolIsSignificant(true);

//		while (true) {
//		st.nextToken();
//		System.out.println(st.ttype);
//		}

        // read the stimulus type
        st.nextToken();
        if (st.ttype != StreamTokenizer.TT_WORD || !st.sval.equals("TYPE")) {
            throw new IOException("Missing TYPE.");
        }
        st.nextToken();
        if (st.ttype != StreamTokenizer.TT_WORD) {
            throw new IOException("Wrong TYPE.");
        }
        String type = st.sval;
        System.out.println("TYPE : " + type);

        // skip to the end of line
        do {
            st.nextToken();
        } while (st.ttype != StreamTokenizer.TT_EOL);

        // read the run lines
        DoubleList directions = new DoubleList();
        DoubleList spatialPeriods = new DoubleList();
        DoubleList temporalPeriods = new DoubleList();

        while (true) {
            // skip spatial period and temporal period
            st.nextToken();
            if (st.ttype != StreamTokenizer.TT_WORD || !st.sval.equals("SPATIAL-PERIOD")) {
                throw new IOException("Missing SPATIAL-PERIOD at line " + st.lineno());
            }
            st.nextToken();
            if (st.ttype != StreamTokenizer.TT_NUMBER) {
                throw new IOException("Missing SPATIAL-PERIOD at line " + st.lineno());
            }
            double spatialPeriod = st.nval;

            st.nextToken();
            if (st.ttype != StreamTokenizer.TT_WORD || !st.sval.equals("TEMPORAL-PERIOD")) {
                throw new IOException("Missing TEMPORAL-PERIOD at line " + st.lineno());
            }
            st.nextToken();
            if (st.ttype != StreamTokenizer.TT_NUMBER) {
                throw new IOException("Missing TEMPORAL-PERIOD at line " + st.lineno());
            }
            double temporalPeriod = st.nval;

            // read the direction of the drift
            st.nextToken();
            if (st.ttype != StreamTokenizer.TT_WORD || !st.sval.equals("DIRECTION")) {
                throw new IOException("Missing DIRECTION at line " + st.lineno());
            }
            st.nextToken();
            if (st.ttype != StreamTokenizer.TT_NUMBER) {
                throw new IOException("Missing DIRECTION at line " + st.lineno());
            }
            double direction = st.nval;

            // read the EOL
            st.nextToken();
            if (st.ttype != StreamTokenizer.TT_EOL) {
                throw new IOException("Too many values at line " + st.lineno());
            }

            directions.add(direction);
            spatialPeriods.add(spatialPeriod);
            temporalPeriods.add(temporalPeriod);

            // see whether we reached the end of the file
            st.nextToken();
            if (st.ttype == StreamTokenizer.TT_EOF) {
                break;
            } else {
                st.pushBack();
            }
        }

        r.close();
        System.out.println("directions " + directions.size());

        double[][] out = new double[3][];
        out[0] = spatialPeriods.toArray();
        out[1] = temporalPeriods.toArray();
        out[2] = directions.toArray();
        return out;
    }


    synchronized public double[][] calculate(int id) throws IOException {
        double[] gratingSpikeRates = new double[directionChoices[0].length *
                                                directionChoices[1].length *
                                                directionChoices[2].length];
        double[] gratingSpikeRatesAveraged = new double[directionChoices[0].length *
                                                        directionChoices[1].length];
        for (int i = 0; i < neuronFile.length; i++) {
            //directionInfo is used to determine number of rows only
            double[] sinResponse = getDriftingSinusoidsResponse(
                    neuronFile[i].getSpikeTimes(id),
                    neuronFile[i].getTTLTimes(), stimulusInfo[i][0]);
            for (int j = 0; j < sinResponse.length; j++) {
                for (int spatial = 0; spatial < this.directionChoices[0].length;
                spatial++) {
                    for (int temporal = 0; temporal < directionChoices[1].length;
                    temporal++) {
                        for (int direction = 0; direction < directionChoices[2].length;
                        direction++) {
                            if (stimulusInfo[i][0][j] == directionChoices[0][spatial] &&
                                    stimulusInfo[i][1][j] == directionChoices[1][temporal] &&
                                    stimulusInfo[i][2][j] == directionChoices[2][direction]) {
                                gratingSpikeRates[spatial * directionChoices[1].length *
                                                  directionChoices[2].length
                                                  + temporal * directionChoices[2].length
                                                  + direction] += sinResponse[j];
                                gratingSpikeRatesAveraged[spatial *
                                                          directionChoices[1].length + temporal]
                                                          += (double) sinResponse[j] /
                                                          (double) directionChoices[2].length;
                            }
                        }
                    }
                }
//				hist[ (int) directionInfo[i][2][j] /
//				45].fill(Math.log(directionInfo[i][0][j]) / Math.log(2),
//				Math.log(directionInfo[i][1][j]) / Math.log(2),
//				sinResponse[j]);
            }

        }

        /*       DoubleHistogram2D hist = new DoubleHistogram2D("", 0,
                   directionChoices[0].length, 0, directionChoices[1].length, 1, 1);
               for (int j = 0; j < gratingSpikeRatesAveraged.length; j++) {
                   hist.fill(j / directionChoices[1].length,
                             j % directionChoices[1].length,
                             gratingSpikeRatesAveraged[j]);
               }*/
        double[][] directionRates = new double[2][];
        directionRates[0] = gratingSpikeRates;
        directionRates[1] = gratingSpikeRatesAveraged;
        return directionRates;
    }


    synchronized static public PlotPanel makeGratingPanel(
            int neuron, ParametersFile paramsFile, boolean showLabels) {

        double[] gratingSpikeRates = paramsFile.getArrayCell(
                neuron, "gratingResponse");
        int dividers = paramsFile.getArrayCell(
                neuron, "spatialFrequencies").length - 1;
        return makeGratingPanel(gratingSpikeRates, showLabels, dividers);
    }


    synchronized static public PlotPanel makeGratingPanel(double[] gratingSpikeRates,
            boolean showLabels, int dividers) {

        double[] frequencyHist = new double[gratingSpikeRates.length];
        for (int i = 0; i < frequencyHist.length; i++) {
            frequencyHist[i] = gratingSpikeRates[frequencyHist.length - 1 - i];
        }
        DoubleHistogram hist = new DoubleHistogram(
                "histogram", frequencyHist, 0, gratingSpikeRates.length);

//		PlotUtil.showData( "Grating Panel", hist,
//		new HistogramStyle(), 400, 400
//		);

        PlotPanel panel = new PlotPanel("GratingResponse");
        panel.addData(hist, histogramStyle);
        if (showLabels) {
            panel.setLabels("Run Number", "Spikes/Run");
        }

        makeDashedLineDividers(dividers, panel, 0, gratingSpikeRates.length, 0,
                MathUtil.max(gratingSpikeRates));
        panel.autoscale();

        return panel;
    }


    synchronized static public PlotPanel makeAverageGratingPanel(int neuron,
            ParametersFile paramsFile, boolean showLabels) {
        double[] gratingSpikeRatesAveraged = paramsFile.getArrayCell(
                neuron, "averageGrating");

        int dividers = paramsFile.getArrayCell(
                neuron, "spatialFrequencies").length - 1;

        double[] frequencyHist = new double[gratingSpikeRatesAveraged.length];
        for (int i = 0; i < frequencyHist.length; i++) {
            frequencyHist[i] = gratingSpikeRatesAveraged[frequencyHist.length - 1 - i];
        }
        DoubleHistogram hist = new DoubleHistogram(
                "gratings", frequencyHist, 0, gratingSpikeRatesAveraged.length);

//		PlotUtilities.showData( "Autocorrelation Plot", hist,
//		new HistogramStyle(), new Rectangle(400,0,400,400)
//		);
        PlotPanel panel = new PlotPanel("AverageGrating");
        panel.addData(hist, histogramStyle);
        if (showLabels) {
            panel.setLabels("Run Number", "Spikes/Run");
        }

        makeDashedLineDividers(dividers, panel, 0, gratingSpikeRatesAveraged.length, 0,
                MathUtil.max(gratingSpikeRatesAveraged));
        panel.autoscale();
        return panel;
    }


    synchronized static public void makeDashedLineDividers(int dividers, PlotPanel panel,
            double minX, double maxX, double minY, double maxY) {

        double spacing = (maxX - minX) / (dividers + 1.0);
        for (int i = 0; i < dividers; i++) {
            double[] xPoints = new double[] {minX + spacing * (i + 1),
                    minX + spacing * (i + 1)};
            double[] yPoints = new double[] {minY, maxY};

            ScatterPlot scatter = new ScatterPlot(xPoints, yPoints, null);
            ScatterPlotStyle style = new ScatterPlotStyle("Divider", SymbolType.NONE, 0,
                    Color.black, true, Color.black, .5f);

            style.setDashPattern("5,5");
            panel.addData(scatter, style);

        }

    }



    synchronized static public PlotPanel makeAverageGratingPanel2D(int neuron,
            ParametersFile paramsFile, boolean showLabels) {
        return makeClassAverageGratingPanel2D(new int[] {neuron}, paramsFile, showLabels);

    }


    synchronized static public PlotPanel makeClassAverageGratingPanel2D(int[]
                                                                            neuronsInClass,
                                                                            ParametersFile paramsFile, boolean showLabels) {
        if (neuronsInClass.length > 0) {
            double[][] directionChoices = new double[3][];
            directionChoices[0] = paramsFile.getArrayCell(
                    neuronsInClass[0], "spatialFrequencies");
            directionChoices[1] = paramsFile.getArrayCell(
                    neuronsInClass[0], "temporalFrequencies");
            directionChoices[2] = paramsFile.getArrayCell(
                    neuronsInClass[0], "directions");
            double lowSpatialFrequency = Math.log(1 /
                    (mmPerPixel *
                            (directionChoices[0][directionChoices[0].length - 1]))) / Math.log(10);
            double lowTemporalFrequency = Math.log(1 /
                    (frameTime * (directionChoices[1][directionChoices[1].length - 1]))) /
                    Math.log(10);
            double spatialBinSize = Math.log(1 /
                    (mmPerPixel *
                            directionChoices[0][directionChoices[0].
                                                length - 2])) / Math.log(10) -
                                                lowSpatialFrequency;
            double temporalBinSize = Math.log(1 /
                    (frameTime *
                            directionChoices[1][directionChoices[1].
                                                length - 2])) / Math.log(10) -
                                                lowTemporalFrequency;
            double highSpatialFrequency = lowSpatialFrequency +
            directionChoices[0].length * spatialBinSize;
            double highTemporalFrequency = lowTemporalFrequency +
            directionChoices[1].length * temporalBinSize;
            double[][] bins = new double[directionChoices[0].length]
                                         [directionChoices[1].length];
            for (int i = 0; i < neuronsInClass.length; i++) {

                double[] gratingSpikeRatesAveraged = paramsFile.getArrayCell(
                        neuronsInClass[i], "averageGrating");
                directionChoices = new double[3][];
                directionChoices[0] = paramsFile.getArrayCell(
                        neuronsInClass[i], "spatialFrequencies");
                directionChoices[1] = paramsFile.getArrayCell(
                        neuronsInClass[i], "temporalFrequencies");
                directionChoices[2] = paramsFile.getArrayCell(
                        neuronsInClass[i], "directions");

                for (int j = 0; j < gratingSpikeRatesAveraged.length; j++) {
                    bins[directionChoices[0].length - 1 -
                         j / directionChoices[1].length][directionChoices[1].length - 1 -
                                                         j % directionChoices[1].length]
                                                         += gratingSpikeRatesAveraged[j];
                }
            }

            DoubleHistogram2D hist = new DoubleHistogram2D("", lowSpatialFrequency,
                    highSpatialFrequency, lowTemporalFrequency, highTemporalFrequency, bins);

//			for (int i = 0; i < hist.getBinCountX(); i++) {
//			for (int j = 0; j < hist.getBinCountY(); j++) {
//			System.out.println("X: " + i + " Y: " + j + " : " + hist.getBin(i, j));
//			}
//			}
//			PlotUtilities.showData(id + "" , hist,  new HistogramStyle(), new Rectangle(800,0,400,400));

            PlotPanel panel = new PlotPanel("AverageGratingPanel2D");
            panel.addData(hist, new HistogramStyle());
            if (showLabels == true) {
                panel.setLabels("<html>log<sub>10</sub> Spatial Frequency(cyc/mm)",
                "<html>log<sub>10</sub> Temporal F.(cyc/s)");
            }
            panel.getXAxis().setFixedTickSpacing(lowSpatialFrequency,
                    2*spatialBinSize + .00001,
                    2);

            panel.getYAxis().setFixedTickSpacing(lowTemporalFrequency,
                    2*temporalBinSize + .00001, 2);
            panel.autoscale();
            return panel;
        } else {
            return new PlotPanel();
        }
    }


    synchronized static public PlotPanel makeSinusoidsPanel(int[] neuronsInClass,
            ParametersFile paramsFile) {
        if (neuronsInClass.length > 0) {
            double[] directions = paramsFile.getArrayCell(neuronsInClass[0], "directions");
            PlotPanel panel = new PlotPanel("overlayedSinusoids");

            for (int id = 0; id < neuronsInClass.length; id++) {
                double[] spikes = (double[]) paramsFile.getArrayCell(neuronsInClass[id],
                        "gratingResponse");
                double[] sinResponse = new double[directions.length];
                int spikeCount = 0;
                for (int i = 0; i < spikes.length; i++) {
                    sinResponse[i % directions.length] += spikes[i];
                    spikeCount += spikes[i];
                }
                for (int i = 0; i < directions.length; i++) {
                    sinResponse[i] /= spikeCount;
                }
                //currently no rotation.  
                double[] sinResponseRotated = new double[sinResponse.length];
                for (int i = 0; i < sinResponse.length; i++) {
                    sinResponseRotated[i] = sinResponse[ (i + 0) % sinResponse.length];
                }
                ScatterPlot h = new ScatterPlot(directions, sinResponseRotated,
                        new double[sinResponse.length]);

                panel.addData(h,
                        new ScatterPlotStyle("", SymbolType.NONE, 1, Color.black, true,
                                Color.black, .25f));

//				DoubleHistogram hist = new DoubleHistogram("Sinusoid Response", sinResponse, 0 , 360-360/(directions.length-1));
//				panel.addData(hist, new HistogramStyle(HistogramStyle.CoordinatesType.POLAR));




            }
            panel.setLabels("Direction (degrees)", "Normalized Spike Rate");
            panel.autoscale();
            return panel;
        } else {
            return new PlotPanel();
        }
    }


    synchronized static public PlotPanel makeSinusoidsPanel(int id,
            ParametersFile paramsFile) {
        double[] sinResponse = paramsFile.getArrayCell(id, "sinResponse");
        return makeSinusoidsPanel(sinResponse);
    }


    //sinResponse has been removed from the params file.  Display functions kept
    //for debugging purposes.
    synchronized static public PlotPanel makeSinusoidsPanel(double[] sinResponse) {

//		DoubleHistogram h = new DoubleHistogram("", -2 * 11.25, 360 + 2 * 11.25,
//		2 * 11.25);
//		for (int i = 0; i < directions.length; i++) {
//		h.fill(directions[i], nSpikes[i]);
//		}

        DoubleHistogram h = new DoubleHistogram("Sinusoid Response", sinResponse, 0, 15);

        PlotPanel p = new PlotPanel();
        HistogramStyle style = new HistogramStyle();
//		style.setCoordinatesType(HistogramStyle.POLAR);
        p.addData(h, style);
        p.autoscale();
        double[] range = p.getRange();
//		System.out.println(range[0] + " " + range[1] + " " + range[2] + " " + range[3]);
        range[2] = 0;
        p.setRange(range);

        return p;
    }


}
