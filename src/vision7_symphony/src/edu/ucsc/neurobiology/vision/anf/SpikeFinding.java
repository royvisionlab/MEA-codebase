package edu.ucsc.neurobiology.vision.anf;

import java.io.*;
import java.text.*;
import java.util.*;

import javax.swing.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.calculations.*;
import edu.ucsc.neurobiology.vision.electrodemap.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.io.chunk.ChunkFile;
import edu.ucsc.neurobiology.vision.io.chunk.GlobalsFile;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.onlinediagnostics.*;
import edu.ucsc.neurobiology.vision.parameters.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * Vision calculation that does spike finding on the raw data. It can be configured to
 * save the spikes and to calculate and save the covariance matrices of each electrode.
 * It is also used during online data taking when it is configured to display a set of
 * diagnostic plots.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class SpikeFinding extends AbstractCalculation {
    public static enum AnalysisToDo { DO_NOTHING, SAVE_SPIKES, SAVE_SPIKES_AND_COVARINCES }

    private static int MINCOVSPIKES = 100;

    private float[] sigma;
    private MultipleCompressedSampleInputStream sampleInputStream;
    private java.util.Timer timer;
    private SpikeFinder spikeFinder;
    private SpikeBuffer spikeBuffer;
    private NoiseFinder noiseFinder;
    private SpikeSaver spikeSaver;
    private JPanel diagnosticPanel;
    private JLabel refreshRateLabel, samplesReadLabel, spikesFoundLabel, spikeRateLabel,
    spikeBufferLabel;
    private JProgressBar inputBufferUsage;
    private JProgressBar[] outputBuffer;
    private RawDataSaver dataSaver;
    private Vision application;
    private int nElectrodes;
    private CovariancesCalculator covarianceCalculator, noiseCovarianceCalculator;
    private int nSamples;

    // diagnostic plots
    OscilloscopeWindow oscilloscope;
    MeansAndSigmasDiagnosticPlot meansAndSigmas;
    SpikeSummary spikeSummary;
    SpikeStatisticsPlot spikeRatePlot;
    TTLPeriodPlot ttlPlot;
    boolean leaveWindowsOpen;

    private boolean setElectrodes = false;
    //default values used when setElectrodes is false.  packedArrayID is overridden with information from raw data file.
    private boolean flipX = false, flipY = false;
    private int arrayID, arrayPart, arrayNParts;
    private int packedArrayID;

    public final static int SPIKE_BUFFER_SIZE = 60;


    public SpikeFinding() {
        if (Vision.isGUIBased()) {
            diagnosticPanel = new JPanel(new TableLayout(0, 0));
        }
    }


    TimerTask timerTask = new TimerTask() {
        NumberFormat f = DecimalFormat.getInstance();

        public void run() {
            int samplesRead = sampleInputStream.getSamplesRead();
            application.setProgress( (int) (100.0 * samplesRead / (nSamples - 1)));

            // show GUI info if appropriate
            if (Vision.isGUIBased()) {
                refreshRateLabel.setText(StringUtil.format(spikeFinder.getRefreshRate(),
                        0));
                samplesReadLabel.setText(f.format(sampleInputStream.getSamplesRead()));
                spikesFoundLabel.setText(f.format(spikeFinder.nSpikes));
                spikeRateLabel.setText(f.format(
                        (int) (spikeFinder.nSpikes / (samplesRead * 5e-5) / nElectrodes)));
                if (spikeBuffer != null) {
                    spikeBufferLabel.setText("" + spikeBuffer.spikeList.size());
                }
                inputBufferUsage.setValue(sampleInputStream.getBufferUsage());
                if (dataSaver != null) {
                    for (int i = 0; i < outputBuffer.length; i++) {
                        outputBuffer[i].setValue( (int) dataSaver.getBufferUsage(i));
                    }
                }
            }

            // finish the calculation
            if (sampleInputStream.isFinished()) {
                timer.cancel();
                application.endProgressBar();

                endCalculation();
                Vision.getInstance().getCalculationManager().calculationDone();
            }
        }
    };


    public void endCalculation() {
        if (!leaveWindowsOpen) {
            if (meansAndSigmas != null) {
                meansAndSigmas.dispose();
            }

            if (oscilloscope != null) {
                oscilloscope.dispose();
            }

            if (spikeSummary != null) {
                spikeSummary.dispose();
            }

            if (spikeRatePlot != null) {
                spikeRatePlot.dispose();
            }

            if (ttlPlot != null) {
                ttlPlot.dispose();
            }
        }
    }


    public void startCalculation() {
        application = Vision.getInstance();
        application.sendMessage("Data Acquisition...");

        application.startProgressBar();

        // start the machinery !
        sampleInputStream.start();
        if (covarianceCalculator != null) {
            covarianceCalculator.start();
        }
        if (noiseCovarianceCalculator != null) {
            noiseCovarianceCalculator.start();
        }

        // start the timer
        timer = new java.util.Timer(true);
        timer.scheduleAtFixedRate(timerTask, 0, 100);
    }


    public JComponent getDiagnosticPanel() {
        return diagnosticPanel;
    }


    public static float[] getSigmas(String fileNameOrValue, int nElectrodes) throws
    IOException {

        float[] sigma = new float[nElectrodes];
        try {
            float s = Float.parseFloat(fileNameOrValue);
            Arrays.fill(sigma, s);
        } catch (NumberFormatException e) {
            double[] s = null;
            try {
                s = IOUtil.loadDoubleArray(fileNameOrValue);
            } catch (IOException ex) {
                throw new IOException(
                        "The Spike Threshold parameter: " + fileNameOrValue +
                ", is not a number nor a file name");
            }

            if (s.length != nElectrodes) {
                throw new Error("s.length != nElectrodes");
            }
            for (int i = 0; i < nElectrodes; i++) {
                sigma[i] = (float) s[i];
            }
        }
        sigma[0] = 100;

        return sigma;
    }


    public void setParameters(HashMap<String, String> parameters) throws Exception {
        application = Vision.getInstance();
        HashMap<String, String> p = parameters;

        String rawDataSource = p.get("Raw_Data_Source");
        float spikeThreshold = Float.parseFloat(p.get("Spike Threshold"));
        float ttlThreshold = Float.parseFloat(p.get("TTL Threshold"));
        float meanTimeConstant = (float) Double.parseDouble(p.get("Mean Time Constant"));
        int bufferSizeInBytes = 1024 * Integer.parseInt(p.get("Buffer Size (Kb)"));
        int nBuffers = Integer.parseInt(p.get("Buffers Count"));
        String st = p.get("Sigma");

        int nSamplesToBuffer = bufferSizeInBytes / 770;
        bufferSizeInBytes = nSamplesToBuffer * 770;

        //        System.out.println("Sample Size: " + 770 + " bytes");
        //        System.out.println("Corrected Buffer Size: " + bufferSizeInBytes + " bytes");
        //        System.out.println("Input/Output Buffer Size: " + StringUtil.format(
        //            bufferSizeInBytes * nBuffers / 1024.0 / 1024.0, 1) + " Mb");

        // create the Sample Input Stream
        application.sendMessage("Spike Finding: create the Sample Input Stream");
        boolean waitForData = parameters.containsKey("waitForData") ? Boolean.valueOf(parameters.get("waitForData")) : false;
        sampleInputStream = new MultipleCompressedSampleInputStream(rawDataSource, bufferSizeInBytes, nBuffers, waitForData);
        RawDataHeader512 header = sampleInputStream.getHeader();
        nSamples = header.getNumberOfSamples();
        
        //data in raw header is not always trustworthy.  setElectrodes overrides it.
        setElectrodes = Boolean.valueOf((String) parameters.get("Set Electrodes"));
        if (setElectrodes) {
            arrayID = Integer.valueOf((String) parameters.get("Set Electrodes.arrayID"));
            arrayPart = Integer.valueOf((String) parameters.get("Set Electrodes.arrayPart"));
            arrayNParts = Integer.valueOf((String) parameters.get("Set Electrodes.arrayNParts"));
            flipX = Boolean.valueOf( (String) parameters.get("Set Electrodes.flipX")).booleanValue();
            flipY = Boolean.valueOf( (String) parameters.get("Set Electrodes.flipY")).booleanValue();
            packedArrayID = arrayID + (arrayPart << 16) + (arrayNParts << 24);
        } else {
            packedArrayID = header.getArrayID();
            arrayID = packedArrayID & 0xFFFF;
            arrayPart = (packedArrayID >> 16) & 0xFF;
            arrayNParts = packedArrayID >> 24;
            if (arrayPart == 0) {
                arrayPart = 1;
                arrayNParts = 1;
            }
        }
        
        ElectrodeMap electrodeMap = ElectrodeMapFactory.getElectrodeMap(packedArrayID);
        nElectrodes = electrodeMap.getNumberOfElectrodes();
        sigma = getSigmas(st, nElectrodes);
        MathUtil.multiply(sigma, spikeThreshold);
        
        // RAW DATA SAVING
        boolean saveRawData = Boolean.valueOf(p.get("saveRawData"));
        if (saveRawData) {
            String commonPath = (p.get("saveRawData.Common Path")).trim();
            String outputFolders = p.get("saveRawData.outputFolders");
            int secondsToStream = Integer.parseInt(p.get("saveRawData.Stop Streaming After (sec)"));
            
            // create the Data Saver and connect it with the Sample Input Stream (if needed)
            if (header.getDatasetIdentifier().equals("9999-99-99-data999")) {
                System.out.println("WILL NOT SAVE/STREAM DATASET " + header.getDatasetIdentifier());
                saveRawData = false;
            } else {
                String[] fileNames = StringUtil.decomposeString(outputFolders, ";");
                if (fileNames.length != 1 && fileNames.length != 2 &&
                        fileNames.length != 4 && fileNames.length != 6 &&
                        fileNames.length != 8) {
                    throw new IllegalArgumentException(
                    "Can split raw data in 1, 2, 4, 6 or 8 files only");
                }
                
                dataSaver = new RawDataSaver(fileNames, commonPath, header,
                        nSamplesToBuffer, nBuffers, secondsToStream);
                sampleInputStream.addSampleListener(dataSaver);
            }
        }

        // create the Spike Finder
        application.sendMessage("Spike Finding: create the Spike Finder");
        spikeFinder = new SpikeFinder(electrodeMap, sigma, ttlThreshold, meanTimeConstant);
        sampleInputStream.addSampleListener(spikeFinder);
        sampleInputStream.setInitializationListener(spikeFinder);

        // ANALYSIS
        boolean analysis = Boolean.valueOf(p.get("Analysis"));
        if (analysis) {
            AnalysisToDo analysisToDo =	AnalysisToDo.values()[ (int) Double.parseDouble(p.get("Analysis.Analysis To Do"))];
            final int nLeftPoints = Integer.parseInt(p.get("Analysis.Left Points"));
            final int nRightPoints = Integer.parseInt(p.get("Analysis.Right Points"));
            final double minimizationError = Double.parseDouble(p.get("Analysis.Minimization Error"));
            final int spikesToUse = Integer.parseInt(p.get("Analysis.Spike To Use"));
            final int minimumNoiseEvents = Integer.parseInt(p.get("Analysis.Minimum Noise Events"));
            ElectrodeUsage electrodeUsage =	ElectrodeUsage.values()[ (int) Double.parseDouble(p.get("Analysis.Electrode Usage"))];
            String outputPath = new File(p.get("Analysis.Output_Path")).getAbsolutePath();

            File outputDir = new File(outputPath);
            SpikeFile spikeFile = null;

            // Save spikes?
            if (analysisToDo == AnalysisToDo.SAVE_SPIKES || analysisToDo == AnalysisToDo.SAVE_SPIKES_AND_COVARINCES) {

                spikeBuffer = new SpikeBuffer(SPIKE_BUFFER_SIZE);
                spikeFinder.addSpikeListener(spikeBuffer);
                sampleInputStream.addSampleListener(spikeBuffer);

                // create the output dir
                outputDir.mkdirs();
                String outFileName =
                    outputDir.getAbsolutePath() + File.separator + outputDir.getName();

                // write the header out
                spikeFile = new SpikeFile(outFileName + ".spikes", packedArrayID,
                        meanTimeConstant, spikeThreshold,
                        header.getNumberOfSamples(),
                        header.getSamplingFrequency());

                spikeSaver = new SpikeSaver(spikeFile);
                spikeBuffer.addSpikeListener(spikeSaver);

                GlobalsFile globalsFile = new GlobalsFile(outFileName + ".globals", ChunkFile.WRITE);
                globalsFile.writeRDH512(header);
                globalsFile.setImageCalibrationParams(electrodeMap.micronsPerPixelX, electrodeMap.micronsPerPixelY, 
                        electrodeMap.centerX, electrodeMap.centerY, flipX, flipY, electrodeMap.angle, arrayID, arrayPart, arrayNParts);
                globalsFile.close();
            }

            if (analysisToDo == AnalysisToDo.SAVE_SPIKES_AND_COVARINCES) {
                String covFileName = outputDir.getAbsolutePath() + File.separator +
                outputDir.getName() + ".cov";

                covarianceCalculator = new CovariancesCalculator(
                        covFileName, nLeftPoints, nRightPoints, electrodeUsage, MINCOVSPIKES,
                        spikesToUse, minimizationError, spikeFile.getHeader(), true);

                spikeBuffer.addSpikeListener(covarianceCalculator);
                sampleInputStream.addSampleListener(covarianceCalculator);
                
                // Whitened?
                if (minimumNoiseEvents > 0) {
                    noiseFinder = new NoiseFinder(nElectrodes, nSamples, minimumNoiseEvents, nLeftPoints+nRightPoints, SPIKE_BUFFER_SIZE);
                    String noiseCovFileName = outputDir.getAbsolutePath() + File.separator +
                    outputDir.getName() + ".ncov";

                    noiseCovarianceCalculator = new CovariancesCalculator(
                            noiseCovFileName, nLeftPoints, nRightPoints, electrodeUsage, MINCOVSPIKES,
                            spikesToUse, minimizationError, spikeFile.getHeader(), false);
                    spikeBuffer.addSpikeListener(noiseFinder); //gives spikes to the noise finder
                    noiseFinder.addNoiseListener(noiseCovarianceCalculator); //gives noise to the the covariance calculator.
                    sampleInputStream.addSampleListener(noiseCovarianceCalculator);
                    sampleInputStream.addSampleListener(noiseFinder);
                }
            }
        }
        
        // DIAGNOSTIC PLOTS
        final boolean diagnosticPlots = Boolean.valueOf(p.get("Diagnostic Plots"));
        if (diagnosticPlots) {
            boolean showOscilloscope = Boolean.valueOf(p.get(
            "Diagnostic Plots.Show Oscilloscope"));
            int nTraces = Integer.parseInt(p.get(
            "Diagnostic Plots.Number of Traces"));
            boolean showSpikeRate = Boolean.valueOf(p.get(
            "Diagnostic Plots.Show Spike Rate"));
            double spikeRateBinning = Double.parseDouble(p.get(
            "Diagnostic Plots.Spike_Rate_Binning"));
            boolean findMeansAndSigmas = Boolean.valueOf(p.get(
            "Diagnostic Plots.Find_Means_Sigmas"));
            boolean showSpikeSummary = Boolean.valueOf(p.get(
            "Diagnostic Plots.Show_Spike_Summary"));
            boolean showTTLPeriodPlot = Boolean.valueOf(p.get(
            "Diagnostic Plots.Show_TTL_Period_Plot"));
            leaveWindowsOpen = Boolean.valueOf(p.get(
            "Diagnostic Plots.Leave Windows Open"));
            
            if (findMeansAndSigmas) {
                meansAndSigmas = new MeansAndSigmasDiagnosticPlot(sampleInputStream,
                        electrodeMap.getNumberOfElectrodes());
            }

            if (showOscilloscope) {
                oscilloscope = new OscilloscopeWindow(sampleInputStream, header, nTraces);
            }

            // create the spike summary plots
            if (showSpikeSummary) {
                spikeSummary = new SpikeSummary(spikeFinder, electrodeMap);
            }

            if (showSpikeRate) {
                spikeRatePlot = new SpikeStatisticsPlot(
                        spikeFinder, electrodeMap, header.getTimePerSample(),
                        spikeRateBinning, nSamples / (double) header.getSamplingFrequency());
            }

            if (showTTLPeriodPlot) {
                ttlPlot = new TTLPeriodPlot(spikeFinder);
            }
        }



        if (Vision.isGUIBased()) {
            // create the diagnostic panel
            diagnosticPanel.add(new JLabel("Input Buffer:"));
            inputBufferUsage = new JProgressBar(JProgressBar.HORIZONTAL, 0, 100);
            diagnosticPanel.add(inputBufferUsage);

            diagnosticPanel.add(new JLabel("Samples read:"));
            samplesReadLabel = new JLabel("0", JLabel.RIGHT);
            diagnosticPanel.add(samplesReadLabel);

            diagnosticPanel.add(new JLabel("Spikes found:"));
            spikesFoundLabel = new JLabel("0", JLabel.RIGHT);
            diagnosticPanel.add(spikesFoundLabel);

            diagnosticPanel.add(new JLabel("Spike Rate (Hz/El)"));
            spikeRateLabel = new JLabel("0", JLabel.RIGHT);
            diagnosticPanel.add(spikeRateLabel);

            diagnosticPanel.add(new JLabel("TTL separation:"));
            refreshRateLabel = new JLabel("0", JLabel.RIGHT);
            diagnosticPanel.add(refreshRateLabel);

            if (saveRawData) {
                outputBuffer = new JProgressBar[dataSaver.getOutputStreamsCount()];
                for (int i = 0; i < outputBuffer.length; i++) {
                    diagnosticPanel.add(new JLabel("Output Buffer " + (i + 1)));
                    outputBuffer[i] = new JProgressBar(JProgressBar.HORIZONTAL, 0, 100);
                    diagnosticPanel.add(outputBuffer[i]);
                }
            }

            if (analysis /* && (analysisToDo == AnalysisToDo.SAVE_SPIKES ||
                 analysisToDo == AnalysisToDo.SAVE_SPIKES_AND_COVARINCES)*/) {
                diagnosticPanel.add(new JLabel("Spike Buffer"));
                spikeBufferLabel = new JLabel("0", JLabel.RIGHT);
                diagnosticPanel.add(spikeBufferLabel);
            }

            if (rawDataSource.startsWith("net://")) {
                System.out.println("Sending commence signal to DAQ computer");
                sampleInputStream.commenceWriting();
            }
        }
    }


    enum en {DO_NOTHING1, SAVE_SPIKES1};
}
