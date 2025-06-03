package edu.ucsc.neurobiology.vision.util;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

import edu.ucsc.neurobiology.vision.Vision;
import edu.ucsc.neurobiology.vision.anf.RawDataNoiseEvaluation;
import edu.ucsc.neurobiology.vision.anf.SpikeFinding;
import edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMap;
import edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory;
import edu.ucsc.neurobiology.vision.io.RawDataHeader;
import edu.ucsc.neurobiology.vision.io.SampleListener;
import edu.ucsc.neurobiology.vision.io.SpikeFile;
import edu.ucsc.neurobiology.vision.io.SpikeListener;
import edu.ucsc.neurobiology.vision.io.SpikeSaver;
import edu.ucsc.neurobiology.vision.math.MathUtil;

/**
 * 
 * 
 * @author Peter H. Li, The Salk Institute
 *
 */
public class StreamingSpikeFinder implements SampleListener {	
    private final RawDataHeader rawDataHeader;
    private final ElectrodeMap electrodeMap;
    private short[][] preBuffer;
    private int preBufferCount = 0;
    private int samplesCount = 0;
    private final int noiseEvalStart = 5;  // 5 seconds in
    private final int noiseEvalLength = 5; // 5 seconds long
    private final int noiseEvalStartSample;
    private final int noiseEvalEndSample;
    private final int samplesPerSecond = 20000;
    private NoiseEvaluationThread noiseEvalThread;
    private SpikeFinder spikeFinder = null;
    private boolean spikeFindingStarted = false;
    private final float spikeThreshold;
    private final float ttlThreshold;
    private final float timeConstant;
    private SpikeListener[] spikeListeners;
    private String outputFolder;
    
    public StreamingSpikeFinder(RawDataHeader rawDataHeader, float spikeThreshold, float ttlThreshold, 
            float timeConstant, String outputFolder, SpikeListener[] spikeListeners) {
        
        this.rawDataHeader = rawDataHeader;
        this.spikeThreshold = spikeThreshold;
        this.ttlThreshold = ttlThreshold;
        this.timeConstant = timeConstant;
        this.outputFolder = outputFolder;
        this.spikeListeners = spikeListeners;

        this.electrodeMap = ElectrodeMapFactory.getElectrodeMap(rawDataHeader.getArrayID());

        preBuffer = new short[noiseEvalLength * samplesPerSecond][];
        noiseEvalStartSample = noiseEvalStart * samplesPerSecond;
        noiseEvalEndSample = (noiseEvalStart + noiseEvalLength) * samplesPerSecond;
    }
    
    public StreamingSpikeFinder(RawDataHeader rawDataHeader, float spikeThreshold, float ttlThreshold, 
            float timeConstant, String outputFolder) {

        this(rawDataHeader, spikeThreshold, ttlThreshold, timeConstant, outputFolder, new SpikeListener[0]);
    }
    
    public void processSample(short[] sample) {
        if (spikeFindingStarted) {
            spikeFinder.processSample(sample);
            samplesCount++;
            return;
        }
        
        
        if (noiseEvalStartSample <= samplesCount && samplesCount < noiseEvalEndSample) {
            preBuffer[preBufferCount] = sample;
            preBufferCount++;
            
        } else if (samplesCount == noiseEvalEndSample) {
            noiseEvalThread = new NoiseEvaluationThread(preBuffer, electrodeMap.getNumberOfElectrodes(), outputFolder);
            noiseEvalThread.start();
            
        } else if (samplesCount > noiseEvalEndSample && noiseEvalThread.sigmasSet() && !spikeFindingStarted) {
//			float[] spikeThresholds = noiseEvalThread.scaleSpikeThresholds(spikeThreshold);
            float[] spikeThresholds = new float[electrodeMap.getNumberOfElectrodes()];
            Arrays.fill(spikeThresholds, spikeThreshold);
            spikeThresholds[0] = 100;
            spikeFinder = new SpikeFinder(electrodeMap, spikeThresholds, ttlThreshold, timeConstant);
            spikeFinder.initialize();
            
            for (SpikeListener sl : spikeListeners) spikeFinder.addSpikeListener(sl);

            SpikeBuffer spikeBuffer = new SpikeBuffer(SpikeFinding.SPIKE_BUFFER_SIZE);
            spikeFinder.addSpikeListener(spikeBuffer);

            String setName = new File(outputFolder).getName();
            String spikesFileName = outputFolder + File.separator + setName + VisionParams.SPIKES_FILE_EXTENSION;
            try {
                SpikeFile spikeFile = new SpikeFile(spikesFileName, rawDataHeader.getArrayID(),
                        timeConstant, spikeThreshold,
                        rawDataHeader.getNumberOfSamples(),
                        rawDataHeader.getSamplingFrequency());
                SpikeSaver spikeSaver = new SpikeSaver(spikeFile);
                spikeBuffer.addSpikeListener(spikeSaver);
            } catch (IOException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
            
            spikeFindingStarted = true;
        } 
        
        samplesCount++;
    }

    public void finishSampleProcessing() throws IOException {
        if (spikeFinder != null) {
            System.out.println("StreamingSpikeFinder: Found " + spikeFinder.numSpikesFound() + " spikes in " + spikeFinder.numSamplesProcessed() + " samples processed out of " + samplesCount + " samples total");
            spikeFinder.finishSampleProcessing();
        }
    }
    
    
    /**
     * 
     * @author Peter H. Li, The Salk Institute
     *
     */
    static class NoiseEvaluationThread extends Thread {
        private Vision app = Vision.getInstance();
        private final int nElectrodes;
        private final int nSamples;
        private final short[][] samples;
        private double[] sigmas = null;
        private final String outputFolder;
        
        public NoiseEvaluationThread(short[][] samples, int nElectrodes, String outputFolder) {
            this.nElectrodes = nElectrodes;
            this.nSamples = samples.length;
            this.outputFolder = outputFolder;
            
            this.samples = samples;
//			// Deep copy samples
//			this.samples = new short[nSamples][];
//			for (int i = 0; i < nSamples; i++) {
//				this.samples[i] = samples[i].clone();
//			}
        }
        
        public void run() {
            // These blocks copied and pasted from RawDataNoiseEvaluation; should be abstracted!
            System.out.println("StreamingSpikeFinder: Calculating rms noise with " + nSamples + " samples over " + nElectrodes + " electrodes...");
            app.startProgressBar();

            double[] sigmas = new double[nElectrodes];
            for (int e = 1; e < nElectrodes; e++) {
                short[] data = new short[nSamples];
                for (int i = 0; i < nSamples; i++) {
                    data[i] = samples[i][e];
                }
                sigmas[e] = RawDataNoiseEvaluation.getNoise(data, nSamples);

                app.setProgress(100 * e / (nElectrodes - 1));
            }
            app.endProgressBar();
            
            String setName = new File(outputFolder).getName();
            String sigmasFileName = outputFolder + File.separator + setName + VisionParams.NOISE_FILE_EXTENSION;
            System.out.println("StreamingSpikeFinder: Saving noise to " + sigmasFileName);

            // This allows spike finding to begin
            synchronized (this) {
                this.sigmas = sigmas;
            }
            
            try {
                File f = new File(sigmasFileName);
                if (f.exists())	f.delete();

                PrintWriter pw = new PrintWriter(new FileWriter(sigmasFileName));
                for (int i = 0; i < nElectrodes; i++) {
                    pw.println(sigmas[i]);
                }
                pw.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        
        public synchronized boolean sigmasSet() {
            return sigmas != null;
        }

        /**
         * If called before thread is run, returns null.
         * @return
         */
        public float[] scaleSpikeThresholds(float spikeThreshold) {
            float[] spikeThresholds = new float[sigmas.length];
            for (int i = 0; i < sigmas.length; i++) {
                spikeThresholds[i] = (float) sigmas[i];
            }
            MathUtil.multiply(spikeThresholds, spikeThreshold);
            return spikeThresholds;
        }
    }
}