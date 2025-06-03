package edu.ucsc.neurobiology.vision.matlab;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import edu.ucsc.neurobiology.vision.Config;
import edu.ucsc.neurobiology.vision.Vision;
import edu.ucsc.neurobiology.vision.io.RawDataFile;
import edu.ucsc.neurobiology.vision.io.RawDataSaver;
import edu.ucsc.neurobiology.vision.stimulus.BinaryFrameGenerator;
import edu.ucsc.neurobiology.vision.stimulus.BufferedMovie;
import edu.ucsc.neurobiology.vision.stimulus.FrameGenerator;
import edu.ucsc.neurobiology.vision.stimulus.FrameGenerator.ColorType;
import edu.ucsc.neurobiology.vision.stimulus.FrameGenerator.RandomNumberGenerator;
import edu.ucsc.neurobiology.vision.stimulus.MovieType;
import edu.ucsc.neurobiology.vision.stimulus.STA;
import edu.ucsc.neurobiology.vision.stimulus.WhiteNoiseMovie;
import edu.ucsc.neurobiology.vision.util.AsyncAccumulator;
import edu.ucsc.neurobiology.vision.util.STACalculator;
import edu.ucsc.neurobiology.vision.util.VisionParams;
import edu.ucsc.neurobiology.vision.util.WaveformCalculator;

/**
 * These methods should be the standard means of generating STAs, white noise movies, and EIs from MATLAB.
 * New static wrapper methods that interface with the core vision code and are intended to interface with MATLAB
 * should be added to this class. More complicated wrapper classes (like the Read* classes) should be added to the
 * MATLAB package instead of elsewhere. Cell-Finder should never be called from MATLAB.
 * 
 * @author tamachado@salk.edu (6/15/08)
 * 
 */
public class Matlab {
    
    public static AsyncAccumulator<float[]> getFloatFrameAccumulator(final FrameGenerator fg, int queueLength) {
        return AsyncAccumulator.build(new AsyncAccumulator.Accumulator<float[]>() {
            public float[] accumulate() throws Throwable {
                return fg.nextFloatFrame();
            }
        }, queueLength);
    }

    /**
     * 
     * (Currently only BinaryFrameGenerator implements nextMGLFrame, but this will be expanded.)
     * 
     * @param fg
     * @param queueLength
     * @return
     */
    public static AsyncAccumulator<char[]> getMGLFrameAccumulator(final BinaryFrameGenerator fg, int queueLength) {
        return AsyncAccumulator.build(new AsyncAccumulator.Accumulator<char[]>() {
            public char[] accumulate() throws Throwable {
                return fg.nextMGLFrame();
            }
        }, queueLength);
    }


    public static BufferedMovie computeMovie(String rawMovieFileName, String globalsFileName) throws IOException {
        boolean generateSeeds = true;
        WhiteNoiseMovie wnm = new WhiteNoiseMovie(rawMovieFileName, globalsFileName, generateSeeds);
        return new BufferedMovie(wnm);
    }
    
    /**
     * DEPRECATED: the new version above is much cleaner and simpler and stays in sync with changes in the movie-xml
     * file format.
     * Compute a white noise movie
     */
    public static BufferedMovie computeMovie(String movieXmlFileName, int[] ttlTimes) throws IOException {

        // Parse the movie xml file
        Config movieXml = null;
        try{ 
            movieXml = new Config(movieXmlFileName);
        } catch (Exception e) {
            Vision.reportFatalException("Error reading movie XML file!", e);
        }

        // Put the parameters into a hash map
        HashMap<String, String> wnParams = movieXml.getParameterList("Make White Noise Movie");
        HashMap<String, String> auxParams = movieXml.getParameterList("Calculate Auxiliary Parameters");

        // Optional params to make sparse
        boolean sparse = false;
        float probability = 0;
        if (wnParams.containsKey("Sparse")) {
            sparse = Boolean.parseBoolean(wnParams.get("Sparse"));
            if (sparse) probability = Float.parseFloat(wnParams.get("Sparse.Probability"));
        }
        
        // Load params for determining dropped frames, if not specified, use defaults from VisionParams
        int refreshInterval = Integer.parseInt(auxParams.get("Set Movie.refreshInterval"));
        double monitorFrequency = VisionParams.DEFAULT_MONITOR_FREQUENCY;
        if (auxParams.containsKey("Set Movie.Monitor Frequency")) 
            monitorFrequency = Double.valueOf(auxParams.get("Set Movie.Monitor Frequency"));

        int framesPerTTL = VisionParams.DEFAULT_FRAMES_PER_TTL;
        if (auxParams.containsKey("Set Movie.Frames Per TTL"))
            framesPerTTL = Integer.valueOf(auxParams.get("Set Movie.Frames Per TTL"));		
        
        // Determine dropped frames
        ArrayList<Integer> droppedFramesList = new ArrayList<Integer> ();
        double refreshPeriod = WhiteNoiseMovie.calculateRefreshPeriod(droppedFramesList,
                ttlTimes, refreshInterval, monitorFrequency, framesPerTTL);
        
        // the frame number of each dropped frame
        int[] droppedFrameNumbers = new int[droppedFramesList.size()];
        for (int i = 0; i < droppedFramesList.size(); i++) {
            droppedFrameNumbers[i] = (int) (droppedFramesList.get(i));
        }

        // figure out nFrames required (from AuxillaryParametersCalculator)
        int nFrames = ( ttlTimes.length * VisionParams.DEFAULT_FRAMES_PER_TTL) / refreshInterval + 1;

        // make the white noise movie using parameters from the xml file
        boolean generateSeeds = true;
        long[] seeds = null;
        WhiteNoiseMovie movie = new WhiteNoiseMovie(
                Integer.parseInt(wnParams.get("Width")), 
                Integer.parseInt(wnParams.get("Height")), 
                Integer.parseInt(auxParams.get("Set Movie.pixelsPerStixelX")), 
                Integer.parseInt(auxParams.get("Set Movie.pixelsPerStixelY")), 
                refreshPeriod,
                nFrames, 
                Long.parseLong(wnParams.get("Seed")), 
                ColorType.values()[(int)Float.parseFloat(wnParams.get("ColorType"))], 
                MovieType.values()[(int)Float.parseFloat(wnParams.get("NoiseType"))],
                sparse, probability,
                RandomNumberGenerator.values()[(int)Float.parseFloat(wnParams.get("RandomNumberGenerator"))], 
                Float.parseFloat(wnParams.get("ContrastSigma")), 
                droppedFrameNumbers,
                seeds, generateSeeds);

        // return buffered version of movie
        return new BufferedMovie(movie);
    }
    
    
    // Compute a single STA
    // staOffset of zero means that spike occurs after last frame.
    // staOffset of -1 means that spike occurs in last frame.
    public static STA computeSTA(int[] spikeTimes, int firstTTL,
                                   BufferedMovie bufferedMovie, int staDepth, int staOffset)
            throws IOException {
        // reset buffer
        bufferedMovie.reset();

        // calculate times per sample (from STACalculation.startCalculation)
        double timePerSample = 1 / VisionParams.SAMPLES_PER_MILLISECOND;

        // create STACalculator
        STACalculator staCalculator = new STACalculator(staDepth, staOffset,
                firstTTL, timePerSample, bufferedMovie);

        // add spikes sequentially
        for (int i = 0; i < spikeTimes.length; i++)
            staCalculator.addSpike(spikeTimes[i]);

        // finish up
        staCalculator.finish();
        return staCalculator.getSTA();
    }


    // RawDataSaver must be created from Matlab
    public static void saveRawData(RawDataSaver saver, short[][] samples) {
        for (int i = 0; i < samples.length; i++)
            saver.processSample(samples[i]);
    }
    
    
    // 2010-02 phli
    // Wrapper for below; this version accepts MatLab cell array of int32 arrays to 
    // match the ragged array of int[] expected by Java.
    public static STA[] computeSTAs(Object[] spikeCA,
            int firstTTL, BufferedMovie bufferedMovie, int staDepth) 
            throws IOException {
        
        int[][] spikeTimes = new int[spikeCA.length][];
        for (int i = 0; i < spikeCA.length; i++) {
            spikeTimes[i] = (int[]) spikeCA[i];
        }
        
        return computeSTAs(spikeTimes, firstTTL, bufferedMovie, staDepth);
    }
    
    
    // 2010-02 phli
    // Compute multiple STAs.  Try to be efficient by only looping through bufferedMovie once.  
    // This can probably be improved; I simply ported it over from cellfinder.ReceptiveField.calculateManySTAs
    public static STA[] computeSTAs(int[][] spikeTimes,
            int firstTTL, BufferedMovie bufferedMovie, int staDepth, int staOffset)
            throws IOException {

        // determine number of neurons
        int neurons = spikeTimes.length;

        // calculate times per sample (from STACalculation.startCalculation)
        double timePerSample = 1 / VisionParams.SAMPLES_PER_MILLISECOND;

        // generate STA calculators for each cluster
        STACalculator[] staCalculators = new STACalculator[neurons];
        for (int i = 0; i < neurons; i++)
            staCalculators[i] = new STACalculator(staDepth, staOffset, firstTTL, timePerSample, bufferedMovie);

        // number of spikes processed
        int[] processed = new int[neurons];
        for (int i = 0; i < neurons; i++) processed[i] = 0;

        // loop through ALL cluster spikes sequentially (thus, interleaved)
        // the idea is to make one pass through the BufferedMovie
        // Code derived from STACalculation.STAThread.run()
        STA[] stas = new STA[neurons];
        while (true) {

            // determine first neuron with spikes remaining
            int min;
            for (min = 0; min < neurons; min++) {
                if (processed[min] < spikeTimes[min].length) break;
            }

            // if no neurons remain, then finished
            if (min == neurons) break;

            // starting with "min" neuron, find neuron with earliest spike
            for (int i = 0; i < neurons; i++) {
                //if neuron has spikes to process
                if (processed[i] < spikeTimes[i].length) {
                    // if neuron has earlier spike then min
                    if (spikeTimes[i][processed[i]] < spikeTimes[min][processed[min]]) {
                        min = i;
                    }
                }
            }

            // add the spike!
            staCalculators[min].addSpike(spikeTimes[min][processed[min]]);
            processed[min]++;
        }

        // finish it up and save
        for (int i = 0; i < neurons; i++) {
            staCalculators[i].finish();
            stas[i] = staCalculators[i].getSTA();
        }

        return stas;
    }
    
    
    // Compute a single electrophysiological image
    public static float[][][] computeElectrophysiologicalImage(String rawDataFile, int[] spikeTimes, int nlPoints, int nrPoints, int nSpikesToAverage) throws IOException {
        // load up raw data
        RawDataFile rawData = new RawDataFile(new File(rawDataFile));
        int nElectrodes = rawData.getHeader().getNumberOfElectrodes();
        int nSamples = rawData.getHeader().getNumberOfSamples();
        int nSpikes = Math.min(nSpikesToAverage, spikeTimes.length);

        // compiling EI statistics for each spike
        WaveformCalculator average = new WaveformCalculator(nElectrodes, nlPoints, nrPoints);
        short[][] samples = average.createCompatibleBuffer();

        // loop through each spike
        for (int i = 0; i < nSpikes; i++) {
            int time = spikeTimes[i];

            final int startSample = time - nlPoints - 1;
            if (startSample > 0 && startSample < nSamples - (samples.length)) {
                rawData.getData(startSample, samples);
                average.addSpike(samples);
            }
        }

        // close the raw data files
        rawData.close();

        // calculate EI
        average.finish();
        final float[][] image = average.getAverage();
        final float[][] error = average.getError();

        // remove the mean
        final int nAverage = image[0].length;  // could be shortened to 10

        for (int electrode = 0; electrode < nElectrodes; electrode++) {

            // calculate mean for particular electrode across time
            double mean = 0;
            for (int i = 0; i < nAverage; i++) {
                mean += image[electrode][i];
            }
            mean /= nAverage;

            // subtract off mean at every moment of time
            for (int i = 0; i < image[electrode].length; i++) {
                image[electrode][i] -= mean;
            }
        }

        return new float[][][] {image, error};
    }
    
    // Calculate the average frame refresh rate for an acquired data set
    // This calculation makes many assumptions and is largely based off of
    // earlier code by M. Grivich with further consultation with EJC.
    // Edited by shlens.
    public static double computeRigFrameTime(final int[] ttlTimes) {
        
        // Initialize the tallying
        double sum = 0.0;
        int nTTL = 0;

        // calculate nominal monitor refresh duration
        double refreshDuration = 1 / VisionParams.DEFAULT_MONITOR_FREQUENCY;


        // loop through the TTL pulses and tally "good" pulse durations
        for (int i = 1; i < ttlTimes.length; i++) {

            // calculate duration of each TTL pulse
            double durationTTL = (double) (ttlTimes[i] - ttlTimes[i - 1]);

            // calculate the number of frames per TTL pulse duration
            long nFramesPerTTL = Math.round( (durationTTL / (VisionParams.SAMPLES_PER_MILLISECOND * 1000))
                    / refreshDuration);

            // only add to tally if no frames were dropped
            if (nFramesPerTTL == VisionParams.DEFAULT_FRAMES_PER_TTL) {
                sum += durationTTL;
                nTTL++;
            }

        }

        // average duration of TTL pulse
        double avgDurationTTL = sum / nTTL;

        // rig frame time (units = samples)
        double rigFrameTime = avgDurationTTL / VisionParams.DEFAULT_FRAMES_PER_TTL;

        // convert: recording samples --> msec
        rigFrameTime = rigFrameTime / VisionParams.SAMPLES_PER_MILLISECOND;

        return rigFrameTime;

    }
}
