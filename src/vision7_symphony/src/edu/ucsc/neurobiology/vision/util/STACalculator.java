package edu.ucsc.neurobiology.vision.util;

import java.io.*;

import edu.ucsc.neurobiology.vision.stimulus.*;


/**
 * A class that implements STA calculation by averaging of movie frames at given
 * spike times.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class STACalculator {
    
    private UpdatingFrame[] staFrames;
    private BufferedMovie bufferedMovie;
    private final double refreshTime, firstImageChangeTime;
    private final double conversionFactor;
    private final int staDepth,staOffset, movieSize;
    
    
    /**
     * The movie has to be reset in advance.
     *
     * @param staDepth int
     * @param firstTTL int
     * @param timePerSample double
     * @param bufferedMovie BufferedMovie
     */
    public STACalculator(
        int staDepth, int staOffset, int firstTTL, double timePerSample, BufferedMovie bufferedMovie) {

        this.staDepth = staDepth;
        this.staOffset = staOffset;
        this.refreshTime = bufferedMovie.getRefreshTime();
        this.conversionFactor = 1 / (refreshTime / timePerSample);
        this.firstImageChangeTime = firstTTL  + refreshTime / timePerSample;
        this.bufferedMovie = bufferedMovie;
        this.movieSize = bufferedMovie.size();

        staFrames = new UpdatingFrame[staDepth];
        for (int i = 0; i < staDepth; i++) {
            if (bufferedMovie.getUnderlyingMovie() instanceof WhiteNoiseMovie) {
                staFrames[i] = ( (WhiteNoiseMovie) bufferedMovie.getUnderlyingMovie()).
                               createCompatibleFrame();
            } else {
                staFrames[i] = ( (RawMovie) bufferedMovie.getUnderlyingMovie()).
                               createCompatibleFrame();
            }
        }
    }


    /**
     * This is now essentially abstracted into {@link #spikeTimesToFirstFrameIndices(int[])}
     * and {@link #addFramesAt(int)}, but leave it together here as it is performance
     * critical.
     *
     * @param time
     * @throws IOException
     */ 
    public void addSpike(int time) throws IOException {
        final int firstFrame = (int) ( (time - firstImageChangeTime) * conversionFactor) - staDepth + 2 + staOffset;
        if (firstFrame < 0 || firstFrame + staDepth >= movieSize) return;
        for (int f = 0; f < staDepth; f++)
            staFrames[f].update(bufferedMovie.getEncodedFrame(firstFrame + f));
    }
    
    public void addFramesAt(int firstFrameIndex) throws IOException {
        for (int f = 0; f < staDepth; f++) {
            staFrames[f].update(bufferedMovie.getEncodedFrame(firstFrameIndex + f));
        }
    }


    public int[] spikeTimesToFirstFrameIndices(int[] spikeTimes) {
        int[] firstFrameIndices = new int[spikeTimes.length];
        for (int i = 0; i < spikeTimes.length; i++) {
            firstFrameIndices[i] = (int) ((spikeTimes[i] - firstImageChangeTime) * conversionFactor) - staDepth + 2 + staOffset;
        }
        return firstFrameIndices;
    }
    


    public void finish() {
        for (int f = 0; f < staDepth; f++) {
            staFrames[f].finish();
        }
    }


    public STA getSTA() throws IOException {
        STA sta = new STA(
            staDepth,
            bufferedMovie.getWidth(),
            bufferedMovie.getHeight(),
            refreshTime,
            bufferedMovie.getFrame(0).getStixelWidth(),
            bufferedMovie.getFrame(0).getStixelHeight()
                  );

        for (int f = 0; f < staDepth; f++) {
            ( (STAFrame) sta.getFrame(f)).set(staFrames[f].getColorBuffer(),
                                              staFrames[f].getErrorBuffer());
        }

        return sta;
    }


    public STA getSTV() throws IOException {
        STA ste = new STA(
            staDepth,
            bufferedMovie.getWidth(),
            bufferedMovie.getHeight(),
            refreshTime,
            bufferedMovie.getFrame(0).getStixelWidth(),
            bufferedMovie.getFrame(0).getStixelHeight()
                  );

        for (int f = 0; f < staDepth; f++) {
            ( (STAFrame) ste.getFrame(f)).set(staFrames[f].getCovarianceBuffer(),
                                              staFrames[f].getErrorBuffer());
            //.getErrorBuffer() gets the wrong errors.  I am leaving it this way, because
            //I am not using the error for the spike triggered variance
        }

        return ste;
    }


    public static STA calculate(int[] time, int staDepth, int staOffset, int firstTTL,
                                double timePerSample, BufferedMovie bufferedMovie) throws
        IOException {

        STACalculator c = new STACalculator(
            staDepth, staOffset, firstTTL, timePerSample, bufferedMovie);

        int percent = -1, oldPercent = -1;

        for (int i = 0; i < time.length; i++) {
            c.addSpike(time[i]);

            percent = (int) (i / (double) time.length);
            if (percent != oldPercent) {
                System.out.println(percent + " %");
                oldPercent = percent;
            }
        }

        c.finish();
        return c.getSTA();
    }
}
