package edu.ucsc.neurobiology.vision.stimulus;

import java.io.*;
import java.util.*;


/**
 * A movie that buffers the frames produced by another movie so that there is only ONE
 * creation of an underlying movie frame.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class BufferedMovie implements AdvancedMovie {

    private final int bufferSize;
    private AdvancedMovie movie;
    private LinkedList<Object> frameBuffer;
    private int bufferEnd, bufferStart;


    /**
     * Creates a BufferedMovie with a buffer size of 1000 frames.
     *
     * @param movie
     */
    public BufferedMovie(AdvancedMovie movie) throws IOException {
        this(movie, 1000);
    }


    public BufferedMovie(AdvancedMovie movie, int bufferSize) throws IOException {
        this.movie = movie;
        bufferSize = Math.min(bufferSize, movie.size());
        this.bufferSize = bufferSize;

        frameBuffer = new LinkedList<Object>();
        reset();
    }

    public void reset() throws IOException {
        jumpToFrame(0);
    }
    
    public void setBufferStartFrame(int frameIndex) throws IOException {
        if (frameIndex >= bufferStart && frameIndex < bufferEnd) {
            rollToFrame(frameIndex);
        } else {
            jumpToFrame(frameIndex);
        }
    }
    
    public void setBufferEndFrame(int frameIndex) throws IOException {
        setBufferStartFrame(frameIndex - bufferSize + 1);
    }
        
    private void rollToFrame(int frameIndex) throws IOException {
        int framesToAdvance = frameIndex - bufferStart;
        for (int i = 1; i <= framesToAdvance; i++) {
            frameBuffer.removeFirst();
            frameBuffer.addLast(movie.getEncodedFrame(bufferEnd + i, null));
        }
        bufferEnd   += framesToAdvance;
        bufferStart += framesToAdvance;
    }
    
    private void jumpToFrame(int frameIndex) throws IOException {
        frameBuffer.clear();
        setBufferStart(frameIndex);
        for (int i = bufferStart; i <= bufferEnd; i++) {
            frameBuffer.addLast(movie.getEncodedFrame(i, null));
        }
    }
    
    private void setBufferStart(int frameIndex) {
        bufferStart = frameIndex;
        bufferEnd = bufferStart + bufferSize - 1;
    }
    
    
//    /**
//     * New version, cleaner but not significantly faster, even when logic is moved inline.
//     * 
//     * @param frameIndex
//     * @return
//     * @throws IOException
//     */
//    public Object getEncodedFrame(int frameIndex) throws IOException {
//        if (frameIndex >= bufferStart && frameIndex <= bufferEnd) {
//        	return frameBuffer.get(frameIndex - bufferStart);
//        } else {
//    		  setBufferEndFrame(frameIndex);
//            return frameBuffer.getLast();
//        }
//    }

    // Old version
    public Object getEncodedFrame(int frameIndex) throws IOException {
        if (frameIndex <= bufferEnd) {
            if (frameIndex >= bufferStart) {
                return frameBuffer.get(frameIndex - bufferStart);
            } else {
                throw new RuntimeException("Too old in the past");
            }
        } else {
            int missingFrames = frameIndex - bufferEnd;
            for (int i = 1; i <= missingFrames; i++) {
                frameBuffer.removeFirst();
                //Matthew changed this line.  In create pink movie the returned frame is valid,
                //while the passed reference is not.
                //movie.getEncodedFrame(bufferEnd + i, frame);
                Object frame = movie.getEncodedFrame(bufferEnd + i, null);
                frameBuffer.addLast(frame);
            }
            bufferEnd += missingFrames;
            bufferStart += missingFrames;

            return frameBuffer.getLast();
        }
    }
    
    
    public Object getEncodedFrame(int frameIndex, Object frame) {
        throw new IllegalAccessError("Method not implemented.");
    }


    public ImageFrame getFrame(int frameIndex) throws IOException {
        return movie.getFrame(frameIndex);
    }


    public String getDescription() {
        return "Buffered Movie";
    }


    public double getRefreshTime() {
        return movie.getRefreshTime();
    }


    public int getWidth() {
        return movie.getWidth();
    }


    public int getHeight() {
        return movie.getHeight();
    }


    public int size() {
        return movie.size();
    }


    public Movie getUnderlyingMovie() {
        return movie;
    }


    public void setFrameEncoding(FrameEncoding frameEncoding) {
        movie.setFrameEncoding(frameEncoding);
    }

}
