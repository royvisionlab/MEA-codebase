package edu.ucsc.neurobiology.vision.stimulus;

import java.io.*;

import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.io.chunk.*;


/**
 * This class defines methods to create a new raw movie or open an old raw movie.
 *
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class RawMovie
implements AdvancedMovie {

    private final int nFrames;
    private final int width;
    private final int height;

    private final double stixelWidth, stixelHeight;
    private double refreshTime;
    private int[] iDroppedFrames;
    public FrameGenerator frameGenerator = null;

    private STAFrame generalFrame;
    private FrameEncoding frameEncoding = FrameEncoding.PACKED_ENCODING;
    RawMovieFile rawFile = null;
    private boolean calculateSTV = true;

    public RawMovie(String fileName, String globalsFileName) throws IOException {

        boolean rtParamsOpen = true;
        rawFile = new RawMovieFile(fileName);

        int nFramesGenerated = rawFile.getFramesGenerated();
        width = rawFile.getWidth();
        height = rawFile.getHeight();

        GlobalsFile.RunTimeMovieParams runTime = null;
        try {
            GlobalsFile globalsFile = new GlobalsFile(globalsFileName, ChunkFile.READ);
            runTime = globalsFile.getRunTimeMovieParams();

        } catch (IOException e) {
            System.out.println(
            "Globals file not found.  Using default values.  If you are calculating STAs, abort and try again.");
            rtParamsOpen = false;
        }

        if (rtParamsOpen) {
            nFrames = runTime.nFramesRequired;
            stixelWidth = runTime.micronsPerStixelX;
            stixelHeight = runTime.micronsPerStixelY ;
            refreshTime = runTime.refreshPeriod;
            iDroppedFrames = runTime.droppedFrames;
        } else {
            nFrames = nFramesGenerated;
            stixelWidth = 5.8;
            stixelHeight = 5.8;
            refreshTime = -1;
            this.iDroppedFrames = new int[0];
        }


        generalFrame = new STAFrame(
                width, height, stixelWidth, stixelHeight, new float[width * height * 3]);
    }





    public void setCalculateSTV(boolean calculateSTV) {
        this.calculateSTV = calculateSTV;
    }


    public void setFrameEncoding(FrameEncoding frameEncoding) {
        this.frameEncoding = frameEncoding;
    }


    public Object getEncodedFrame(int frameIndex, Object frame) throws IOException {
        int actualFrame = frameIndex;
        for (int i = 0; i < iDroppedFrames.length; i++) {
            if (frameIndex >= iDroppedFrames[i]) {
                actualFrame--;
            }
        }
        switch (frameEncoding) {
        case INT_ENCODING:
            frame = rawFile.readFrameInt(actualFrame);
            break;

        case PACKED_ENCODING:
            frame = rawFile.readFrameLong(actualFrame, calculateSTV);
            break;

        case FLOAT_ARRAY_ENCODING:

//			frame = rawFile.readFrameFloat(actualFrame);
            throw new Error("FLOAT_ARRAY_ENCODING");
//			break;
        }

        return frame;
    }


    public UpdatingFrame createCompatibleFrame() {
        return new GaussUpdatingFrame(width, height, calculateSTV);
    }


    public ImageFrame getFrame(int frameIndex) throws IOException {
        int actualFrame = frameIndex;
        for (int i = 0; i < iDroppedFrames.length; i++) {
            if (frameIndex > iDroppedFrames[i]) {
                actualFrame--;
            }
        }

        generalFrame.setBuffer(rawFile.readFrameFloat(actualFrame));
        return generalFrame;
    }


    public String getDescription() {
        return "Raw Movie(" + size() + ")";
    }


    public double getRefreshTime() {
        return refreshTime;
    }


    public int getWidth() {
        return width;
    }


    public int getHeight() {
        return height;
    }


    @Override
    public int size() {
        return nFrames;
    }
}
