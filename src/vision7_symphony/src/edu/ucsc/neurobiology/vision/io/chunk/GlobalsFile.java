package edu.ucsc.neurobiology.vision.io.chunk;

import java.io.*;
import edu.ucsc.neurobiology.vision.io.RawDataHeader512;
import edu.ucsc.neurobiology.vision.util.*;

public class GlobalsFile extends ChunkFile {
    static final boolean DEBUG = false;
    
    static final long fileID = 237L + ( 202L << 8L) + (168L << 16L)	+ (33L << 24L) + (26L << 32L) + (125L << 40L) + (6L << 48L) + (101L << 56L);
    
    // Choices for created by
    public static final int DEFAULT = 0;
    public static final int UPDATER = 1;
    
    public static final int ICP_TAG = 0;
    public static final int RTMP_TAG = 1;
    private static final int CREATED_BY_TAG = 2;
    public static final int VERSION_TAG = 3;
    private static final int MONITOR_FREQ_TAG = 4; // Has own TAG, but XML and Java fields actually live inside RTMP
    private static final int RDH512_TAG = 5;

    public GlobalsFile(String fileName) throws IOException {
        this(fileName, new File(fileName).exists() ? ChunkFile.READ : ChunkFile.WRITE);
    }
    
    public GlobalsFile(String fileName, int mode) throws IOException {
        super(fileName, mode);
        if (mode == ChunkFile.WRITE) writeVisionVersion(VisionParams.VERSION);
    }
    
    public long getFileID() { return fileID; }
    
    public String toString() {
        String nl = System.getProperty("line.separator");
        String toReturn = "Globals File:" + nl + nl;
        try {
            if (!imageCalibrationParamsExists())
                toReturn = "Image Calibration Parameters do not exist." + nl +nl;
            else toReturn = new ImageCalibrationParams() + nl;
            
            if (!runTimeMovieParamsExists())
                toReturn += "Run Time Image Calibration Parameters do not exist." + nl + nl;
            else toReturn += new RunTimeMovieParams() + nl;
            
            if (!rawDataHeader512Exists())
                toReturn += "Raw Data Header information does not exist." + nl + nl;
            else toReturn += "Raw Data Header:" + readRDH512() + nl + nl;
            
            toReturn += "createdBy: " + getCreatedBy() +nl;
            toReturn += "Vision Version: " + VisionParams.versionString(getVisionVersion());
        } catch (IOException e) { e.printStackTrace(); }

        return toReturn;
    }
    
    public boolean runTimeMovieParamsExists()     throws IOException { return getChunk(RTMP_TAG)   != null; }
    public boolean imageCalibrationParamsExists() throws IOException { return getChunk(ICP_TAG)    != null; }
    public boolean rawDataHeader512Exists()       throws IOException { return getChunk(RDH512_TAG) != null;	}
    
    public RawDataHeader512 readRDH512() throws IOException {
        byte[] buffer = getChunk(RDH512_TAG);
        DataInputStream dis = new DataInputStream(new ByteArrayInputStream(buffer));
        RawDataHeader512 rdh512 = new RawDataHeader512(dis);
        dis.close();
        
        if (DEBUG) System.out.println(rdh512);
        return rdh512;
    }
    
    public void writeRDH512(RawDataHeader512 rdh512) throws IOException {
        addChunk(RDH512_TAG, rdh512.getBinaryRepresentation());
        if (DEBUG) System.out.println(toString());		
    }
    
    
    public ImageCalibrationParams getImageCalibrationParams() throws IOException {
        return new ImageCalibrationParams();
    }
    
    public void setImageCalibrationParams(ImageCalibrationParams calParams) throws IOException {
        new ImageCalibrationParams(calParams.micronsPerPixelX, calParams.micronsPerPixelY, calParams.centerX,
                calParams.centerY, calParams.flipX, calParams.flipY, calParams.angle, calParams.arrayID, 
                calParams.arrayPart, calParams.arrayNParts);
    }
    
    public void setImageCalibrationParams(double micronsPerPixelX, double micronsPerPixelY, 
            double centerX, double centerY, boolean flipX, boolean flipY, 
            double angle, int arrayID, int arrayPart, int arrayNParts) throws IOException {
        new ImageCalibrationParams(micronsPerPixelX, micronsPerPixelY, 
                centerX, centerY, flipX, flipY, angle, arrayID, arrayPart, arrayNParts);
    }
        
    public class ImageCalibrationParams {
        public double micronsPerPixelX, micronsPerPixelY;
        public double centerX, centerY; //microns
        public boolean flipX, flipY; 
        public double angle;//radians
        public int arrayID, arrayPart, arrayNParts;		
        
        //For writing
        ImageCalibrationParams(double micronsPerPixelX, double micronsPerPixelY, 
                double centerX, double centerY, boolean flipX, boolean flipY, 
                double angle, int arrayID, int arrayPart, int arrayNParts) throws IOException {
            
            ByteArrayOutputStream bos = new ByteArrayOutputStream(0);
            DataOutputStream dos = new DataOutputStream(bos);
            dos.writeDouble(micronsPerPixelX);
            dos.writeDouble(micronsPerPixelY);
            dos.writeDouble(centerX);
            dos.writeDouble(centerY);
            dos.writeBoolean(flipX);
            dos.writeBoolean(flipY);
            dos.writeDouble(angle);
            dos.writeInt(arrayID);
            dos.writeInt(arrayPart);
            dos.writeInt(arrayNParts);
            addChunk(ICP_TAG, bos.toByteArray());	
            
            if (DEBUG) System.out.println(toString());
        }
        
        //For reading
        ImageCalibrationParams() throws IOException {
            byte[] buffer = getChunk(ICP_TAG);
            
            DataInputStream dis = new DataInputStream(new ByteArrayInputStream(buffer));
            micronsPerPixelX = dis.readDouble();
            micronsPerPixelY = dis.readDouble();
            centerX = dis.readDouble();
            centerY = dis.readDouble();
            flipX = dis.readBoolean();
            flipY = dis.readBoolean();
            angle = dis.readDouble();
            arrayID = dis.readInt();
            arrayPart = dis.readInt();
            arrayNParts = dis.readInt();
            dis.close();
            
            if (DEBUG) System.out.println(toString());
        }
        
        /**
         * ImageCalibrationParams
         */
        public String toString() {
            String newline = System.getProperty("line.separator");
            String toReturn = 
                "micronsPerPixelX: " + micronsPerPixelX + newline +
                "micronsPerPixelY: "  + micronsPerPixelY + newline +
                "centerX: " + centerX + newline + 
                "centerY: " + centerY + newline +
                "flipX: " + flipX+ newline + 
                "flipY: " + flipY + newline +
                "angle:" + angle + newline + 
                "arrayID: " + arrayID + newline + 
                "arrayPart: " + arrayPart + newline + 
                "arrayNParts: " + arrayNParts + newline;
            
            return toReturn;
        }
    }
    
    
    public RunTimeMovieParams getRunTimeMovieParams() throws IOException {
        return new RunTimeMovieParams();
    }
    
    public void setRunTimeMovieParams(RunTimeMovieParams rtParams) throws IOException {
        new RunTimeMovieParams(rtParams.pixelsPerStixelX, rtParams.pixelsPerStixelY, rtParams.width,
                rtParams.height, rtParams.micronsPerStixelX, rtParams.micronsPerStixelY, 
                rtParams.xOffset, rtParams.yOffset, rtParams.interval, rtParams.monitorFrequency, rtParams.framesPerTTL, 
                rtParams.refreshPeriod,	rtParams.nFramesRequired, rtParams.droppedFrames);
    }
    
    public void setRunTimeMovieParams(int pixelsPerStixelX, int pixelsPerStixelY, double width, 
            double height, double micronsPerStixelX, double micronsPerStixelY,
            double xOffset, double yOffset, int interval, double monitorFrequency, int framesPerTTL, 
            double refreshPeriod, int nFramesRequired, int[] droppedFrames) throws IOException {
        new RunTimeMovieParams(pixelsPerStixelX,  pixelsPerStixelY, width, 
                height,  micronsPerStixelX, micronsPerStixelY,
                xOffset, yOffset, interval, monitorFrequency, framesPerTTL,
                refreshPeriod, nFramesRequired, droppedFrames);
    }
    
    public class RunTimeMovieParams{
        public int pixelsPerStixelX, pixelsPerStixelY;
        public double width, height; //stixels
        public double micronsPerStixelX, micronsPerStixelY; 
        public double xOffset, yOffset; //microns, from bottom left
        public int interval;
        public double monitorFrequency;
        public int framesPerTTL;
        public double refreshPeriod; //ms
        public int nFramesRequired;
        public int[] droppedFrames;
        
        // For writing
        RunTimeMovieParams(int pixelsPerStixelX, int pixelsPerStixelY, double width, 
                double height, double micronsPerStixelX, double micronsPerStixelY,
                double xOffset, double yOffset, int interval, double monitorFrequency, int framesPerTTL,
                double refreshPeriod, int nFramesRequired, int[] droppedFrames) throws IOException {
            
            ByteArrayOutputStream bos = new ByteArrayOutputStream(0);
            DataOutputStream dos = new DataOutputStream(bos);
            dos.writeInt(pixelsPerStixelX);
            dos.writeInt(pixelsPerStixelY);
            dos.writeDouble(width);
            dos.writeDouble(height);
            dos.writeDouble(micronsPerStixelX);
            dos.writeDouble(micronsPerStixelY);
            dos.writeDouble(xOffset);
            dos.writeDouble(yOffset);
            dos.writeInt(interval);
            dos.writeDouble(refreshPeriod);
            dos.writeInt(nFramesRequired);
            dos.writeInt(droppedFrames.length);
            for (int i = 0; i < droppedFrames.length; i++)
                dos.writeInt(droppedFrames[i]);
            addChunk(RTMP_TAG, bos.toByteArray());	
            
            ByteArrayOutputStream mfBOS = new ByteArrayOutputStream(0);
            DataOutputStream mfDOS = new DataOutputStream(mfBOS);
            mfDOS.writeDouble(monitorFrequency);
            mfDOS.writeInt(framesPerTTL);
            addChunk(MONITOR_FREQ_TAG, mfBOS.toByteArray());

            if (DEBUG) System.out.println(toString());
        }
        
        // For reading
        RunTimeMovieParams() throws IOException {
            DataInputStream dis = new DataInputStream(new ByteArrayInputStream(getChunk(RTMP_TAG)));
            pixelsPerStixelX = dis.readInt();
            pixelsPerStixelY = dis.readInt();
            width = dis.readDouble();
            height = dis.readDouble();
            micronsPerStixelX = dis.readDouble();
            micronsPerStixelY = dis.readDouble();
            xOffset = dis.readDouble();
            yOffset = dis.readDouble();
            interval = dis.readInt();
            refreshPeriod = dis.readDouble();
            nFramesRequired = dis.readInt();
            
            int nDroppedFrames = dis.readInt();
            droppedFrames = new int[nDroppedFrames];
            for (int i = 0; i < nDroppedFrames; i++)
                droppedFrames[i] = dis.readInt();
            
            monitorFrequency = VisionParams.DEFAULT_MONITOR_FREQUENCY; // Defaults for older globals files
            framesPerTTL = VisionParams.DEFAULT_FRAMES_PER_TTL; 	   // 
            byte[] mfBuffer = getChunk(MONITOR_FREQ_TAG);
            if (mfBuffer != null) {
                DataInputStream mfDIS = new DataInputStream(new ByteArrayInputStream(mfBuffer));
                monitorFrequency = mfDIS.readDouble();
                framesPerTTL = mfDIS.readInt();
            }
            
            if (DEBUG) System.out.println(toString());
        }
        
        /**
         * RunTimeMovieParams
         */
        public String toString() {
            String newline = System.getProperty("line.separator");
            String toReturn = 
                "pixelsPerStixelX: "  + pixelsPerStixelX  + newline +
                "pixelsPerStixelY: "  + pixelsPerStixelY  + newline +
                "width: " 			  + width             + newline + 
                "height: " 			  + height            + newline +
                "micronsPerStixelX: " + micronsPerStixelX + newline + 
                "micronsPerStixelY: " + micronsPerStixelY + newline +
                "xOffset: " 		  + xOffset           + newline + 
                "yOffset: " 		  + yOffset           + newline + 
                "interval: " 		  + interval          + newline + 
                "monitorFrequency: "  + monitorFrequency  + newline +
                "framesPerTTL: "      + framesPerTTL      + newline +
                "refreshPeriod: "     + refreshPeriod     + newline +
                "nFramesRequired: "   + nFramesRequired   + newline + 
                "droppedFrames:" 	  + newline;
            
            for (int i = 0; i < droppedFrames.length; i++)
                toReturn = toReturn + "   " + droppedFrames[i] + newline;
            
            return toReturn;
        }
        
        /**
         * In units for STA plotting
         * @return
         */
        public double[] getSTARange() {
            return new double[]{xOffset/1000, (micronsPerStixelX*width  + xOffset)/1000,
                                yOffset/1000, (micronsPerStixelY*height + yOffset)/1000};
        }

    }
    
    public int getCreatedBy() throws IOException {
        byte[] buffer = getChunk(CREATED_BY_TAG);
        if (buffer == null) return 0;
        
        DataInputStream dis = new DataInputStream(new ByteArrayInputStream(buffer));		
        return dis.readInt();
    }
    
    public void writeCreatedBy(int createdBy) throws IOException {
        ByteArrayOutputStream bos = new ByteArrayOutputStream(0);
        DataOutputStream dos = new DataOutputStream(bos);
        dos.writeInt(createdBy);
        addChunk(CREATED_BY_TAG, bos.toByteArray());	
    }
    
    
    public int getVisionVersion() throws IOException {
        byte[] buffer = getChunk(VERSION_TAG);
        if (buffer == null) return 0;
        
        DataInputStream dis = new DataInputStream(new ByteArrayInputStream(buffer));
        return dis.readInt();
    }

    public void writeVisionVersion(int version) throws IOException {
        ByteArrayOutputStream bos = new ByteArrayOutputStream(0);
        DataOutputStream dos = new DataOutputStream(bos);
        dos.writeInt(version);
        addChunk(VERSION_TAG, bos.toByteArray());	
    }
        
}
