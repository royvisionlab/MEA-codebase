package edu.ucsc.neurobiology.vision.io;

import java.io.*;

import edu.ucsc.neurobiology.vision.stimulus.*;


/*This class creates or opens a raw movie file.
 * This format is characterized by an ascii header followed by
 * three bytes per stixel to give the color information. It does not
 * depend on the algorithm used to create the frames.
 * Extend this class to add information to the header.  See PinkMovieFile
 * for an example.
 File will look like:
 width + fieldSeparator + value + newline
 height + field Separator + value + newline
 frames-generated + fieldSeparator + value + newline //OR frames + fieldSeparator + value + newline
 algorithm + field Separator + value + newline
 additional lines, that depend upon the movie in question.
 newline
 frame1 + frame2 + .... frameN

 Values may be ints or doubles, and may be an array.
 If they are boolean, the should have field values TRUE or FALSE.
 The + symbol is not reproduced in the output.
 fieldSeparator = '\t' = the standard tab.
 newline = '\r\n', '\r', or '\n'  The reader must determine which was used.
 Each frame has the format: (0,0) + (0,1) + ... (0, width-1) + (1, 0) + (1, 1) +
     (height - 1, width-1).    (0,0) is the upper left.  Each (#,#) is three bytes,
     with values Red, Green, Blue, and can range from 0 to 255.  0 is black.  255
     is maximum.
 Sample:
 width	4
 height	2
 frames-generated	2  //OR frames	2
 algorithm	PowerLaw-v1
 seed	11111
 independent	TRUE
 distribution	GAUSSIAN
 standard-deviation	0.16
 white-spatial-period-cutoff	4
 white-temporal-period-cutoff	4
 spatial-power	-1.0
 temporal-power	-2.0

 binary data.....

 * @author Matthew Grivich, The Salk Institute
 */
public class RawMovieFile
    extends AsciiHeaderFile {


    public String algorithm;
    private byte[] buffer;


    //This constructor is used for creating a file
    public RawMovieFile(
        String fileName, int width, int height, int framesGenerated, String algorithm) throws
        IOException {

        super(fileName, true);

        this.algorithm = algorithm;

        addToHeader("width", width);
        addToHeader("height", height);
        addToHeader("frames-generated", framesGenerated);
        addToHeader("algorithm", algorithm);
        this.buffer = new byte[width * height * 3];
    }


    //This constructor is used for creating a file
    public RawMovieFile(
        String fileName, int width, int height, int framesGenerated, String algorithm,
        int dummy, String ...params) throws IOException {

        super(fileName, true);

        this.algorithm = algorithm;

        addToHeader("width", width);
        addToHeader("height", height);
        addToHeader("frames-generated", framesGenerated);
        addToHeader("algorithm", algorithm);

        for (String p : params) {
            addToHeader(p);
        }

        writeHeader();
        this.buffer = new byte[width * height * 3];
    }


    //This constructor is used for opening an existing file.
    public RawMovieFile(String fileName) throws IOException {
        super(fileName, false);
        this.buffer = new byte[getWidth() * getHeight() * 3];
    }


    public int getWidth() {
        return (new Integer(getValues("width")[0])).intValue();
    }


    public int getHeight() {
        return (new Integer(getValues("height")[0])).intValue();
    }


    public int getFramesGenerated() {
        //added extra branch to support ej's file format fork.  Done over my protests.  --mgrivich
        int frames = -1;
        try {
            frames = new Integer(getValues("frames")[0]).intValue();
        } catch(Exception e) {
            frames = new Integer(getValues("frames-generated")[0]).intValue();
        }
        return frames;
    }


    public String getAlgorithm() {
        return getValues("algorithm")[0];
    }


    synchronized public void writeFrame(float[] colors) throws IOException {
//        file.seek(0);

        ByteArrayOutputStream bos = new ByteArrayOutputStream();
        DataOutputStream dos = new DataOutputStream(bos);

        for (int i = 0; i < colors.length; i++) {
            //Writes the bytes as unsigned, even though unsigned does not exist in java.
            //The bits show up in the right place.
            float value = Math.round( (257.0 * colors[i]) - 1.0);
            if (value < 0) {
                value = 0;
            }
            if (value > 255) {
                value = 255;
            }
            dos.writeByte( (byte) value);
        }

        dos.close();
        bos.close();
        file.write(bos.toByteArray());
    }


    //read a frame.  Return it as an array of floats
    synchronized public float[] readFrameFloat(int frame) throws IOException {
        float[] colors = new float[buffer.length];
        file.seek( (long) headerLength + (long) buffer.length * (long) frame);
        file.readFully(buffer);

        //read byte integers as signed, even though they are really unsigned,
        //because java does not have unsigned integers.
        for (int i = 0; i < colors.length; i++) {
            byte b = buffer[i];
            if (b >= 0) {
                colors[i] = (b + 1.000f) / 257.0f;
            } else {
                colors[i] = (b + 257.0f) / 257.0f;
            }
        }
        return colors;
    }


    synchronized public int[] readFrameInt(int frame) throws IOException {
        int[] colors = new int[buffer.length];
        file.seek( (long) headerLength + (long) buffer.length * (long) frame);
        file.readFully(buffer);

        //read byte integers as signed, even though they are really unsigned,
        //because java does not have unsigned integers.
        for (int i = 0; i < colors.length; i++) {
            if (buffer[i] >= 0) {
                colors[i] = buffer[i];
            } else {
                colors[i] = buffer[i] + 256;
            }
        }

        return colors;
    }


    synchronized public long[][] readFrameLong(int frame, boolean calculateSTV) throws
        IOException {

        file.seek( (long) headerLength + (long) buffer.length * (long) frame);
        file.readFully(buffer);

        //read byte integers as signed, even though they are really unsigned,
        //because java does not have unsigned integers.
        long[][] colors = new long[2][buffer.length / 3]; // i.e. w*h

        for (int i = 0, j = 0; i < colors[0].length; i++) {
            long r = (buffer[j] >= 0) ? buffer[j] : buffer[j] + 256;
            j++;
            long g = (buffer[j] >= 0) ? buffer[j] : buffer[j] + 256;
            j++;
            long b = (buffer[j] >= 0) ? buffer[j] : buffer[j] + 256;
            j++;

            colors[0][i] =
                (b << GaussFrameGenerator.TWICE_BITS_PER_COLOR) +
                (g << GaussFrameGenerator.BITS_PER_COLOR) +
                r;

            if (calculateSTV) {
                colors[1][i] =
                    ( ( (b * b)) << GaussFrameGenerator.TWICE_BITS_PER_COLOR) +
                    ( ( (g * g)) << GaussFrameGenerator.BITS_PER_COLOR) +
                    ( (r * r));
            }

//            colors[1][i] =
//                ( ( (b * b) >> 8) << GaussFrameGenerator.TWICE_BITS_PER_COLOR) +
//                ( ( (g * g) >> 8) << GaussFrameGenerator.BITS_PER_COLOR) +
//                ( (r * r) >> 8);

        }

        return colors;
    }


    synchronized public void close() throws IOException {
        file.close();
    }

}
