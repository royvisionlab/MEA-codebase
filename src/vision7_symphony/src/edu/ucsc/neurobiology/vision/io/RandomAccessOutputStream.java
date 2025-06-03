package edu.ucsc.neurobiology.vision.io;

import java.io.*;


/**
 * Adapter to transform a RandomAccessFile into an InputStream.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class RandomAccessOutputStream
    extends OutputStream {

    private RandomAccessFile f;


    public RandomAccessOutputStream(RandomAccessFile f) {
        this.f = f;
    }


    public void write(int b) throws java.io.IOException {
        f.write(b);
    }


    public void write(byte[] b) throws java.io.IOException {
        f.write(b);
    }


    public void flush() throws java.io.IOException {
    }


    public void write(byte[] b, int off, int len) throws java.io.IOException {
        f.write(b, off, len);
    }


    public void close() throws java.io.IOException {
        throw new Error("Method not implemented");
    }

}
