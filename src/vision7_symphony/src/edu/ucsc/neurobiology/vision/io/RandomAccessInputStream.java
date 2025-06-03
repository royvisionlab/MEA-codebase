package edu.ucsc.neurobiology.vision.io;

import java.io.*;


/**
 * Adapter to transform a RandomAccessFile into an OutputStream.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class RandomAccessInputStream
    extends InputStream {

    private RandomAccessFile f;


    public RandomAccessInputStream(RandomAccessFile f) {
        this.f = f;
    }


    public int read() throws java.io.IOException {
        return f.read();
    }


    public long skip(long n) throws java.io.IOException {
        return f.skipBytes( (int) n);
    }


    public int read(byte[] b, int off, int len) throws java.io.IOException {
        return f.read(b, off, len);
    }


    public int read(byte[] b) throws java.io.IOException {
        return f.read(b);
    }


    public synchronized void reset() throws java.io.IOException {
        throw new Error("Method not implemented");
    }


    public synchronized void mark(int readlimit) {
        throw new Error("Method not implemented");
    }


    public void close() throws java.io.IOException {
        throw new Error("Method not implemented");
    }


    public boolean markSupported() {
        throw new Error("Method not implemented");
    }


    public int available() throws java.io.IOException {
        throw new Error("Method not implemented");
    }
}
