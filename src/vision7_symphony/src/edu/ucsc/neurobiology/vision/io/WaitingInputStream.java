package edu.ucsc.neurobiology.vision.io;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

/**
 * This is an InputStream that waits for a file to grow, rather than exiting when EOF is reached.
 * 
 * @author Peter H. Li, The Salk Institute
 */
public class WaitingInputStream extends InputStream {
    private final InputStream input;
    private final long waitMillis = 3000; // This could be made user-facing if necessary
    private final int iterations = 3;     // This could be made user-facing if necessary

    public WaitingInputStream(InputStream input) {
        if (!(input instanceof FileInputStream)) 
            throw new Error("WaitingInputStream is only designed to be used for underyling InputStreams of type FileInputStream");

        this.input = input;
    }

    @Override
    public synchronized int read() throws IOException {
        throw new Error("WaitingInputStream#read() with no arguments has not been tested; see source for testable implementation");

        // Untested, but should work
//		byte[] b = new byte[1];
//		if (read(b, 0, 1) == -1) return -1;
//		return (int) b[0];
    }
    
    @Override
    public synchronized int read(byte[] b) throws IOException {
        return read(b, 0, b.length);
    }
    
    @Override
    public synchronized int read(byte[] b, int off, int len) throws IOException {
        int bytesRead = 0;
        
        for (int i = 0; i < iterations; i++) {
            bytesRead = input.read(b, off, len);
            if (bytesRead != -1) break; // read succeeded; we're done.

            // Otherwise wait for more data to come in and try again
            try {
                Thread.sleep(waitMillis);
            } catch (InterruptedException e) { e.printStackTrace();	}
        }

        return bytesRead;
    }

}