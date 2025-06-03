package edu.ucsc.neurobiology.vision.io;

import java.io.*;
import java.util.*;


/**
 * This class implements a Java InputStream that reads data from an underlying
 * InputStream asynchronously by doing all reads in another thread.
 * The data is also buffered after being read.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public final class AsynchronousInputStream
    extends InputStream implements Runnable {

    private final InputStream input;
    private LinkedList<Buffer> emptyList;
    private LinkedList<Buffer> fullList;
    private final int nBuffers;
    private final int bufferSize;
    private Buffer readBuffer, writeBuffer;
    private boolean finished;
    private boolean alive = true;
    private IOException exceptionHappened = null;


    public AsynchronousInputStream(InputStream input, int bufferSize, int nBuffers) throws
        IOException {

        this.input = input;
        this.bufferSize = bufferSize;
        this.nBuffers = nBuffers;

        this.emptyList = new LinkedList<Buffer> ();
        this.fullList = new LinkedList<Buffer> ();
        for (int i = 0; i < nBuffers; i++) {
            emptyList.add(new Buffer());
        }
    }


    public void run() {
        int bytesRead;
        run:while (alive) {
        
            // if there is no current buffer for writing: wait until there is a buffer
            // in the empty list, extract it and make it current
            if (writeBuffer == null) {
                
                synchronized (this) {
                    if (!alive) {
                        break run;
                    }
                    while (emptyList.isEmpty()) {
                        notifyAll();
                        try {
                            wait();
                        } catch (InterruptedException e) {
                            e.printStackTrace();
                        }
        
                    }
  

                    writeBuffer = emptyList.removeFirst();
                    writeBuffer.pointer = 0;
                    writeBuffer.readLimit = 0;
                }
            }

            // try to fill the buffer with data
            if (writeBuffer.writeLimit - writeBuffer.pointer <= 0) {
                throw new Error("aaa");
            }
            try {
                bytesRead = input.read(writeBuffer.bytes, writeBuffer.pointer,
                                       writeBuffer.writeLimit - writeBuffer.pointer);
            } catch (IOException e) {
                synchronized (this) {
                    exceptionHappened = e;
                }
                break run;
            }

            if (bytesRead == -1 ) {
                // we have reached the end of the stream, insert the current buffer
                // as it is to the end of the full list
                synchronized (this) {
                    finished = true;
                    fullList.addLast(writeBuffer);
                    notifyAll();
                }

                break run;
            } else {
                // we haven't read anything, go back and try again
                if (bytesRead == 0) {      
                    continue;
                }
            }

            // we have read something, increment the pointer and the readLimit
            writeBuffer.pointer += bytesRead;
            writeBuffer.readLimit += bytesRead;
            // if we have filled the buffer with data, append the current buffer to the
            // end of the full list and cancel the current buffer
            if (writeBuffer.pointer == writeBuffer.writeLimit) {
                synchronized (this) {
                    fullList.addLast(writeBuffer);
                    writeBuffer = null;
                    notifyAll();
                }
            }
        } // while
    }


    public synchronized int read() throws IOException {
        throw new IOException("Method not supported");
    }


    public synchronized int read(byte[] b) throws IOException {
        return read(b, 0, b.length);
    }


    public synchronized int read(byte[] b, int offset, int length) throws IOException {
        if (exceptionHappened != null) {
            System.err.println("THIS EXCEPTION IS THROWN FROM THE read(...) METHOD.");
            throw exceptionHappened;
        }

        if (finished && fullList.isEmpty() && readBuffer == null) {
            return -1;
        }

        int bytesToRead;
        if (readBuffer == null) {
            // there is no current buffer for reading: wait until there is one and
            // make it current
            while (fullList.isEmpty()) {
                notifyAll();
                try {
                    wait();
                } catch (InterruptedException e) {}
            }
            readBuffer = fullList.removeFirst();
            readBuffer.pointer = 0;
        }

        // copy the requested data to the user's buffer
        bytesToRead = Math.min(readBuffer.readLimit - readBuffer.pointer, length);
        System.arraycopy(readBuffer.bytes, readBuffer.pointer, b, offset, bytesToRead);
        // increment the pointer
        readBuffer.pointer += bytesToRead;
        // if we have read all the buffer's content attach it to the empty list and
        // cancel it
        if (readBuffer.pointer == readBuffer.readLimit) {
            emptyList.addLast(readBuffer);
            readBuffer = null;
            notifyAll();
        }

        return bytesToRead;
    }


    public synchronized int getBufferUsage() {
        return 100 * fullList.size() / nBuffers;
    }


    public synchronized int available() throws IOException {
        return readBuffer.readLimit - readBuffer.pointer;
    }


    public synchronized void close() throws IOException {
        alive = false;
       
 //       while (!runExited) {
 //           notifyAll();
 //       }
       
        input.close();

        emptyList.clear();
        fullList.clear();
    }


    public synchronized void mark(int readlimit) {
        throw new IllegalAccessError("Method not supported");
    }


    public synchronized boolean markSupported() {
        return false;
    }


    public synchronized void reset() throws IOException {
        throw new IOException("Method not supported");
    }


    public synchronized long skip(long n) throws IOException {
        throw new IllegalAccessError("Method not supported");
    }


    public void start() {
        Thread t = new Thread(this, "AsynchronousInputStream Thread");
        t.setPriority(Thread.MAX_PRIORITY);
        t.start();
    }


    private class Buffer {
        final byte[] bytes = new byte[bufferSize];
        final int writeLimit = bytes.length;
        int pointer = 0;
        int readLimit = 0;
    }


}
