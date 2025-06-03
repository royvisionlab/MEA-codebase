package edu.ucsc.neurobiology.vision.io;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.*;


/**
 * This class implements a Java OutputStream that writes data to an underlying
 * OutputStream asynchronously by doing all writes in another thread.
 * The data is also buffered prior to writing.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public final class AsynchronousOutputStream
    extends OutputStream implements Runnable {

    private final OutputStream output;
    private final int nBuffers;
    private final int bufferSize;
    private LinkedList<Buffer> emptyList;
    private LinkedList<Buffer> fullList;
    private Buffer writeBuffer;
    private boolean closeCalled;
    private boolean died = false;

    private Exception exceptionOccured = null;


    public AsynchronousOutputStream(
        OutputStream output, int bufferSize, int nBuffers) throws IOException {

        this.output = output;
        this.bufferSize = bufferSize;
        this.nBuffers = nBuffers;

        this.emptyList = new LinkedList<Buffer>();
        this.fullList = new LinkedList<Buffer>();
        for (int i = 0; i < nBuffers; i++) {
            emptyList.add(new Buffer());
        }
    }


    public void run() {
        Buffer saveBuffer = null;

        while (true) {
            if (closeCalled && fullList.isEmpty()) {
                try {
                    output.close();
                } catch (IOException e) {
                    Vision.reportException("Could not close streams. ", e);
                } finally {
                    break;
                }
            }

            // if there is no current buffer to save: wait until there is one
            // in the empty list, extract it and make it current
            if (saveBuffer == null) {
                synchronized (this) {
                    while (fullList.isEmpty()) {
                        try {
                            wait();
                        } catch (InterruptedException e) {}
                    }
                    saveBuffer = (Buffer) fullList.removeFirst();
                }
            }

            // try to write out the buffer
            try {
//                try {
                output.write(saveBuffer.bytes, 0, saveBuffer.pointer);
//                }catch (NullPointerException e) {
//                    System.err.println(output);
//                    System.err.println(saveBuffer);
//                    e.printStackTrace();
//                }
            } catch (IOException e) {
                synchronized (this) {
                    exceptionOccured = e;
                }
                break;
            }

            synchronized (this) {
                emptyList.addLast(saveBuffer);
                saveBuffer = null;
                notifyAll();
            }
        }

        died = true;
        // tell the listeners (the close method) that the thread dies
        synchronized (this) {
            notifyAll();
        }
    }


    public synchronized void write(int i) throws IOException {
      // Convert the integer to a byte array.
      byte[] result = new byte[4];
      result[0] = (byte) (i >> 24);
      result[1] = (byte) (i >> 16);
      result[2] = (byte) (i >> 8);
      result[3] = (byte) (i /*>> 0*/);
      write(result);

        // throw new IllegalArgumentException("Method not implemented");

    }

    public byte[] toByteArray(int value) {
      return new byte[] {
              (byte)(value >> 24),
              (byte)(value >> 16),
              (byte)(value >> 8),
              (byte)value};
    }


    public synchronized void write(byte[] b) throws IOException {
        write(b, 0, b.length);
    }


    synchronized public void write(byte[] b, int offset, int length) throws IOException {
        if (closeCalled) {
            throw new IOException("The stream was closed");
        }

        if (exceptionOccured != null) {
            throw new ConditionedIOException(exceptionOccured.getMessage(),
                                             exceptionOccured);
        }

        while (length != 0) {
            if (writeBuffer == null) {
                // there is no current buffer for writing: wait until there is one and
                // make it current
                while (emptyList.isEmpty()) {
                    try {
                        wait();
                    } catch (InterruptedException e) {}
                }
                writeBuffer = (Buffer) emptyList.removeFirst();
                writeBuffer.pointer = 0;
            }

            // copy the user's data to the write buffer
            int bytesToCopy = Math.min(
                writeBuffer.writeLimit - writeBuffer.pointer, length);
            System.arraycopy(b, offset, writeBuffer.bytes, writeBuffer.pointer,
                             bytesToCopy);

            // increment the pointers
            offset += bytesToCopy;
            length -= bytesToCopy;
            writeBuffer.pointer += bytesToCopy;

            // if we have filled the buffer attach it to the full list and cancel it
            if (writeBuffer.pointer == writeBuffer.writeLimit) {
                fullList.addLast(writeBuffer);
                writeBuffer = null;
                notifyAll();
            }
        }
    }


    public double getBufferUsage() {
        return 100.0 * fullList.size() / nBuffers;
    }


    public synchronized void flush() throws IOException {
    }


    public synchronized void close() throws IOException {
        fullList.addLast(writeBuffer);
        closeCalled = true;
        notifyAll();

        // wait until the thread dies
        while (!died) {
            try {
                wait();
            } catch (InterruptedException e) {}
        }
    }


    public void start() {
        Thread t = new Thread(this, "BufferingOutputStream Thread");
        t.setPriority(Thread.MAX_PRIORITY);
        t.start();
    }


    class Buffer {
        final byte[] bytes = new byte[bufferSize];
        final int writeLimit = bytes.length;
        int pointer = 0;
    }


}
