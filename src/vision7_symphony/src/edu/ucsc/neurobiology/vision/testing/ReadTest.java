package edu.ucsc.neurobiology.vision.testing;

import java.io.*;

import edu.ucsc.neurobiology.vision.*;


/**
 * @author nobody, anyone can change
 */
public class ReadTest
    extends Thread {

    final long nBytes;
    final int bufferSize;
    byte[] buffer;
    int N;
    RandomAccessFile fos;
    String fileName;


    public ReadTest(String fileName, long nBytes, int bufferSize) throws IOException {
        this.nBytes = nBytes;
        this.bufferSize = bufferSize;
        buffer = new byte[bufferSize];
        N = (int) (nBytes / bufferSize);
        fos = new RandomAccessFile(fileName, "r");
        this.fileName = fileName;
    }


    public void run() {
        try {
            long t1 = System.currentTimeMillis();
            for (int i = 0; i < N; i++) {
                fos.read(buffer);
            }
            long t2 = System.currentTimeMillis();
            System.out.println("Read Performance : " + fileName + "  : " +
                               (nBytes / 1024 / 1024) * 1000.0 / (t2 - t1));
            fos.close();
        } catch (IOException e) {
            Vision.reportFatalException("", e);
        }
    }


    public static void main(String[] args) throws Exception {
        ReadTest t1 = new ReadTest(
            "d:\\2005-04-26-2\\data002\\data002.rem",
            100 * 1024 * 1024, 2 * 1024);
        t1.start();
    }

}
