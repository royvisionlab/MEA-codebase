package edu.ucsc.neurobiology.vision.testing;

import java.io.*;

import edu.ucsc.neurobiology.vision.*;


/**
 * @author nobody, anyone can change
 */
public class WriteTest
    extends Thread {

    final long nBytes;
    final int bufferSize;
    byte[] buffer;
    int N;
    RandomAccessFile fos;
    String fileName;

    public WriteTest(String fileName, long nBytes, int bufferSize) throws IOException {
        this.nBytes = nBytes;
        this.bufferSize = bufferSize;
        buffer = new byte[bufferSize];
        N = (int) (nBytes / bufferSize);
        fos = new RandomAccessFile(fileName, "rw");
        this.fileName = fileName;
    }


    public void run() {
        try {
            long t1 = System.currentTimeMillis();
            for (int i = 0; i < N; i++) {
                fos.write(buffer);
            }
            long t2 = System.currentTimeMillis();
            System.out.println("Write Performance : " + fileName + "  : " +
                               (nBytes / 1024 / 1024) * 1000.0 / (t2 - t1));
            fos.close();
        } catch (IOException e) {
            Vision.reportFatalException("", e);
        }
    }


    public static void main(String[] args) throws Exception {
//        WriteTest t1 = new WriteTest(
//            "/private/var/automount/snle/home/snl-e/Desktop/test.test", 1*1024*1024*1024, 1*1024);
        WriteTest t1 = new WriteTest(
            "c:\\test.test", 100 * 1024 * 1024, 512);

//        WriteTest t2 = new WriteTest("f:\\test2");
        t1.start();
//        t2.start();
    }

}
