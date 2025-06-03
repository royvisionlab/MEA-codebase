package edu.ucsc.neurobiology.vision.io;

import java.io.*;


/**
 * This class allows reading and writing raw data in electrode major format. After a standard
 * 512 header the data from each electrode is appended after the previous electode. The data
 * is compressed. Used by the snf (Serial Neuron Finder) package.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ElectrodeMajorRawDataFile {
    private RandomAccessFile rawDataFile;
    private RawDataHeader header;
    private final int nSamples;
    private byte[] byteBuffer;
    private final int blockSize = 1 * 1024;


    public ElectrodeMajorRawDataFile(String rawFileName) throws IOException {
        FileInputStream is = null;
        try {
            is = new FileInputStream(rawFileName);
            header = new RawDataHeader512(is);
        } catch (IOException e) {
            if (is != null) {
                is.close();
            }
            throw e;
        }
        nSamples = header.getNumberOfSamples();
        byteBuffer = new byte[ (int) (3L * nSamples / 2L)];

        rawDataFile = new RandomAccessFile(rawFileName, "rw");
    }


    public RawDataHeader getHeader() {
        return header;
    }


    public void writeRawData(int electrode, short[] data, int length) throws IOException {
        rawDataFile.seek(
            (long) header.getHeaderSize() + 3L * electrode * nSamples / 2L);
        for (int i = 0, byteIndex = 0; i < length; ) {
            int s1 = data[i++] + 2048;
            int s2 = data[i++] + 2048;
            byteBuffer[byteIndex++] = (byte) (s1 >> 4);
            byteBuffer[byteIndex++] = (byte) ( ( (s1 & 0x000f) << 4) + (s2 >> 8));
            byteBuffer[byteIndex++] = (byte) (s2 & 0x00ff);
        }
        final int n = byteBuffer.length / blockSize;
        int offset = 0;
        for (int i = 0; i < n; i++) {
            rawDataFile.write(byteBuffer, offset, blockSize);
            offset += blockSize;
        }
        rawDataFile.write(byteBuffer, offset, byteBuffer.length - n * blockSize);
    }


    /*
        // written by Matt to fix the memory bug
     public void readRawData(int electrode, short[] data, int length) throws IOException {
            int samplesLeft = length;
            //Reading data in blocks prevents an out of memory hour with 1 hr long data sets.
            int samplesAtATime = 2000000;
            int dataIndex = 0;
            boolean dataLeft = true;
            Arrays.fill(data, (short) 0);
            rawDataFile.seek(
                (long) header.getHeaderSize() + 3L * electrode * nSamples / 2L);

            while (dataLeft == true) {
                int samplesToRead = Math.min(samplesAtATime, samplesLeft);
                samplesLeft -= samplesToRead;
                rawDataFile.readFully(byteBuffer, 0, (int) (3L * samplesToRead / 2L));
                int b1, b2, b3;
                for (int byteIndex = 0, i = 0; i < samplesToRead; ) {
                    b1 = byteBuffer[byteIndex++] & 0xff;
                    b2 = byteBuffer[byteIndex++] & 0xff;
                    b3 = byteBuffer[byteIndex++] & 0xff;
                    data[dataIndex++] = (short) ( (b1 << 4) + (b2 >> 4) - 2048);
                    data[dataIndex++] = (short) ( ( (b2 & 0x0f) << 8) + b3 - 2048);
                    i += 2;
                }
                if (samplesToRead == 0) {
                    dataLeft = false;
                }
            }
        }
     */


    /*
     public void readRawData(int electrode, short[] data, int length) throws IOException {
        rawDataFile.seek(
            (long) header.getHeaderSize() + 3L * electrode * nSamples / 2L);
        rawDataFile.readFully(byteBuffer, 0, (int) (3L * length / 2L));

        Arrays.fill(data, (short) 0);
        int b1, b2, b3;
        for (int i = 0, byteIndex = 0; i < length; ) {
            b1 = byteBuffer[byteIndex++] & 0xff;
            b2 = byteBuffer[byteIndex++] & 0xff;
            b3 = byteBuffer[byteIndex++] & 0xff;
            data[i++] = (short) ( (b1 << 4) + (b2 >> 4) - 2048);
            data[i++] = (short) ( ( (b2 & 0x0f) << 8) + b3 - 2048);
        }
     }
     */


    public void readRawData(int electrode, short[] data, int samplesToRead) throws
        IOException {

        rawDataFile.seek(
            (long) header.getHeaderSize() + 3L * electrode * nSamples / 2L);
        int bytesToRead = (int) (3L * samplesToRead / 2L);
        final int n = bytesToRead / blockSize;
        int offset = 0;
        for (int i = 0; i < n; i++) {
            rawDataFile.readFully(byteBuffer, offset, blockSize);
            offset += blockSize;
        }
        rawDataFile.readFully(byteBuffer, offset, bytesToRead - n * blockSize);

        int b1, b2, b3;
        for (int i = 0, byteIndex = 0; i < nSamples; ) {
            b1 = byteBuffer[byteIndex++] & 0xff;
            b2 = byteBuffer[byteIndex++] & 0xff;
            b3 = byteBuffer[byteIndex++] & 0xff;
            data[i++] = (short) ( (b1 << 4) + (b2 >> 4) - 2048);
            data[i++] = (short) ( ( (b2 & 0x0f) << 8) + b3 - 2048);
        }
    }


    /*
        public void benchmarkReadWrite() throws IOException {
            double t1, t2, tt1, tt2, readT = 0, writeT = 0;
            short[] data = new short[nSamples];
            int n = header.getNumberOfElectrodes();
            tt1 = System.currentTimeMillis();
            for (int electrode = 1; electrode < n; electrode++) {
                t1 = System.currentTimeMillis();
                readRawData(electrode, data);
                t2 = System.currentTimeMillis();
                readT += (t2 - t1);
                t1 = System.currentTimeMillis();
                writeRawData(electrode, data);
                t2 = System.currentTimeMillis();
                writeT += (t2 - t1);
                System.out.println(electrode);
            }
            tt2 = System.currentTimeMillis();
            double t = (tt2 - tt1) / 1000;
            readT /= 1000;
            writeT /= 1000;
            System.out.println("Read Speed: " +
                               1.0 * n * byteBuffer.length / 1024.0 / 1024.0 / readT);
            System.out.println("Write Speed: " +
                               1.0 * n * byteBuffer.length / 1024.0 / 1024.0 / writeT);
            System.out.println("Overall: " +
                               2.0 * n * byteBuffer.length / 1024.0 / 1024.0 / t +
                               ", " + t + " seconds.");
        }
     */

    public void close() throws IOException {
        rawDataFile.close();
    }

}
