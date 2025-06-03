package edu.ucsc.neurobiology.vision.test;

import java.io.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.io.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public final class HuffmanSampleInputStream
    extends Thread {

    private BitInputStream bitInput;
    private SampleListener sampleListener1, sampleListener2, sampleListener3;
    private final int nElectrodes;
    private int samplesRead;

    HuffmanCoding.Node[][] nodes;
    int[] rootNodeIndex;


    /**
     * Creates a HuffmanSampleInputStream which will read data from a specific imput
     * stream.
     *
     * @param input The InputStream which will be used as a source of data
     * @param nElectrodes The number of electrodes the raw data are comming from
     */
    public HuffmanSampleInputStream(
        InputStream input, int nElectrodes) throws IOException {

        super("HuffmanSampleInputStream Thread");
        this.setPriority(Thread.MAX_PRIORITY);
        this.nElectrodes = nElectrodes;

        // read the Huffman tree
        DataInputStream dis = new DataInputStream(input);
        nodes = new HuffmanCoding.Node[nElectrodes]
                [ (HuffmanCoding.max - HuffmanCoding.min) * 2 - 1];
        rootNodeIndex = new int[nElectrodes];
        for (int electrode = 0; electrode < nElectrodes; electrode++) {
            rootNodeIndex[electrode] = dis.readShort();

//            System.out.println(rootNodeIndex[electrode] + " root");

            for (int i = 0; i < nodes[electrode].length; i++) {
                nodes[electrode][i] = new HuffmanCoding.Node(
                    -1, i - 2048, dis.readShort(), dis.readShort());

//               System.out.println(nodes[electrode][i].leftIndex + " left");
//               System.out.println(nodes[electrode][i].leftIndex + " right");
            }
        }

        this.bitInput = new BitInputStream(input);

        System.out.println("constucted");
    }


    /**
     * Registers a new DataListener to this Input Stream. This function should be
     * called before the stream is started.
     */
    public void addSampleListener(SampleListener listener) {
        if (sampleListener1 == null) {
            sampleListener1 = listener;
        } else if (sampleListener2 == null) {
            sampleListener2 = listener;
        } else if (sampleListener3 == null) {
            sampleListener3 = listener;
        } else {
            throw new IllegalArgumentException("Too many SampleListeners (3 allowed)");
        }
    }


    private final void endOfData() {
        if (sampleListener1 != null) {
            try {
                sampleListener1.finishSampleProcessing();
            } catch (IOException e) {
                Vision.reportException("Could not finish processing for " +
                                       sampleListener1.getClass().getName(), e);
            }
        }
        if (sampleListener2 != null) {
            try {
                sampleListener2.finishSampleProcessing();
            } catch (IOException e) {
                Vision.reportException("Could not finish processing for " +
                                       sampleListener2.getClass().getName(), e);
            }
        }
        if (sampleListener3 != null) {
            try {
                sampleListener3.finishSampleProcessing();
            } catch (IOException e) {
                Vision.reportException("Could not finish processing for " +
                                       sampleListener3.getClass().getName(), e);
            }
        }

        // close the stream
        try {
            bitInput.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    public synchronized int getSamplesRead() {
        return samplesRead;
    }


    /**
     * Expanding compressed data is a little harder than the compression
     * phase.  As each new symbol is decoded, the tree is traversed,
     * starting at the root node, reading a bit in, and taking either the
     * iLeftChild or iRightChild path.  Eventually, the tree winds down to a
     * leaf node, and the corresponding symbol is output.  If the symbol
     * is the END_OF_STREAM symbol, it doesn't get written out, and
     * instead the whole process terminates.
     *
     * @throws IOException
     */
    public void run() {
        short[] sample = new short[nElectrodes];
        int currentNode;

        System.out.println("running...");

        while (true) {
//            System.out.println("sample " + samplesRead);

            for (int electrode = 0; electrode < nElectrodes; electrode++) {
                currentNode = rootNodeIndex[electrode];

                // We will read the input file bit by bit and move the node pointer one
                // level in the tree for each bit read.  We will stop when we get to a
                // leaf node, write the uncompressed character to output, reset the node
                // pointer to the root, and repeat the process.
                do {
                    // If we read a 1 bit then traverse right else
                    // it's a "0" so left in the Huffman tree..
                    try {
                        if (bitInput.readBitFlag()) {
                            currentNode = nodes[electrode][currentNode].rightIndex;
                        } else {
                            currentNode = nodes[electrode][currentNode].leftIndex;
                        }
                    } catch (IOException e) {
                        e.printStackTrace();
                    }

//                    System.out.println(nodes[electrode][currentNode].leftIndex);
//                    System.out.println("--");
                } while (nodes[electrode][currentNode].leftIndex != -1);
                // i.e. while we're pointing to an internal node.. stop at a leaf.

                //System.out.println("electrode " + electrode);

                sample[electrode] = (short) nodes[electrode][currentNode].value;

//                System.out.println(nodes[][currentNode].value + " read");
            }

            // tell the listeners we have a sample
            if (sampleListener1 != null) {
                sampleListener1.processSample(sample);
            }
            if (sampleListener2 != null) {
                sampleListener2.processSample(sample);
            }
            if (sampleListener3 != null) {
                sampleListener3.processSample(sample);
            }

            samplesRead++;
        }

//        while (true) {
        // shift the remaining bytes to the begining of the buffer
//            System.arraycopy(buffer, bytesUsed, buffer, 0, bytesAvailable - bytesUsed);
//            bytesAvailable -= bytesUsed;

        // read data from the underlying stream until we have at least one full sample
//            do {
//                try {
//                    lastBytesRead =
//                        input.read(buffer, bytesAvailable, buffer.length - bytesAvailable);

        // if we reached the end of the stream: 1) anounce the end of data;
        // and 2) renurn from the run() method (the thread will dye).
        // The thread ends cleanly reading up to the last COMPLETE sample.
//                    if (lastBytesRead <= 0) {
//                        endOfData();
//                        return;
//                    }

//                    bytesAvailable += lastBytesRead;
//                    nAvailableSamples = bytesAvailable / sampleSize;
//                } catch (IOException e) {
        // if an error happens: 1) anounce the end of data; 2) inform the
        // user about the error and 3) renurn from the run() method
        // (the thread will dye).
//                    endOfData();
//                    e.printStackTrace();
//                    JOptionPane.showMessageDialog(
//                        null, e, "IOException in HuffmanSampleInputStream",
//                        JOptionPane.ERROR_MESSAGE);
//                    return;
//                }
//            } while (nAvailableSamples == 0);

        // calculate the number of bytes we are going to use now
//            bytesUsed = nAvailableSamples * sampleSize;

        // send the data to SampleListeners
//            byteIndex = 0;
//            int electrode = 0, b1, b2, b3;
//            for (sampleIndex = 0; sampleIndex < nAvailableSamples; sampleIndex++) {
//                if (nElectrodes % 2 == 1) {
//                    // the TTL is not encoded
//                    sample[0] = (short) ( (buffer[byteIndex++] << 8) +
//                                         (buffer[byteIndex++] & 0xff));
//
//                    for (electrode = 1; electrode < nElectrodes; ) {
//                        b1 = buffer[byteIndex++] & 0xff;
//                        b2 = buffer[byteIndex++] & 0xff;
//                        b3 = buffer[byteIndex++] & 0xff;
//
//                        sample[electrode++] = (short) ( (b1 << 4) + (b2 >> 4) - 2048);
//                        sample[electrode++] = (short) ( ( (b2 & 0x0f) << 8) + b3 - 2048);
//                    }
//                } else {
//                    // the TTL is encoded
//                    for (electrode = 0; electrode < nElectrodes; ) {
//                        b1 = buffer[byteIndex++] & 0xff;
//                        b2 = buffer[byteIndex++] & 0xff;
//                        b3 = buffer[byteIndex++] & 0xff;
//
//                        sample[electrode++] = (short) ( (b1 << 4) + (b2 >> 4) - 2048);
//                        sample[electrode++] = (short) ( ( (b2 & 0x0f) << 8) + b3 - 2048);
//                    }
//                }
//
//            if (sampleListener1 != null) {
//                sampleListener1.processSample(sample);
//            }
//            if (sampleListener2 != null) {
//                sampleListener2.processSample(sample);
//            }
//            if (sampleListener3 != null) {
//                sampleListener3.processSample(sample);
//            }
//
//            samplesRead++;
//            }

//        } //while(true)
    }


    public static void main(String[] args) throws IOException {
        InputStream fis = new BufferedInputStream(new FileInputStream(
            "c:\\vision4\\1.bin"), 1 * 1024 * 1024);
        RawDataHeader header = new RawDataHeader512(fis);

//        DataInputStream dis = new DataInputStream(fis);
//        System.out.println(dis.readShort());
//        System.out.println(dis.readInt());
//        System.exit(1);

        final int nSamples = header.getNumberOfSamples();

        HuffmanSampleInputStream h = new HuffmanSampleInputStream(fis, 129);
        h.addSampleListener(new SampleListener() {
            long sampleIndex, percentage, oldPercentage;

            public void processSample(short[] sample) {
                // calculate progress
                sampleIndex++;
                percentage = 100 * sampleIndex / nSamples;
                if (percentage != oldPercentage) {
                    System.out.println(percentage + "%");
                    oldPercentage = percentage;
                }
            }


            public void finishSampleProcessing() {
            }
        });
        h.start();
    }
}
