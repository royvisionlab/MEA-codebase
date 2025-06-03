package edu.ucsc.neurobiology.vision.test;

import java.io.*;

import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class HuffmanCoding
    implements SampleListener {

    public static int min = -2048, max = 2048;

    SynchronizationObject sync;
    int nElectrodes;
    long sampleIndex, percentage, oldPercentage = -1;
    int nSamples;
    BitOutputStream out;

    Node[][] nodes;
    int[] nodesSize;
    int[] rootNodeIndex;
    int[][] codes;
    int[][] codeLength;

    static double t1;

    public static class Node {
        double weight;
        public int value;
        public int leftIndex;
        public int rightIndex;

        public Node(double weight, int value, int leftIndex, int rightIndex) {
            this.weight = weight;
            this.value = value;
            this.leftIndex = leftIndex;
            this.rightIndex = rightIndex;
        }


        public String toString() {
            return weight + ", " + value + ", " + leftIndex + ", " + rightIndex;
        }
    }


    public int getMinimumWeightNode(int electrode) {
        double minWeight = Double.POSITIVE_INFINITY;
        int minIndex = -1;
        double w;
        for (int i = 0; i < nodesSize[electrode]; i++) {
            w = nodes[electrode][i].weight;
            if (w >= 0 && w < minWeight) {
                minWeight = w;
                minIndex = i;
            }
        }
        return minIndex;
    }


    public void getCodes(int electrode, int nodeIndex, int iCode, int nBits) {
        if (nodes[electrode][nodeIndex].leftIndex == -1) {
            // I = 2048 + V
            // V = I - 2048
            codes[electrode][nodes[electrode][nodeIndex].value - min] = iCode;
            codeLength[electrode][nodes[electrode][nodeIndex].value - min] = nBits;
        } else {
            nBits++;
            iCode <<= 1;
            getCodes(electrode, nodes[electrode][nodeIndex].leftIndex, iCode, nBits);
            getCodes(electrode, nodes[electrode][nodeIndex].rightIndex, iCode | 1, nBits);
        }
    }


    public void processSample(short[] sample) {
        // calculate progress
        sampleIndex++;
        percentage = 100 * sampleIndex / nSamples;
        if (percentage != oldPercentage) {
            System.out.println(percentage + "%");
            oldPercentage = percentage;
        }
    }


    public void finishSampleProcessing() throws IOException {
        out.close();
        sync.done();
    }


    public void compress(String inputFileName, String outputFileName) throws IOException {
        // open the file for random access, to make the probability distribution
        RawDataFile rawDataFile = new RawDataFile(new File(inputFileName));
        RawDataHeader header512 = rawDataFile.getHeader();
        nElectrodes = header512.getNumberOfElectrodes();
        nSamples = header512.getNumberOfSamples();

        nodes = new Node[nElectrodes][ (max - min) * 2 - 1];
        nodesSize = new int[nElectrodes];
        rootNodeIndex = new int[nElectrodes];
        codes = new int[nElectrodes][max - min];
        codeLength = new int[nElectrodes][max - min];
        short[] rawData = new short[1 * 20000];
        DoubleHistogram h = new DoubleHistogram("", min, max, 1);

        t1 = System.currentTimeMillis();

        for (int electrode = 0; electrode < nElectrodes; electrode++) {
            rawDataFile.getData(electrode, 20000, rawData);

            h.clear();
            for (int i = 0; i < rawData.length; i++) {
                h.fill(rawData[i], 1);
            }

//        PlotUtilities.showData("", h, new HistogramStyle());
            // build the Huffman tree
            for (int i = 0; i < h.getBinCount(); i++) {
                nodes[electrode][i] = new Node(h.getBin(i), (int) h.getMin() + i, -1, -1);
                nodesSize[electrode]++;
            }

            //while (getRemainingNodesCount(electrode) > 1) {
            for (int nn = nodesSize[electrode]; nn > 1; nn--) {
                int n1 = getMinimumWeightNode(electrode);
                double weight1 = nodes[electrode][n1].weight;
                nodes[electrode][n1].weight = -1;
                int n2 = getMinimumWeightNode(electrode);
                double weight2 = nodes[electrode][n2].weight;
                nodes[electrode][n2].weight = -1;

                nodes[electrode][nodesSize[electrode]] =
                    new Node(weight1 + weight2, -1, n1, n2);
                rootNodeIndex[electrode] = nodesSize[electrode];
                nodesSize[electrode]++;
            }

            // compute the codes
            getCodes(electrode, rootNodeIndex[electrode], 0, 0);

            // get the min/max code length
            int minCodeLength = Integer.MAX_VALUE;
            int maxCodeLength = Integer.MIN_VALUE;
            for (int i = 0; i < codeLength[electrode].length; i++) {
                if (codeLength[electrode][i] < minCodeLength) {
                    minCodeLength = codeLength[electrode][i];
                }
                if (codeLength[electrode][i] > maxCodeLength) {
                    maxCodeLength = codeLength[electrode][i];
                }
            }
            System.out.println(
                electrode + " codeLength " + minCodeLength + " - " + maxCodeLength);
        }

        System.out.println( (System.currentTimeMillis() - t1) / 1000.0 + " sec");
//        System.exit(1);

//                int size = 0;
//                h.clear();
//                for (int i = 0; i < rawData.length; i++) {
//                    h.fill(rawData[i], 1);
//                }
//                for (int i = 0; i < h.getBinCount(); i++) {
//                    size += h.getBin(i) * codeLength[i];
//                }
//                System.out.println(size / 8.0 / 1024.0 / 1024.0 + " Mb");
//                System.out.println(rawData.length * 1.5 / 1024.0 / 1024.0 + " Mb");
//                System.out.println(rawData.length * 1.5 / (size / 8.0) + " ratio");
//
//                // calculate Shannon codelength
//                h.scale(1 / h.getBinSum());
//                double H = 0;
//                for (int i = 0; i < h.getBinCount(); i++) {
//                    double pi = h.getBin(i);
//                    if (pi != 0) {
//                        H += -pi * Math.log(pi) / Math.log(2);
//                    }
//                }
//                System.out.println(H);

        t1 = System.currentTimeMillis();

        // open the input file for sample streaming
        MultipleCompressedSampleInputStream sis =
            new MultipleCompressedSampleInputStream(inputFileName);
        sis.addSampleListener(this);

        // open and prepare the output file
        DataOutputStream outputStream = new DataOutputStream(new FileOutputStream(
            outputFileName));
//        RawDataHeader512 newHeader = new RawDataHeader512(
//            0,
//            nElectrodes,
//            header512.getSamplingFrequency(),
//            header512.getNumberOfSamples(),
//            header512.getArrayID(),
//            RawDataHeader.FORMAT_HUFFMAN_COMPRESSED,
//            header512.getComment()
//            );
//        outputStream.write(newHeader.getBinaryRepresentation());

        // write the Huffman trees to the output file
        for (int electrode = 0; electrode < nElectrodes; electrode++) {
            outputStream.writeShort(rootNodeIndex[electrode]);

//            System.out.println(rootNodeIndex[electrode] + " root");

            for (int i = 0; i < nodes[electrode].length; i++) {
                outputStream.writeShort(nodes[electrode][i].leftIndex);
                outputStream.writeShort(nodes[electrode][i].rightIndex);

//                System.out.println(nodes[electrode][i].leftIndex + " left");
//                System.out.println(nodes[electrode][i].leftIndex + " right");
            }
        }
        outputStream.close();

        out = new BitOutputStream(new BufferedOutputStream(
            new FileOutputStream(outputFileName, true), 10 * 1024 * 1024));
//        out = new BitOutputStream(new FileOutputStream(outputFileName, true));

        sync = new SynchronizationObject();
        sync.setWorking();
        sis.start();
        sync.waitUntilDone();
        System.out.println("now exit the method");
    }


    public static void main(String[] args) throws IOException {
//        new HuffmanCoding().compress("d:\\data000-1.bin", "1.bin");
        new HuffmanCoding().compress("F:\\2005-05-04-0\\data007", "d:\\qa.bin");

        double t2 = System.currentTimeMillis();
        System.out.println("t: " + (t2 - t1) / 1000);
    }

}
