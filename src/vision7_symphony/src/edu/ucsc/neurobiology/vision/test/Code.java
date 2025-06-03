package edu.ucsc.neurobiology.vision.test;

import java.io.*;

import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.math.*;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class Code {
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


    public static int getMinimumWeightNode(Node[] nodes, int nodesSize) {
        double minWeight = Double.POSITIVE_INFINITY;
        int minIndex = -1;
        double w;
        for (int i = 0; i < nodesSize; i++) {
            w = nodes[i].weight;
            if (w >= 0 && w < minWeight) {
                minWeight = w;
                minIndex = i;
            }
        }
        return minIndex;
    }


    public static void getCodes(
        Node[] nodes, int nodeIndex, int[] codeLength, int iCode, int nBits) {

        if (nodes[nodeIndex].leftIndex == -1) {
            // I = 2048 + V
            // V = I - 2048
//            codes[nodes[nodeIndex].value] = iCode;
            codeLength[nodes[nodeIndex].value] = nBits;
        } else {
            nBits++;
            iCode <<= 1;
            getCodes(nodes, nodes[nodeIndex].leftIndex, codeLength, iCode, nBits);
            getCodes(nodes, nodes[nodeIndex].rightIndex, codeLength, iCode | 1, nBits);
        }
    }


    public static int[] calculateCodeLength(int[] ri2pi) {
        // build the Huffman tree
        Node[] nodes = new Node[ri2pi.length * 2 - 1];
        int nodesSize = 0;
        int rootNodeIndex = -1;
        for (int ri = 0; ri < ri2pi.length; ri++) {
            nodes[ri] = new Node(ri2pi[ri], ri, -1, -1);
            nodesSize++;
        }

        for (int nn = nodesSize; nn > 1; nn--) {
            int n1 = getMinimumWeightNode(nodes, nodesSize);
            double weight1 = nodes[n1].weight;
            nodes[n1].weight = -1;
            int n2 = getMinimumWeightNode(nodes, nodesSize);
            double weight2 = nodes[n2].weight;
            nodes[n2].weight = -1;

            nodes[nodesSize] =
                new Node(weight1 + weight2, -1, n1, n2);
            rootNodeIndex = nodesSize;
            nodesSize++;

            System.out.println("joining " + weight1 + ", " + weight2);
        }

        // compute the codes
        int[] codeLength = new int[ri2pi.length];
        getCodes(nodes, rootNodeIndex, codeLength, 0, 0);

        return codeLength;
    }


    public static void main(String[] args) throws IOException {
//        String[] si = {"bat", "cat", "eat", "fat", "hat", "mat", "oal", "pat", "rat",
//            "sat", "wat"};
        int[] si = new int[] {2, 14, 3, 5, 41, 6, 18, 9, 10, 1, 11};
        int[] pi = new int[] {8, 21, 8, 9, 23, 3, 10, 7, 21, 5, 6};

        // convert to an 1..n alphabet
        int[] si2ri = new int[50];
        int[] ri2si = new int[si.length];
        int[] ri2pi = new int[si.length];

        for (int ri = 0; ri < si.length; ri++) {
            int minSi = MathUtil.minIndex(pi);
            si2ri[si[minSi]] = ri;
            ri2si[ri] = si[minSi];
            ri2pi[ri] = pi[minSi];

            pi[minSi] = Integer.MAX_VALUE;
//            ri++;
        }

        IOUtil.printArray(ri2si);
        IOUtil.printArray(ri2pi);

        // calculating codeword length
        int[] codeLength = calculateCodeLength(ri2pi);
        // get the min/max code length
        int minCodeLength = Integer.MAX_VALUE;
        int maxCodeLength = Integer.MIN_VALUE;
        double K = 0;
        for (int i = 0; i < codeLength.length; i++) {
            if (codeLength[i] < minCodeLength) {
                minCodeLength = codeLength[i];
            }
            if (codeLength[i] > maxCodeLength) {
                maxCodeLength = codeLength[i];
            }
            K += Math.pow(2, -codeLength[i]);
        }

        // calculate the canonical codes
        int[] codes = new int[codeLength.length];
        for (int i = 0; i < codeLength.length; i++) {
            for (int j = 0; j <= i - 1; j++) {
                codes[i] += Math.pow(2, maxCodeLength - codeLength[j]);
            }
            codes[i] /= Math.pow(2, maxCodeLength - codeLength[i]);
        }

        System.out.println("codeLength = " + minCodeLength + " .. " + maxCodeLength);
        System.out.println("Kraft's number = " + K);
        for (int i = minCodeLength; i <= maxCodeLength; i++) {
            System.out.println("code " + i + ": " + codes[i]);
        }

        IOUtil.printArray(codeLength);
        IOUtil.printArray(codes);

    }

}
