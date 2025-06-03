package edu.ucsc.neurobiology.vision.io;

import java.io.*;

import static edu.ucsc.neurobiology.vision.math.MathUtil.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class HistogramFile {
    public static final int BASIC_VERSION = 0;

    private final int headerLength = 4 + 8 + 8 + 4;
    private final int histogramLength;

    public final double minValue;
    public final double maxValue;
    public final int nBins;
    private RandomAccessFile raf;
    private IntegerList idList = new IntegerList();


    public HistogramFile(String fileName, double minValue, double maxValue, int nBins) throws
        IOException {

        if (maxValue < minValue) {
            throw new IllegalArgumentException("Wrong Histogram Limits");
        }

        if (nBins < 1) {
            throw new IllegalArgumentException("nBins should be at least 1");
        }

        this.minValue = minValue;
        this.maxValue = maxValue;
        this.nBins = nBins;

        raf = new RandomAccessFile(fileName, "rw");
        raf.setLength(0);
        raf.seek(0);

        raf.writeInt(BASIC_VERSION);
        raf.writeDouble(minValue);
        raf.writeDouble(maxValue);
        raf.writeInt(nBins);

        histogramLength = nBins * 8 * 2;
    }


    public HistogramFile(String fileName) throws
        IOException {

        raf = new RandomAccessFile(fileName, "rw");
        int version = raf.readInt();
        minValue = raf.readDouble();
        maxValue = raf.readDouble();
        nBins = raf.readInt();

        System.out.println("version " + version);
        System.out.println("minValue " + minValue);
        System.out.println("maxValue " + maxValue);
        System.out.println("nBins " + nBins);

        histogramLength = nBins * 8 * 2;

        int i = 0;
        while (headerLength + i * (histogramLength + 4) < raf.length()) {
            raf.seek(headerLength + i * (histogramLength + 4));
            int id = raf.readInt();
            idList.add(id);
            i++;
        }
    }


    synchronized public int[] getIDList() {
        return idList.toArray();
    }


    synchronized public void addHistogram(int id, DoubleErrorHistogram h) throws
        IOException {
        //        if (h.getMin() != minValue || h.getMax() != maxValue ||
        if (!areEqual(h.getMin(), minValue) ||
            !areEqual(h.getMax(), maxValue) ||
            !areEqual(h.getBinCount(), nBins)) {

            throw new IllegalArgumentException("Invalid Histogram");
        }

        if (idList.contains(id)) {
            throw new IllegalArgumentException("The file allready contains id " + id);
        }

        // go to the end of the file
        raf.seek(raf.length());
        raf.writeInt(id);
        for (int i = 0; i < nBins; i++) {
            raf.writeDouble(h.getBin(i));
            raf.writeDouble(h.getBinError(i));
        }
    }


    synchronized public DoubleErrorHistogram getHistogram(int id) throws IOException {
        if (!idList.contains(id)) {
            throw new IllegalArgumentException("The file does not contain id " + id);
        }

        int index = -1;
        for (int i = 0; i < idList.size(); i++) {
            if (idList.get(i) == id) {
                index = i;
                break;
            }
        }

        raf.seek(headerLength + index * (histogramLength + 4));
        int _id = raf.readInt();
        double[] value = new double[nBins];
        double[] error = new double[nBins];
        for (int i = 0; i < nBins; i++) {
            value[i] = raf.readDouble();
            error[i] = raf.readDouble();
        }

        DoubleHistogram valueHist = new DoubleHistogram("", value, minValue, maxValue);
        DoubleHistogram errorHist = new DoubleHistogram("", error, minValue, maxValue);

        return new DoubleErrorHistogram(valueHist, errorHist);
    }


    synchronized public void close() throws IOException {
        raf.close();
    }

}
