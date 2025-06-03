package edu.ucsc.neurobiology.vision.io;

import java.io.*;

import edu.ucsc.neurobiology.vision.anf.*;
import edu.ucsc.neurobiology.vision.electrodemap.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ProjectionsFile_old {
    private RandomAccessFile projectionsFile;

    public final static int HEADER_SIZE = 7 * 5;
    public final static int ENTRY_SIZE = 8 + 4;

    public final int writeBufferSize = 5000;

    public final int arrayID;
    public final int nSamples;
    public final int samplingFrequency;
    public final ElectrodeUsage electrodeUsage;
    public final int nlPoints;
    public final int nrPoints;
    public final int nDimensions;
    public final int maxSpikesPerElectrode;
    public final int nElectrodes;

    private byte[] readBuffer;
    private byte[] writeBuffer;
    private long[] startLocation;
    private int[] writtenSpikes;


    public ProjectionsFile_old(String prjFileName) throws IOException {
        projectionsFile = new RandomAccessFile(prjFileName, "r");

        this.arrayID = projectionsFile.readInt();
        this.nSamples = projectionsFile.readInt();
        this.samplingFrequency = projectionsFile.readInt();
        int electrodeUsageInt = projectionsFile.readInt();
        this.nlPoints = projectionsFile.readInt();
        this.nrPoints = projectionsFile.readInt();
        this.nDimensions = projectionsFile.readInt();

        if (electrodeUsageInt < 0 || electrodeUsageInt >= ElectrodeUsage.values().length) {
            throw new IOException("Corrupt projections file: " + this.toString());
        }
        this.electrodeUsage = ElectrodeUsage.values()[electrodeUsageInt];

        if (!areValidParams()) {
            throw new IOException("Corrupt projections file: " + this.toString());
        }

        this.nElectrodes = ElectrodeMapFactory.getElectrodeMap(arrayID).
                           getNumberOfElectrodes();

        int maxSpikesPerElectrode = Integer.MIN_VALUE;
        for (int electrode = 0; electrode < nElectrodes; electrode++) {
            projectionsFile.seek(HEADER_SIZE + electrode * ENTRY_SIZE);
            long seekLocation = projectionsFile.readLong(); // do not remove
            int nSpikes = projectionsFile.readInt();

            if (nSpikes > maxSpikesPerElectrode) {
                maxSpikesPerElectrode = nSpikes;
            }
        }
        this.maxSpikesPerElectrode = maxSpikesPerElectrode;

        readBuffer = new byte[maxSpikesPerElectrode * (nDimensions + 1) * 4];
    }


    public ProjectionsFile_old(
        String prjFileName, int arrayID, int nSamples,
        int samplingFrequency, ElectrodeUsage electrodeUsage, int nlPoints, int nrPoints,
        int nDimensions, int[] spikeCount) throws IOException {

        if (!areValidParams()) {
            throw new IOException("Corrupt projections file: " + this.toString());
        }

        this.arrayID = arrayID;
        this.nSamples = nSamples;
        this.samplingFrequency = samplingFrequency;
        this.electrodeUsage = electrodeUsage;
        this.nlPoints = nlPoints;
        this.nrPoints = nrPoints;
        this.nDimensions = nDimensions;

        projectionsFile = new RandomAccessFile(prjFileName, "rw");
        projectionsFile.writeInt(arrayID);
        projectionsFile.writeInt(nSamples);
        projectionsFile.writeInt(samplingFrequency);
        projectionsFile.writeInt(electrodeUsage.ordinal());
        projectionsFile.writeInt(nlPoints);
        projectionsFile.writeInt(nrPoints);
        projectionsFile.writeInt(nDimensions);

        this.nElectrodes = ElectrodeMapFactory.getElectrodeMap(arrayID).
                           getNumberOfElectrodes();

        startLocation = new long[nElectrodes];
        long index = HEADER_SIZE + ENTRY_SIZE * nElectrodes;
        for (int i = 0; i < nElectrodes; i++) {
            if (spikeCount[i] != 0) {
                startLocation[i] = index;
                index += spikeCount[i] * (nDimensions + 1) * 4;
            } else {
                startLocation[i] = -1;
            }
        }

        writtenSpikes = new int[nElectrodes];

        this.maxSpikesPerElectrode = -1;

        writeBuffer = new byte[writeBufferSize * (nDimensions + 1) * 4];
    }


    public int[] getTTLTimes() {
        throw new Error("NOT IMPLEMENTED");
    }


    private boolean areValidParams() {
        if (arrayID < 0 || nSamples < 0 || samplingFrequency < 0 || nlPoints < 0 ||
            nrPoints < 0 || nDimensions < 0) {
            return false;
        } else {
            return true;
        }
    }


    synchronized public int readProjections(int electrode, float[][] data, int[] times) throws
        IOException {

        projectionsFile.seek(HEADER_SIZE + electrode * ENTRY_SIZE);
        long seekLocation = projectionsFile.readLong();
        if (seekLocation == -1) {
            return -1;
        }
        int nSpikes = projectionsFile.readInt();

        projectionsFile.seek(seekLocation);
        projectionsFile.readFully(readBuffer, 0, nSpikes * (nDimensions + 1) * 4);

        int byteIndex = 0, ch1, ch2, ch3, ch4;
        for (int i = 0; i < nSpikes; i++) {
            // reconstruct the spike time
            ch1 = readBuffer[byteIndex++] & 0xff;
            ch2 = readBuffer[byteIndex++] & 0xff;
            ch3 = readBuffer[byteIndex++] & 0xff;
            ch4 = readBuffer[byteIndex++] & 0xff;
            if ( (ch1 | ch2 | ch3 | ch4) < 0) {
                throw new EOFException();
            }
            times[i] = (ch1 << 24) + (ch2 << 16) + (ch3 << 8) + (ch4 << 0);

            // reconstruct the projections
            for (int d = 0; d < nDimensions; d++) {
                ch1 = readBuffer[byteIndex++] & 0xff;
                ch2 = readBuffer[byteIndex++] & 0xff;
                ch3 = readBuffer[byteIndex++] & 0xff;
                ch4 = readBuffer[byteIndex++] & 0xff;
                if ( (ch1 | ch2 | ch3 | ch4) < 0) {
                    throw new EOFException();
                }
                data[d][i] = Float.intBitsToFloat(
                    (ch1 << 24) + (ch2 << 16) + (ch3 << 8) + (ch4 << 0));
            }
        }

        return nSpikes;
    }


    synchronized public void saveData(
        final int electrode, int nSpikesToSave, int[][] spikeTimes,
        float[][][] spikeProjections) throws IOException {

        // prepare the data to be written
        int byteIndex = 0, intBits;
        for (int i = 0; i < nSpikesToSave; i++) {
            // convert to bytes the spike time
            intBits = spikeTimes[electrode][i];
            writeBuffer[byteIndex++] = (byte) ( (intBits >>> 24) & 0xFF);
            writeBuffer[byteIndex++] = (byte) ( (intBits >>> 16) & 0xFF);
            writeBuffer[byteIndex++] = (byte) ( (intBits >>> 8) & 0xFF);
            writeBuffer[byteIndex++] = (byte) ( (intBits >>> 0) & 0xFF);

            // convert to bytes the projections
            for (int d = 0; d < nDimensions; d++) {
                intBits = Float.floatToIntBits(spikeProjections[electrode][i][d]);
                writeBuffer[byteIndex++] = (byte) ( (intBits >>> 24) & 0xFF);
                writeBuffer[byteIndex++] = (byte) ( (intBits >>> 16) & 0xFF);
                writeBuffer[byteIndex++] = (byte) ( (intBits >>> 8) & 0xFF);
                writeBuffer[byteIndex++] = (byte) ( (intBits >>> 0) & 0xFF);
            }
        }

        // write the data
        projectionsFile.seek(startLocation[electrode] +
                             writtenSpikes[electrode] * (nDimensions + 1) * 4);
        projectionsFile.write(writeBuffer, 0, byteIndex);
        writtenSpikes[electrode] += nSpikesToSave;
    }


    public void close() throws IOException {
        // write the header
        projectionsFile.seek(HEADER_SIZE);
        for (int i = 0; i < nElectrodes; i++) {
            projectionsFile.writeLong(startLocation[i]);
            projectionsFile.writeInt(writtenSpikes[i]);
        }

        projectionsFile.close();
    }


    public String toString() {
        return
//            "\n arrayID = " + arrayID +
//            "\n nSamples = " + nSamples +
//            "\n samplingFrequency = " + samplingFrequency +
//            "\n electrodeUsage = " + electrodeUsage +
//            "\n nlPoints = " + nlPoints +
//            "\n nrPoints = " + nrPoints +
//            "\n nDimensions = " + nDimensions;
            " arrayID = " + arrayID +
            " nSamples = " + nSamples +
            " samplingFrequency = " + samplingFrequency +
            " electrodeUsage = " + electrodeUsage +
            " nlPoints = " + nlPoints +
            " nrPoints = " + nrPoints +
            " nDimensions = " + nDimensions;
    }


    public static void main(String[] args) throws IOException {
        ProjectionsFile f = new ProjectionsFile("Y:\\2000-12-14-0\\data051\\data051.prj");
//        Thread
    }
}
