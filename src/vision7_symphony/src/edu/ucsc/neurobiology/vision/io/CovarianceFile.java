package edu.ucsc.neurobiology.vision.io;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.math.*;


/**
 * An interface to the covariance matrices calculated from the spike times and the raw data
 * for every electrode. The covariance matrices are saved in reduced form, that is,
 * only half the matrix + the diagonal. See the CovarianceMatrix code for details.
 *
 * @see CovarianceMatrix
 * @author Dumitru Petrusca, University of California, Santa Cruz
 * @owner Matthew Grivich, The Salk Institute
 */
public class CovarianceFile {
    private static final int MAGIC = 0xAAAAAA;

    private RandomAccessFile covFile;
    private HashMap<Integer, Long> seekLocations = new HashMap<Integer, Long>();
    private byte[] buffer;
    private String covFileName;

    private VisionHeader header;



    public CovarianceFile(String covFileName, VisionHeader header) throws IOException {
        this.covFileName = covFileName;
        header.magic = MAGIC;
        header.version = 1;

        this.header = header;

        // write the header
        covFile = new RandomAccessFile(covFileName, "rw");
        IOUtil.writePublicFields(header, covFile);
    }


    public CovarianceFile(String covFileName) throws IOException {
        this.covFileName = covFileName;
        covFile = new RandomAccessFile(covFileName, "rw");
        header = new VisionHeader();
        IOUtil.readPublicFields(header, covFile);

        if (header.magic != MAGIC) {
            throw new IOException("Requested covariance file is not of the proper type.");
        }

        long seekIndex = header.size();
        while (seekIndex < covFile.length()) {
            covFile.seek(seekIndex);
            int electrode = covFile.readInt();
            int matrixSize = covFile.readInt();

            seekLocations.put(electrode, seekIndex);

            seekIndex += 4 + 4 + 4 * matrixSize;
        }
    }


    public int getArrayID() {
        return header.arrayID;
    }


    public void addCovariance(CovarianceMatrix covMatrix, int electrode) throws
        IOException {

        if (header == null) {
            throw new NullPointerException("header is null");
        }

        if (covMatrix == null || covMatrix.getNPoints() < header.minCovarianceSpikes) {
            return;
        }

        covFile.writeInt(electrode);
        float[] c = covMatrix.getCovariance();
        covFile.writeInt(c.length);

        for (int i = 0; i < c.length; i++) {
            covFile.writeFloat(c[i]);
        }
    }


    private void ensureBufferSize(int size) {
        if (buffer == null || buffer.length < size) {
            buffer = new byte[size];
        }
    }


    synchronized public float[] getCovarianceMatrix(int electrode) throws IOException {
        if (!seekLocations.containsKey(electrode)) {
            return null;
        }

        long seekIndex = seekLocations.get(electrode);

        covFile.seek(seekIndex);
        int _electrode = covFile.readInt();
        if (_electrode != electrode) {
            throw new Error("Internal CovarinaceFile class error");
        }
        int matrixSize = covFile.readInt();

        ensureBufferSize(matrixSize * 4);
        covFile.readFully(buffer, 0, matrixSize * 4);
        return IOUtil.convertToFloatArray(buffer, matrixSize);
    }

/**
 * 
 * @param electrode
 * @param matrix
 * @throws IOException
 * 
 * Used by Tim Machado's code.
 */
    synchronized public void saveCovarianceMatrix(int electrode, float[] matrix) throws
        IOException {
        if (seekLocations.containsKey(electrode)) {
            throw new IOException("Covariance already exists for electrode " + electrode);
        }

        long seekIndex = covFile.length();
        seekLocations.put(electrode, seekIndex);

        covFile.seek(seekIndex);
        covFile.writeInt(electrode);
        covFile.writeInt(matrix.length);

        ensureBufferSize(matrix.length * 4);
        IOUtil.convertToByteArray(matrix, buffer);
        covFile.write(buffer, 0, matrix.length * 4);
    }
    
    public String getFileName() {
        return covFileName;
    }
    

    // ToDo: 2011-01-18: Shouldn't max start at -1?
    public int getMaxElectrode() {
        int max = 0;
        for (Integer key:seekLocations.keySet()) {
            int val = key.intValue();
            if (val > max) {
                max = val;
            }
        }
        return max;
    }


    public void close() throws IOException {
        covFile.close();
    }


    public VisionHeader getHeader() {
        return (VisionHeader) header.clone();
    }

}
