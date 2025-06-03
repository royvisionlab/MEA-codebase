package edu.ucsc.neurobiology.vision.io;

import java.io.*;
import java.util.*;
import edu.ucsc.neurobiology.vision.util.VisionParams;


/**
 * This class provides random access to a raw data file or raw data directory.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 * @owner Matthew Grivich, The Salk Institute
 */
public class RawDataFile {      
    private RawDataHeader512 header;
    private RandomAccessFile[] fileList;
    private long[] fileSizes;
    private int nElectrodes;
    private int nFiles;
    File file;


    /**
     * Creates a Raw Data File.
     * If <code>file</code> represents a file then the instance of this class will read
     * just from that file.
     * If <code>file</code> represents a directory then the instance of this class
     * will read sequentially from all the raw data files in the directory. In this case
     * the whole directory will appear as a big raw data file.
     */
    public RawDataFile(String fileName) throws IOException {
        openRawDataFile(new File(fileName));
    }

    public RawDataFile(File file) throws IOException {
        openRawDataFile(file);
    }

    public void openRawDataFile(File file) throws IOException {
        this.file = file;

        //        if (file.length() == 0) {
        //            throw new IOException("Sorry, the file " + file.getName() + " is empty.");
        //        }

        if (file.isFile()) {
            // the list will contain only one item
            nFiles = 1;
            fileList = new RandomAccessFile[nFiles];
            fileList[0] = new RandomAccessFile(file, "r");
            this.header = new RawDataHeader512(new FileInputStream(file));
            
            // Patching for correct length of file during processing of
            // interrupted streaming
            // Vincent Deo - Stanford University - 10/15/2015
            header.setNumberOfSamples((int) ((fileList[0].length() - header.getHeaderSize()) / header.getSampleSize()));
                                
        } else if (file.isDirectory()) {
            // Use the content of the directory
            String[] list = file.list();

            // Filter out non-standard files (generally hidden files caused by OS X)
            ArrayList<String> listArrayList = new ArrayList<String>();
            for (int i = 0; i < list.length; i++) {
                if (list[i].endsWith(VisionParams.BIN_FILE_EXTENSION_512)) {
                    listArrayList.add(list[i]);
                }
            }
            list = new String[listArrayList.size()];
            for (int i=0; i<list.length; i++) {
                list[i] = listArrayList.get(i);
            }

            Arrays.sort(list);
            nFiles = list.length;
            fileList = new RandomAccessFile[nFiles];
            String fileSeparator = System.getProperty("file.separator");
            int nSamples = 0;
            
            if (nFiles == 0) {
                throw new IOException("Folder " + file.getAbsolutePath() + " has no raw data files in it.");
            }
            
            for (int i = 0; i < nFiles; i++) {
                String name = file.getAbsolutePath() + fileSeparator + list[i];
                if (i == 0) {
                    this.header = new RawDataHeader512(new FileInputStream(name));
                }

                fileList[i] = new RandomAccessFile(name, "r");
                if (i == 0) {
                    nSamples += (fileList[i].length() - header.getHeaderSize()) /
                        header.getSampleSize();
                } else {
                    nSamples += fileList[i].length() / header.getSampleSize();
                }
            }
            header.setNumberOfSamples(nSamples);
        } else {
            throw new IOException(
                                  file.getName() + " - is not a file nor a directory");
        }
        
        // get the file length
        fileSizes = new long[nFiles];
        for (int i = 0; i < nFiles; i++) {
            fileSizes[i] = fileList[i].length();
        }
        
        // initialize parameters
        this.nElectrodes = header.getNumberOfElectrodes();
    }


    public String getName() {
        return file.getAbsolutePath();
    }


    /**
     * Returns the raw data header information.
     * In the case of a single file the header is the one of that file.
     * In the case of a directory the header is read from the first file in the directory.
     */
    public RawDataHeader512 getHeader() {
        return header;
    }


    /**
     * Reads the data from all electrodes starting at sample <code>startSample</code>
     * until the array <code>data</code> is completely filled. If there is not enough
     * data to fill the array a <code>IllegalStateException</code> exception gets thrown.
     * <br>The method is thread safe.
     */
    synchronized public void getData(final long startSample, short[][] data) 
        throws IOException {

        //        throw new Error("not implemented");

        int sampleSize;
        if (nElectrodes % 2 == 1) {
            sampleSize = 2 + (nElectrodes - 1) * 3 / 2;
        } else {
            sampleSize = (nElectrodes) * 3 / 2;
        }
        byte[] b = new byte[sampleSize];
        int b1, b2, b3;
        long s = header.getHeaderSize() + startSample * sampleSize;
        int fileIndex = 0;
        for (fileIndex = 0; fileIndex < nFiles; fileIndex++) {
            if (s < fileSizes[fileIndex]) {
                break;
            } else {
                s -= fileSizes[fileIndex];
            }
        }
        fileList[fileIndex].seek(s);
        int i, electrode;
        for (int sampleIndex = 0; sampleIndex < data.length; sampleIndex++) {
            while (fileList[fileIndex].read(b) == -1) {
                fileIndex++;
                if (fileIndex == nFiles) {
                    throw new IllegalStateException("end of data: " + startSample);
                }
            }

            if (nElectrodes % 2 == 1) {
                data[sampleIndex][0] = (short) ( (b[0] << 8) + (b[1] & 0xff)); //TTL
                for (i = 2, electrode = 1; i < sampleSize; ) {
                    b1 = b[i++] & 0xff;
                    b2 = b[i++] & 0xff;
                    b3 = b[i++] & 0xff;
                    data[sampleIndex][electrode++] =
                        (short) ( (b1 << 4) + (b2 >> 4) - 2048);
                    data[sampleIndex][electrode++] =
                        (short) ( ( (b2 & 0x0f) << 8) + b3 - 2048);
                }
            } else {
                for (i = 0, electrode = 0; i < sampleSize; ) {
                    b1 = b[i++] & 0xff;
                    b2 = b[i++] & 0xff;
                    b3 = b[i++] & 0xff;
                    data[sampleIndex][electrode++] =
                        (short) ( (b1 << 4) + (b2 >> 4) - 2048);
                    data[sampleIndex][electrode++] =
                        (short) ( ( (b2 & 0x0f) << 8) + b3 - 2048);
                }
            }
        }

    }


    /**
     * Same as getData(final long startSample, short[][] data),
     * except that the data array is initialized locally.
     * Written for calling Java from Matlab compatibility.
     */
    synchronized public short[][] getData(final long startSample, final int samples) 
        throws IOException {
        
        short[][] data = new short[samples][nElectrodes];
        getData(startSample, data);
        return data;
    }

    /**
     * Reads the data from the given <code>electrodes</code> starting at sample
     * <code>startSample</code> until the array <code>data</code> is completely filled.
     * If there is not enough data to fill the array a <code>IllegalStateException</code>
     * exception gets thrown.
     * <br>The method is thread safe.
     */
    synchronized public void getData(final int electrode, final long startSample, short[] data) 
        throws IOException {

        //        System.out.println("startSample = " + startSample);
        //        System.out.println("data.lenght = " + data.length);

        byte[] b;

        //      switch (header.getFormat()) {
        //         case RawDataHeader.FORMAT_12BIT_COMPRESSED:
        final int sampleSize;
        final int ind;

        int blockIndex = (electrode - 1) / 2;
        if (nElectrodes % 2 == 1) {
            sampleSize = 2 + (nElectrodes - 1) * 3 / 2;
            ind = 2 + 3 * blockIndex;
        } else {
            sampleSize = (nElectrodes) * 3 / 2;
            ind = 3 * blockIndex;
        }

        long s = header.getHeaderSize() + startSample * (long) sampleSize;
        int fileIndex = 0;
        for (fileIndex = 0; fileIndex < nFiles; fileIndex++) {
            if (s > fileSizes[fileIndex]) {
                s -= fileSizes[fileIndex];
            } else {
                break;
            }
        }
        fileList[fileIndex].seek(s);

        b = new byte[sampleSize];
        int b1, b2, b3;

        for (int n = 0; n < data.length; ) {
            while (fileList[fileIndex].read(b) == -1) {
                fileIndex++;
                if (fileIndex == nFiles) {
                    throw new IllegalStateException("end of data");
                }
            }

            if (nElectrodes % 2 == 1) {
                if (electrode == 0) {
                    data[n++] = (short) ( (b[0] << 8) + (b[1] & 0xFF));
                } else {
                    b1 = b[ind + 0] & 0xff;
                    b2 = b[ind + 1] & 0xff;
                    b3 = b[ind + 2] & 0xff;

                    if (electrode % 2 == 1) {
                        // odd electrode number
                        data[n++] = (short) ( (b1 << 4) + (b2 >> 4) - 2048);
                    } else {
                        // even electrode number
                        data[n++] = (short) ( ( (b2 & 0x0f) << 8) + b3 - 2048);
                    }
                }
            } else {
                //                        System.err.println("--");

                b1 = b[ind + 0] & 0xff;
                b2 = b[ind + 1] & 0xff;
                b3 = b[ind + 2] & 0xff;

                if (electrode % 2 == 1) {
                    // odd electrode number
                    data[n++] = (short) ( (b1 << 4) + (b2 >> 4) - 2048);
                } else {
                    // even electrode number
                    data[n++] = (short) ( ( (b2 & 0x0f) << 8) + b3 - 2048);
                }
            }
        }
        //             break;

        //            case RawDataHeader.FORMAT_8BIT_UNCOMPRESSED:
        //                System.out.println("FORMAT_8BIT_UNCOMPRESSED");
        //
        //                fileList[0].seek(
        //                    (long) header.getHeaderSize() +
        //                    3L * (electrode * header.getNumberOfSamples() + startSample) / 2L);
        //                b = new byte[ (int) (3L * data.length / 2L)];
        //                fileList[0].readFully(b);
        //
        //                for (int i = 0, byteIndex = 0; i < data.length; ) {
        ////                    data[i] = (short) ( (b[byteIndex++] << 8) + (b[byteIndex++] & 0xff));
        //
        //                    b1 = b[byteIndex++] & 0xff;
        //                    b2 = b[byteIndex++] & 0xff;
        //                    b3 = b[byteIndex++] & 0xff;
        //
        //                    data[i++] = (short) ( (b1 << 4) + (b2 >> 4) - 2048);
        //                    data[i++] = (short) ( ( (b2 & 0x0f) << 8) + b3 - 2048);
        //                }
        //                break;

        //          default:

        //                System.out.println("DEFAULT");
        //               break;
        //      }
    }


    /**
     * Same as getData(final int electrode, final long startSample, short[] data),
     * except that the data array is initialized locally.
     * Written for calling Java from Matlab compatibility.
     */
    synchronized public short[] getData(final int electrode,
                                        final long startSample,
                                        int samples) throws IOException {
        short[] data = new short[samples];
        getData(electrode, startSample, data);
        return data;
    }


    /**
     * This should be designed so that it is safe to run it more than once.  RandomAccessFile#close() implements
     * Closeable, which is specified to have no effect if run more than once.
     * 
     * @throws IOException
     */
    public void close() throws IOException {
        for (int i = 0; i < fileList.length; i++) {
            fileList[i].close();
        }
    }

    @Override
    /**
     * This is a backup safeguard, mostly here because people seemed to forget to close STAFiles and 
     * EI files appropriately from Matlab.  Do not rely on this to close files; you should always run close
     * explicitly.
     */
    protected void finalize() throws Throwable {
        try {
            close();
        } catch (Exception e) {} finally {
            super.finalize();
        }
    }
    
    @Override
    public String toString() {
        return getHeader().toString();
    }
}
