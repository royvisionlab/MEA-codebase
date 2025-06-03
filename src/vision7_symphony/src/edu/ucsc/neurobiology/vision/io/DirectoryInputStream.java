package edu.ucsc.neurobiology.vision.io;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.util.*;


/**
 * This class represents an <code>InputStream</code> implementation which gets the
 * input by reading the files from a given directory. All the files for which the
 * name that start with the same name as the directory name are read in alphabetical
 * order.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public final class DirectoryInputStream
    extends InputStream {
    LinkedList<InputStream> streams = new LinkedList<InputStream>();
    private InputStream input;


    /**
     * Creates a <code>DirectroryInputStreeam</code> which will read the files from the
     * given directory.
     */
    public DirectoryInputStream(String directoryName) throws IOException {
        File dir = new File(directoryName);
        if (!dir.isDirectory()) {
            throw new IllegalArgumentException(directoryName + " - not a directory");
        }

        String[] fileName = dir.list();
        if (fileName.length == 0) {
            throw new IllegalArgumentException(directoryName + " - is empty");
        }

        Arrays.sort(fileName);
        for (int i = 0; i < fileName.length; i++) {
            String fullFilePath = directoryName + File.separator + fileName[i];
            if (fullFilePath.endsWith(VisionParams.BIN_FILE_EXTENSION_512))
                streams.add(new FileInputStream(fullFilePath));
        }

        input = (InputStream) streams.remove(0);
    }


    static public long getDirectorySize(String directoryName) {
        File dir = new File(directoryName);
        if (!dir.isDirectory()) {
            throw new IllegalArgumentException(directoryName + " - not a directory");
        }

        String[] fileName = dir.list();
        if (fileName.length == 0) {
            throw new IllegalArgumentException(directoryName + " - is empty");
        }

        Arrays.sort(fileName);
        long length = 0;
        for (int i = 0; i < fileName.length; i++) {    	
            if (fileName[i].endsWith(VisionParams.BIN_FILE_EXTENSION_512)) {
                File f = new File(directoryName + File.separator + fileName[i]);
                length += f.length();
            }
        }
        return length;
    }


    /**
     * Reads <code>length</code> bytes and puts them into <code>b</code> starting at the
     * given <code>offset</code>. If the end of one file is reached the next one will be
     * read. If the end of all files are reached -1 gets returned, otherwise the method
     * returns the number of bytes actually read.
     */
    public int read(byte[] b, int offset, int length) throws IOException {
        int totalBytesRead = 0;

        while (length != 0) {
            int bytesRead = input.read(b, offset, length);

            if (bytesRead == length) {
                return totalBytesRead + length;
            } else {
                if (bytesRead == -1) {
                    if (streams.size() == 0) {
                        return (totalBytesRead == 0) ? -1 : totalBytesRead;
                    } else {
                        input = (InputStream) streams.remove(0);
                    }
                } else {
                    offset += bytesRead;
                    length -= bytesRead;
                    totalBytesRead += bytesRead;
                }
            }
        }

        return totalBytesRead;
    }


    /**
     * Reads a single byte from the underlying data source.
     *
     * @return the byte or -1 if no more data
     */
    public int read() throws IOException {
        while (true) {
            int b = input.read();

            if (b != -1) {
                return b;
            } else {
                if (streams.size() == 0) {
                    return -1;
                } else {
                    input = (InputStream) streams.remove(0);
                }
            }
        }
    }


    public int available() throws IOException {
        return input.available();
    }


    public void close() throws IOException {
        input.close();
    }


    public void mark(int readlimit) {
        throw new IllegalAccessError("Method not supported");
    }


    public synchronized boolean markSupported() {
        return false;
    }


    /**
     * Reads until the array <code>b</code> is filled.
     * @see #read(byte[], int, int)
     */
    public int read(byte[] b) throws IOException {
        return read(b, 0, b.length);
    }


    public void reset() throws IOException {
        throw new IOException("Method not supported");
    }

//NEED TO IMPLEMENT
    public long skip(long n) throws IOException {
        throw new IllegalAccessError("Method not supported");
    }

}
