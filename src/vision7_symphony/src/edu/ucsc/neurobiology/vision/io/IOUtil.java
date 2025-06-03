package edu.ucsc.neurobiology.vision.io;

import java.io.*;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.net.*;
import java.text.*;
import java.util.*;
import java.util.regex.PatternSyntaxException;

import edu.ucsc.neurobiology.vision.io.tags.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * This class is intended to contain miscellaneous IO utility
 * methods used all over the program.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public final class IOUtil {

    private IOUtil() {}


    public static float[] convertToFloatArray(byte[] buffer, int nFloats) {
        int bIndex = 0;
        int ch1, ch2, ch3, ch4;
        float[] array = new float[nFloats];

        for (int i = 0; i < nFloats; i++) {
            ch1 = buffer[bIndex++] & 0xff;
            ch2 = buffer[bIndex++] & 0xff;
            ch3 = buffer[bIndex++] & 0xff;
            ch4 = buffer[bIndex++] & 0xff;
            array[i] = Float.intBitsToFloat(
                (ch1 << 24) + (ch2 << 16) + (ch3 << 8) + (ch4 << 0));
        }

        return array;
    }


    public static void convertToByteArray(float[] array, byte[] buffer) {
        for (int i = 0, bIndex = 0; i < array.length; i++) {
            int v = Float.floatToIntBits(array[i]);
            buffer[bIndex++] = (byte) ( (v >>> 24) & 0xFF);
            buffer[bIndex++] = (byte) ( (v >>> 16) & 0xFF);
            buffer[bIndex++] = (byte) ( (v >>> 8) & 0xFF);
            buffer[bIndex++] = (byte) ( (v >>> 0) & 0xFF);
        }
    }


    public static ArrayList<File> findInClasspath(String extension) {
        ArrayList<File> files = new ArrayList<File> ();

        String[] classPath = StringUtil.decomposeString(System.getProperty(
            "java.class.path"), ";");
        for (String path : classPath) {
            File f = new File(path);
            if (f.isFile() && f.canRead() && path.endsWith(extension)) {
                files.add(f);
            } else if (f.isDirectory()) {
                File[] list = f.listFiles();
                for (int i = 0; i < list.length; i++) {
                    if (list[i].isFile() && list[i].canRead() &&
                        list[i].getName().endsWith(extension)) {
                        files.add(list[i]);
                    }
                }
            }
        }

        return files;
    }


    public static double readNumber(StreamTokenizer st, String name) throws IOException {
        st.nextToken();
        if (st.ttype != st.TT_WORD || !st.sval.toUpperCase().equals(name)) {
            throw new IOException("Missing: " + name + " at line " + st.lineno());
        }
        st.nextToken();
        
        if (st.ttype != st.TT_NUMBER) {
            throw new IOException("Wrong type: number expected at line " + st.lineno());
        }
        return st.nval;
    }


    public static final int readFully(
        InputStream input, byte b[], int off, int len) throws IOException {

        int n = 0;
        do {
            int count = input.read(b, off + n, len - n);
            if (count < 0) {
                break;
            }
            n += count;
        } while (n < len);

        return n;
    }


    public static long getInputStreamSize(String rawDataSource) {
        if (rawDataSource.startsWith("net://")) {
            return -1;
        } else {
            File file = new File(rawDataSource);
            if (file.isFile()) {
                return file.length();
            } else if (file.isDirectory()) {
                return DirectoryInputStream.getDirectorySize(rawDataSource);
            } else {
                throw new IllegalArgumentException(
                    "Not file nor directory: " + rawDataSource);
            }
        }
    }


    /**
     * Calculates teh length of any Tag by writing it into a buffer and
     * looking for the buffers size.
     *
     * @param tag the tag to calculate the length of
     */
    public static int getTagLength(Tag tag) {
        try {
            ByteArrayOutputStream bos = new ByteArrayOutputStream();
            NeuroOutputStream nos = new NeuroOutputStream(bos);
            nos.writeTag(tag);
            nos.flush();
            nos.close();
            return bos.size();
        } catch (IOException e) {
            return -1;
        }
    }


    public static boolean isValidFile(File f) {
        if (f.exists() && f.isFile() && f.canRead()) {
            return true;
        } else {
            return false;
        }
    }


    public static boolean isValidFile(String name) {
        return isValidFile(new File(name));
    }


    public static boolean isValidFolder(File f) {
        if (f.exists() && f.isDirectory()) {
            return true;
        } else {
            return false;
        }
    }


    public static boolean isValidFolder(String name) {
        File f = new File(name);
        if (f.exists() && f.isDirectory()) {
            return true;
        } else {
            return false;
        }
    }
    
    public static boolean isValidDataset(String name) {
        DataFileStringParser parser = null;
        try {
         parser = new DataFileStringParser(name);
        } catch(PatternSyntaxException ex) {
            return false;
        }
 
        String[] datasets = parser.getDatasets();
        if(datasets.length == 0) {
            return false;
        }
        for(int i=0; i<datasets.length; i++) {
            File f = new File(datasets[i]);
            if(!(f.exists() && (f.isDirectory() || f.isFile()))) {
                return false;
            }
        }
//		//Check to see if it is a legit data file
//		boolean good = true;
//		MultipleCompressedSampleInputStream stream = null;
//		try {
//			stream = new MultipleCompressedSampleInputStream(name);
//
//		} catch(Exception e) {
//			good = false;
//		}
//		if(stream != null) {
//			stream.finishSampleProcessing();
//		}
//		return good;
        return true;
    }





    /**
     * Loads the content of the file with the given <i>fileName</i> as a double array.
     * The file should contain a single IEEE double on each row.
     * @return the array with the loaded numbers
     * @throws IOException if the file cannot be opened/read
     */
    public static double[] loadDoubleArray(String fileName) throws IOException {
        LineNumberReader lnr = null;
        try {
            lnr = new LineNumberReader(new FileReader(fileName));

            DoubleList list = new DoubleList();
            String line;
            while ( (line = lnr.readLine()) != null) {
                try {
                    list.add(Double.parseDouble(line));
                } catch (NumberFormatException e) {
                    throw new IOException("Unknown format on line: " + lnr.getLineNumber());
                }
            }

            return list.toArray();
        } finally {
            if (lnr != null) {
                lnr.close();
            }
        }
    }


    /**
     * Loads the content of the file with the given <i>fileName</i> as a double array.
     * The file should contain a single IEEE double on each row.
     * @return the array with the loaded numbers
     * @throws IOException if the file cannot be opened/read
     */
    public static float[] loadFloatArray(String fileName) throws IOException {
        LineNumberReader lnr = null;
        try {
            lnr = new LineNumberReader(new FileReader(fileName));

            FloatList list = new FloatList();
            String line;
            while ( (line = lnr.readLine()) != null) {
                try {
                    list.add(Float.parseFloat(line));
                } catch (NumberFormatException e) {
                    throw new IOException("Unknown format on line: " + lnr.getLineNumber());
                }
            }

            return list.toArray();
        } finally {
            if (lnr != null) {
                lnr.close();
            }
        }
    }


    /**
     * Loads the content of the file with the given <i>fileName</i> as an integer array.
     * The file should contain a single int on each row.
     * @return the array with the loaded numbers
     * @throws IOException if the file cannot be opened/read
     */
    public static int[] loadIntArray(String fileName) throws IOException {
        LineNumberReader lnr = null;
        try {
            lnr = new LineNumberReader(new FileReader(fileName));

            IntegerList list = new IntegerList();
            String line;
            while ( (line = lnr.readLine()) != null) {
                try {
                    list.add(Integer.parseInt(line));
                } catch (NumberFormatException e) {
                    throw new IOException("Unknown format on line: " + lnr.getLineNumber());
                }
            }

            return list.toArray();
        } finally {
            if (lnr != null) {
                lnr.close();
            }
        }
    }


    /**
     * Loads the content of the file with the given <i>fileName</i> as a
     * two-dimensional array of doubles (matrix). The file should contain a single row
     * of the matrix on each row in the file.
     *
     * @return the matrix
     * @throws IOException if the file cannot be opened/read
     */
    public static double[][] loadMatrix(String fileName) throws IOException {
        LineNumberReader lnr = null;
        try {
            lnr = new LineNumberReader(new FileReader(fileName));
            ArrayList<DoubleList> rows = new ArrayList<DoubleList>();
            String line;

            while ( (line = lnr.readLine()) != null) {
                DoubleList row = new DoubleList();
                StringTokenizer st = new StringTokenizer(line, " ", false);
                while (st.hasMoreTokens()) {
                    try {
                        row.add(Double.parseDouble(st.nextToken()));
                    } catch (NumberFormatException e) {
                        throw new IOException("Unknown format on line: " +
                                              lnr.getLineNumber());
                    }
                }
                rows.add(row);
            }

            double[][] matrix = new double[rows.size()][];
            for (int i = 0; i < matrix.length; i++) {
                matrix[i] = ( (DoubleList) rows.get(i)).toArray();
            }

            return matrix;
        } finally {
            if (lnr != null) {
                lnr.close();
            }
        }
    }


    public static boolean isValidNetworkAdress(String name) {
        if (!name.startsWith("net://")) return false;

        String address = name.substring(6);
        int indexOfSlash = address.indexOf('/');
        if (indexOfSlash == -1) return false;

//        String host = address.substring(0, indexOfSlash);
        try {
            Integer.parseInt(address.substring(indexOfSlash + 1));
        } catch (NumberFormatException e) {
            return false;
        }

        return true;
    }


    public static Object[] obtainStreams(String rawDataSource) throws IOException {
        InputStream input = null;
        OutputStream output = null;
        File file = new File(rawDataSource);
        
        if (rawDataSource.startsWith("net://")) {
            String address = rawDataSource.substring(6);
            int indexOfSlash = address.indexOf('/');
            String host = address.substring(0, indexOfSlash);
            String portString = address.substring(indexOfSlash + 1);
            int port;
            try {
                port = Integer.parseInt(portString);
            } catch (NumberFormatException e) {
                throw new NumberFormatException("Illegal port number: " +
                                                portString);
            }
            boolean connected = false;
            while (!connected) {
                try {
                    System.out.println("Attempting to Connect to DAQ Computer.");
                    connected = true;
                    Socket socket = new Socket(host, port);
                    input = socket.getInputStream();
                    output = socket.getOutputStream();
                } catch (SocketException e) {
                    connected = false;
                }
            }
            System.out.println("Connected to DAQ Computer.");
        } else if (file.isFile()) {
            input = new FileInputStream(file);
        } else if (file.isDirectory()) {
            input = new DirectoryInputStream(rawDataSource);
        } else {
            throw new IllegalArgumentException(
                rawDataSource +
                " is not a valid file, directory or network location");
        }
        
        Object[] streams = new Object[2];
        streams[0] = input;
        streams[1] = output;
        return streams;
    }


    /**
     * Prints the array to System.err in a single line.
     */
    public static void printArray(Object[] array) {
        for (int i = 0; i < array.length; i++) System.out.print(array[i] + ", ");
        System.out.println();
    }


    /**
     * Prints the array to System.err in a single line.
     */
    public static void printArray(double[] array) {
        for (int i = 0; i < array.length; i++) {
//            System.out.print(format(array[i], 3) + ", ");
            System.out.print(array[i] + ", ");
        }
        System.out.println();
    }


    public static void printArray(double[] array, int nFractDigits) {
        DecimalFormat f = new DecimalFormat();
        f.setMinimumFractionDigits(nFractDigits);

        for (int i = 0; i < array.length; i++) System.out.print(f.format(array[i]) + ", ");
        System.out.println();
    }


    /**
     * Prints the array to System.err in a single line.
     */
    public static void printArray(float[] array) {
        for (int i = 0; i < array.length; i++) System.out.print(array[i] + ", ");
        System.out.println();
    }


    /**
     * Prints the array to System.err in a single line.
     */
    public static void printArray(int[] array) {
        for (int i = 0; i < array.length; i++) System.out.print(array[i] + ", ");
        System.out.println();
    }


    public static String[] readNameValuePair(LineNumberReader r) throws IOException {
        String line = r.readLine();
        if (line == null) return new String[] {"", ""};
        
        StringTokenizer st = new StringTokenizer(line, "=");
        if (st.countTokens() != 2)
            throw new IOException("The format of each line sould be: name = value");

        return new String[] {st.nextToken().trim(), st.nextToken().trim()};
    }


    public static void saveTags(String description, String fileName, Tag[] tags) throws
        IOException {
        FileOutputStream fos = null;
        BufferedOutputStream bos = null;
        NeuroOutputStream nos = null;

        try {
            fos = new FileOutputStream(fileName);
            bos = new BufferedOutputStream(fos);
            nos = new NeuroOutputStream(bos);
            nos.writeString(description);
            for (int i = 0; i < tags.length; i++) nos.writeTag(tags[i]);
        } finally {
            if (nos != null) nos.close();
            if (bos != null) bos.close();
            if (fos != null) fos.close();
        }
    }


    public static final void writePublicFields(Object o, DataOutput out) throws IOException {

        Field[] f = o.getClass().getFields();
        for (int i = 0; i < f.length; i++) {
            try {				
                if((f[i].getModifiers() & (Modifier.FINAL + Modifier.PUBLIC)) 
                    == (Modifier.FINAL + Modifier.PUBLIC)) {

            } else if (f[i].getType().isEnum()) {
                    Enum en = (Enum) f[i].get(o);
                    out.writeInt(en != null ? en.ordinal() : -1);
                } else if (f[i].getType().equals(int.class)) {
                    int v = f[i].getInt(o);
                    out.writeInt(v);
                } else if (f[i].getType().equals(float.class)) {
                    float v = f[i].getFloat(o);
                    out.writeFloat(v);
                } else if (f[i].getType().equals(double.class)) {
                    double v = f[i].getDouble(o);
                    out.writeDouble(v);
                }
            } catch (IllegalAccessException ex) {
                throw new Error("Cannot access field " + f[i]);
            } catch (IllegalArgumentException ex) {
                throw new Error("Cannot access field " + f[i]);
            }
        }
    }

    /**
     * Read fields from in, assuming they come in the order determined introspectively from
     * the public fields of o.
     * 
     * @param o
     * @param in
     * @throws IOException
     */
    public static final void readPublicFields(Object o, DataInput in) throws IOException {
        Field[] f = o.getClass().getFields();
        for (int i = 0; i < f.length; i++) {
            try {
                if((f[i].getModifiers() & (Modifier.FINAL + Modifier.PUBLIC)) 
                        == (Modifier.FINAL + Modifier.PUBLIC)) {

                } else if (f[i].getType().isEnum()) {
                    int ordinal = in.readInt();
                    if (ordinal == -1) {
                        f[i].set(o, null);
                    } else {
                        f[i].set(o, f[i].getType().getEnumConstants()[ordinal]);
                    }
                } else if (f[i].getType().equals(int.class)) {
                    int v = in.readInt();
                    f[i].setInt(o, v);
                } else if (f[i].getType().equals(float.class)) {
                    float v = in.readFloat();
                    f[i].setFloat(o, v);
                } else if (f[i].getType().equals(double.class)) {
                    double v = in.readDouble();
                    f[i].setDouble(o, v);
                }
            } catch (IllegalAccessException ex) {
                throw new Error("Cannot access field " + f[i]);
            } catch (IllegalArgumentException ex) {
                throw new Error("Cannot access field " + f[i]);
            }
        }
    }


    public static final void printPublicFields(Object o) {
        Field[] f = o.getClass().getFields();
        for (int i = 0; i < f.length; i++) {
            try {
                System.err.println(f[i].getName() + ": " + f[i].get(o));
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }
    }
}