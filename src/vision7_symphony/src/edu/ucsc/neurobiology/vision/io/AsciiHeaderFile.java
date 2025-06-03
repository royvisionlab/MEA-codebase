package edu.ucsc.neurobiology.vision.io;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.HashMap;


/**
 * This class defined a generic file type with an Ascii Header.
 * Extend it to add functionality.
 *
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class AsciiHeaderFile {

    public static final String fieldSeparator = "\t";
    private static final int MAX_HEADER_SIZE = 100000;

    RandomAccessFile file;

    public int headerLength = 0;


    String header = "";
    String newline = System.getProperty("line.separator");
    HashMap<String, String[]> linesHash; //for reading, a HashMap where each element is a record from the file


    public AsciiHeaderFile(String fileName, boolean writing) throws IOException {
        if (writing == true) {
            createNewFile(fileName);
        } else {
            loadOldFile(fileName);
        }
    }


    public void createNewFile(String fileName) throws IOException {
        file = new RandomAccessFile(fileName, "rw");
        //If you overwrite a previous file with the same name, this line is necessary.
        file.setLength(0);
    }


    public int getHeaderLength() {
        return headerLength;
    }


    synchronized public void close() throws IOException {
        file.close();
    }


    public void addToHeader(String recordName, int value) {
        header = header + recordName + fieldSeparator + value + newline;
    }


    public void addToHeader(String nameValuePair) {
        header = header + nameValuePair + newline;
    }


    public void addToHeader(String recordName, float value) {
        header = header + recordName + fieldSeparator + value + newline;
    }


    public void addToHeader(String recordName, double value) {
        header = header + recordName + fieldSeparator + value + newline;
    }


    public void addToHeader(String recordName, String value) {
        header = header + recordName + fieldSeparator + value + newline;
    }


    public void addToHeader(String recordName, boolean value) {
        String stringValue;
        stringValue = (value) ? "TRUE" : "FALSE";
        header = header + recordName + fieldSeparator + stringValue + newline;

    }


    public void addToHeader(String recordName, int[] values) {
        String record = recordName;
        for (int value : values) {
            record = record + fieldSeparator + value;
        }
        header = header + record + newline;
    }


    //This function must be called after all lines have been added to the header.
    public void writeHeader() throws IOException {
        header = header + newline;

        file.writeBytes(header);
        //used by get Values, which is used by functions that want values
        //from new as well as old files.
        makeLinesHash(header);

    }


    public void makeLinesHash(String header) {
        String[] lines = header.split(newline);
        linesHash = new HashMap<String, String[]>();
        for (String line : lines) {
            System.out.println(line);
            String[] record = line.split(fieldSeparator);

            String[] fields = new String[record.length - 1];
            for (int j = 0; j < record.length - 1; j++) {
                fields[j] = record[j + 1];
            }
            linesHash.put(record[0], fields);
        }

    }



    public void loadOldFile(String fileName) throws IOException {
        BufferedInputStream in = new BufferedInputStream(new FileInputStream(new File((fileName))));
//		boolean done = false;
        char c = '0';
        headerLength = 0;
        do {
            c = (char) in.read();
            headerLength++;
            if(c =='\r') {
                c = (char) in.read();
                headerLength++;
                
                if(c=='\r') {
                    newline = "\r";  //Mac OS 9 end of line character
                    break;
                } 
                
                if(c=='\n') {
                    c = (char) in.read();
                    headerLength++;
                    if(c == '\r') {
                        c = (char) in.read();
                        headerLength++;
                        if(c=='\n') {
                            newline = "\r\n"; //Windows end of line character
                            System.out.println("Windows");
                            break;
                        }
                    }
                }
            }

            if(c =='\n') {
                c = (char) in.read();
                headerLength++;
                if(c=='\n') {
                    newline = "\n";  //Unix and Mac OS X end of line character
                    break;
                }
            }


        } while(headerLength< MAX_HEADER_SIZE);
        if(headerLength >= MAX_HEADER_SIZE) {
            throw new IllegalStateException("File is not recognized.  Max header size exceeded.");
        }
        
        in.close();
        
        file = new RandomAccessFile(fileName, "r");
        //       byte buffer[] = new byte[headerHeaderLength];
    

        byte[] buffer = new byte[headerLength];
        file.seek(0);
        file.readFully(buffer);
        String header = new String(buffer);
        makeLinesHash(header);
    }


    //returns the String[] values for the String recordName.
    //The calling function must cast the value into the correct type.
    //String[] has the length of the number of arguments in the record.
    String[] getValues(String recordName) {
        String[] value = linesHash.get(recordName);
        if (value == null) {
            throw new IllegalStateException("Record, " + recordName +
            ", not found in file.");
        }
        return value;
    }


    int getIntParam(String recordName) {
        String value = linesHash.get(recordName)[0];
        if (value == null) {
            throw new IllegalStateException("Record, " + recordName +
            ", not found in file.");
        }
        return Integer.parseInt(value);
    }

}
