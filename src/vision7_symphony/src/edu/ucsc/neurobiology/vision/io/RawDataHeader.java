package edu.ucsc.neurobiology.vision.io;

import java.io.IOException;
import java.io.InputStream;
import java.util.Date;

import com.martiansoftware.jsap.FlaggedOption;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Switch;
import com.martiansoftware.jsap.UnflaggedOption;

import edu.ucsc.neurobiology.vision.calculations.CalculationManager;
import edu.ucsc.neurobiology.vision.util.VisionJSAP;


/**
 * A container for the 512 (or new) raw data header.
 * Methods for reading a header from an InputStream
 * and writing a header to an OutputStream are also provided.
 * See CVS/formats/RawDataFileFormat.txt for a description of the fields.
 *
 * @author Charles A. Loomis, University of California, Santa Cruz
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public abstract class RawDataHeader {
    public static final int FORMAT_16BIT_UNCOMPRESSED = 0;
    public static final int FORMAT_12BIT_COMPRESSED = 1;
    public static final int FORMAT_8BIT_UNCOMPRESSED = 2;
    public static final int FORMAT_HUFFMAN_COMPRESSED = 3;


    public static int composeArrayID(int id, int part, int nParts) {
        return id + (part << 16) + (nParts << 24);
    }


    public static int extractArrayID(int id) {
        return id & 0xFFFF;
    }


    public static int extractElectrodeMapPart(int id) {
        return (id >> 16) & 0xFF;
    }


    public static int extractElectrodeMapPartsCount(int id) {
        return id >> 24;
    }


    public abstract int getSampleSize();


//    public abstract int getFileType();


    public abstract int getHeaderSize();


    public abstract int getNumberOfElectrodes();


    public abstract int getArrayID();


//    public int getRawArrayID() {
//        return
//            getArrayID() +
//            (getElectrodeMapPart() << 16) +
//            (getElectrodeMapPartsCount() << 24);
//    }


//    public abstract int getElectrodeMapPart();


//    public abstract int getElectrodeMapPartsCount();


    public abstract int getSamplingFrequency();


    public abstract int getNumberOfSamples();


    public abstract byte[] getBinaryRepresentation() throws IOException;


    public abstract String getComment();


    public double getTimePerSample() {
        return 1000.0 / getSamplingFrequency();
    }


    public abstract long getTime();


    public abstract int getFormat();


    public abstract void setFormat(int format);


    public abstract void setNumberOfSamples(int nSamples);


    public abstract String getDatasetIdentifier();


    static int unknownIndex = 1;
    public String getExperimentIdentifier() {
        String s = getDatasetIdentifier();
        if (s == null) {
            return null;
        } else {
            int i = s.lastIndexOf("-");
            if (i == -1) {
                return "unknown experiment";
            }
            return s.substring(0, i);
        }
    }


    public String getDatasetName() {
        String s = getDatasetIdentifier();
        if (s == null) {
            return null;
        } else {
            int i = s.lastIndexOf("-");
            if (i == -1) {
                return "unknown dataset " + unknownIndex++;
            }
            return s.substring(i + 1, s.length());
        }
    }

    
    /**
     * For debugging purposes.  Basic functionality of reading in a raw data stream and printing the header
     * should not be changed.
     * @param args
     * @author Peter H. Li, The Salk Institute
     * @throws IOException 
     */
    public static void main(String[] args) throws JSAPException, IOException {    	
        VisionJSAP jsap = new VisionJSAP("RawDataHeader", 
            new com.martiansoftware.jsap.Parameter[] {	
                new UnflaggedOption("rawdatapath", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.REQUIRED, JSAP.NOT_GREEDY, "Raw Data Path"),
            });
        JSAPResult parsedArgs = jsap.parseWithMangledNegs(args);
        
        String rawdataPath = parsedArgs.getString("rawdatapath");
        RawDataFile rawdata = new RawDataFile(rawdataPath);
        RawDataHeader header = rawdata.getHeader();
        System.out.println(header);
    }
}
