package edu.ucsc.neurobiology.vision.convert;

import java.io.*;

import edu.ucsc.neurobiology.vision.io.RawDataFile;

/**
 * This app repairs corrupt data.  It can take a header from one file, and paste it into another.  It
 * also copy ttls from one section of a file and paste it in a different location in the same file.
 * 
 * 
 * @author mgrivich
 *
 */
public class FrankenFile {

    public static void replaceHeader(String donorName, String recipientName, String newFileName, int headerSize) 
    throws IOException {

        File temp = new File(newFileName);
        if(temp.exists()) temp.delete();
        
        BufferedInputStream donor = new BufferedInputStream(
                new FileInputStream(donorName));
        BufferedInputStream recipient = new BufferedInputStream(
                new FileInputStream(recipientName));
        BufferedOutputStream  newFile = new BufferedOutputStream(
                new FileOutputStream(newFileName));

        byte[] b = new byte[headerSize];
        donor.read(b);
        newFile.write(b);
        recipient.read(b); //trash

        int bufferSize = 10*1024*1024;

        b = new byte[bufferSize];
        int bytesRead = bufferSize;
        while (bytesRead == bufferSize) {

            bytesRead = recipient.read(b);
            if (bytesRead > 0) {
                newFile.write(b, 0, bytesRead);
            }
        }
        newFile.flush();
        newFile.close();
        donor.close();
        recipient.close();
    }

    public static void pasteTTLs(String fromFileName, String toFileName, int headerSize, int bytesPerSample) throws IOException {
        RawDataFile rdf = new RawDataFile(fromFileName);
        short[] ttlData = new short[700*20000]; //hard coded where the good ttls are 
        rdf.getData(0, 200*20000, ttlData); //gets all the ttl data from 200 to 900.  This is good data.
        rdf.close();

        File temp = new File(toFileName);
        if(temp.exists()) temp.delete();

        BufferedInputStream fromFile = new BufferedInputStream(
                new FileInputStream(fromFileName));
        BufferedOutputStream toFile = new BufferedOutputStream(
                new FileOutputStream(toFileName));

        int ttlIndex = 13815;
        byte[] b = new byte[headerSize];
        fromFile.read(b);
        toFile.write(b);
        
        b = new byte[bytesPerSample];
        for(int i=0; i<20000*900; i++, ttlIndex++) {
            if(i<20000*10 || (i>152*20000 && i<156*20000)) { //hard coded where to replace ttls
                byte[] replacement = new byte[bytesPerSample];
                replacement[0] =  (byte)(ttlData[ttlIndex] >> 8);
                replacement[1] = (byte)(ttlData[ttlIndex] & 0xFF);
                toFile.write(replacement);
                fromFile.read(b);
            } else {
                fromFile.read(b);
                toFile.write(b);
            }

        }

        toFile.flush();
        fromFile.close();
        toFile.close();


    }
    public static void main(String args[]) throws Exception {  
/*		FrankenFile.replaceHeader("/Volumes/DiskG/data003/2009-10-23/data002/data002000.bin",
                "/Volumes/DiskG/data003/2009-10-23/data003/data003000.bin",
                "/Volumes/Rat/Data/Grivich/raw/2009-10-23/temp.bin", 134);
*/		
        FrankenFile.pasteTTLs("/Volumes/Rat/Data/Grivich/raw/2009-10-23/temp.bin",
                "/Volumes/Rat/Data/Grivich/raw/2009-10-23/data004/data004000.bin", 134, 98);
        /*
        int bytesPerSample = 1960000/20000;

        RawDataFile rdf = new RawDataFile("/Volumes/DiskG/data003/2009-10-23/data003/data003000.bin");
        short[] ttlData = new short[700*20000];
        rdf.getData(0, 200*20000, ttlData); //gets all the ttl data from 200 to 900.  This is good data.
        rdf.close();

        int ttlIndex = 0;  //determines where we are in the ttl array.  Must be incremented every sample.


        File f = new File("/Volumes/Rat/Data/Grivich/raw/2009-10-23/data004/data004000.bin");
        f.delete();
        BufferedInputStream donor = new BufferedInputStream(
                new FileInputStream("/Volumes/DiskG/data003/2009-10-23/data002/data002000.bin"));
        BufferedInputStream recipient = new BufferedInputStream(
                new FileInputStream("/Volumes/DiskG/data003/2009-10-23/data003/data003000.bin"));

        BufferedOutputStream frank = new BufferedOutputStream(
                new FileOutputStream("/Volumes/Rat/Data/Grivich/raw/2009-10-23/data004/data004000.bin"));




        byte[] b = new byte[headerSize];

        donor.read(b);
        frank.write(b); //write donor header
        recipient.read(b); //move pointer, trash data

        b = new byte[bytesPerSample];
        for(int i=0; i<20000*900; i++, ttlIndex++) {
            if(i<20000*10 || (i>152*20000 && i<156*20000)) {
                byte[] replacement = new byte[bytesPerSample];
                replacement[0] =  (byte)(ttlData[ttlIndex] >> 8);
                replacement[1] = (byte)(ttlData[ttlIndex] & 0xFF);
                frank.write(replacement);
            } else {
                recipient.read(b);
                frank.write(b);
            }

        }

        frank.flush();
        frank.close();
        donor.close();
        recipient.close();
        System.out.println("Done writing file.");
         */
    }
}
