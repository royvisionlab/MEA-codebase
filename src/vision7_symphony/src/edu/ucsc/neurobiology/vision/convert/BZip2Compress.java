package edu.ucsc.neurobiology.vision.convert;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.calculations.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Matthew Grivich, University of California, Santa Cruz
 *
 *
 Compress file or folder using bzip2 (String fromPath, String toPath, boolean deleteOldFiles, boolean fullTree);
 If given a file, it will compress it to bzip2.  fileName.ext -> fileName.ext.bz2
 If given a directory, it will compress the contents of the directory, file by file.
     fromDirectory/fileName1.ext -> toDirectory/fileName1.ext.bz2
     fromDirectory/fileName2.ext -> toDirectory/fileName2.ext.bz2
    etc.

 fromPath and toPath can be the same.  If toPath is blank, it will be assumed to be equal to fromPath
 dot files (hidden Mac OS files) will be ignored.  Files less than 1 megabyte in size will be ignored.
 If deleteOldFiles is true, the uncompressed files will be deleted as each file is compressed.
 If fullTree is true, then all (non-hidden) files in the directory tree will be compressed.
 The tree structure will be left intact.

 Decompress file or folder using bzip2 (String fromPath, String toPath, boolean deleteOldFiles, boolean recursive);
 This function does the same as Compress file or folder, but in reverse.

 As data is compressed on jackstraw or SRB, it will have the form directory structure
 as now, except that a .bz2 will be appended to each file.

 These functions can be replicated using mac scripts, but please follow the same
 (or some compatible) specification.  You should note that these scripts will
 work with processed data or anything else, as desired.

 This functions can be called from automatic Java or Mac scripts, as desired.

 Current relevant test numbers:
 Copying a 120 second file, from one disk to another, using java: 87s
 Compressing a file: 1490s
 Decompressing a file: 726s
 Compression rate: 35%

 */
public class BZip2Compress
    extends AbstractCalculation {

    private String inputPath, outputPath;
    private Boolean deleteOldFiles, fullTree;
    private long minimumSize;


    public void startCalculation() throws IOException {
        compressFileOrDirectory(inputPath, outputPath);

        Vision.getInstance().getCalculationManager().calculationDone();
    }


    public void compressFileOrDirectory(String inputPath, String outputPath) throws
        IOException {
        File inFile = new File(inputPath);
        if (inFile.isFile() && !inFile.isHidden() && inFile.length() > minimumSize) {
            File outFile = null;
            if (outputPath.length() == 0) {
                outFile = new File(inputPath + ".bz2");
            } else {
                if (StringUtil.getExtension(outputPath).matches(".bz2")) {
                    outFile = new File(outputPath);
                } else {
                    outFile = new File(outputPath + ".bz2");
                }
            }
            compressFile(inFile, outFile);
            if (deleteOldFiles) {
                inFile.delete();
            }
            //Is a file that is too small, compressed already, or hidden.
        } else if (inFile.isFile()) {
            if (outputPath.length() == 0) {
                outputPath = inputPath;
            }
            if (StringUtil.getExtension(outputPath).matches(".bz2")) {
                outputPath = StringUtil.removeExtension(outputPath);
            }
            File fromFile = new File(inputPath);
            File toFile = new File(outputPath);
            if (! (fromFile.compareTo(toFile) == 0)) {
                if (deleteOldFiles) {
                    System.out.println("Moving: " + fromFile.getAbsolutePath() + " to " +
                                       toFile.getAbsolutePath());
                    inFile.renameTo(toFile);
                } else {
                    copyFile(fromFile, toFile);
                }
            }
        } else if (inFile.isDirectory()) {
            if (outputPath.length() == 0) {
                outputPath = inputPath;
            }
            File outFile = new File(outputPath);
            if (outFile.isFile()) {
                throw new IllegalStateException(
                    "If input path is a directory, output path must be as well");
                //else if directory does not exist.
            } else if (!outFile.isDirectory()) {
                outFile.mkdirs();
            }
            String[] fileNames = inFile.list();
            for (int i = 0; i < fileNames.length; i++) {

                if ( (new File(inputPath + "//" + fileNames[i])).isFile()) {

                    compressFileOrDirectory(inputPath + "//" + fileNames[i],
                                            outputPath + "//" + fileNames[i]);
                    //It is a directory, recurse if requested.
                } else if (fullTree) {
                    (new File(outputPath + "//" + fileNames[i])).mkdir();
                    compressFileOrDirectory(inputPath + "//" + fileNames[i],
                                            outputPath + "//" + fileNames[i]);
                    if (deleteOldFiles) {
                        inFile.delete();
                    }

                }
            }

        } else {
            System.out.println("Error: " + inFile + " not found.");
        }

    }


    public static void compressFile(File fromFile, File toFile) throws IOException {
        System.out.println("Compressing: " + fromFile.getAbsolutePath() + " to " +
                           toFile.getAbsolutePath());
        int bufferSize = 10 * 1024 * 1024;
        BufferedInputStream input = new BufferedInputStream(new FileInputStream(fromFile),
            bufferSize);
        CBZip2OutputStream output = new CBZip2OutputStream(
            new BufferedOutputStream(new FileOutputStream(toFile), bufferSize));

        byte[] b = new byte[bufferSize];
        int bytesRead = bufferSize;
        while (bytesRead == bufferSize) {

            bytesRead = input.read(b);
            output.write(b, 0, bytesRead);

        }

        output.flush();
        output.close();
        input.close();

    }


    public static void copyFile(File fromFile, File toFile) throws IOException {
        System.out.println("Copying: " + fromFile.getAbsolutePath() + " to " +
                           toFile.getAbsolutePath());
        int bufferSize = 10 * 1024 * 1024;
        BufferedInputStream input = new BufferedInputStream(new FileInputStream(fromFile),
            bufferSize);
        BufferedOutputStream output = new BufferedOutputStream(new FileOutputStream(
            toFile), bufferSize);

        byte[] b = new byte[bufferSize];
        int bytesRead = bufferSize;
        while (bytesRead == bufferSize) {

            bytesRead = input.read(b);
            if (bytesRead > 0) {
                output.write(b, 0, bytesRead);
            }
        }
        output.flush();
        output.close();
        input.close();

    }


    public void setParameters(HashMap<String, String> parameters) {
        inputPath = ( (String) parameters.get("inputPath")).trim();
        outputPath = ( (String) parameters.get("outputPath")).trim();
        deleteOldFiles = Boolean.valueOf( (String) parameters.get("deleteOldFiles")).
                         booleanValue();
        fullTree = Boolean.valueOf( (String) parameters.get("fullTree")).
                   booleanValue();
        minimumSize = Integer.valueOf( (String) parameters.get("minimumSize")).intValue();
    }

}
