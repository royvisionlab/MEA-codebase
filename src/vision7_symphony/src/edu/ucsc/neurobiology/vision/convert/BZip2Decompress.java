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
 * See BZip2Compress for detailed comments.
 */
public class BZip2Decompress
    extends AbstractCalculation {

    private String inputPath, outputPath;
    private Boolean deleteOldFiles, fullTree;


    public void startCalculation() throws IOException {
        decompressFileOrDirectory(inputPath, outputPath);
        Vision.getInstance().getCalculationManager().calculationDone();
    }


    public void decompressFileOrDirectory(String inputPath, String outputPath) throws
        IOException {
        File inFile = new File(inputPath);
        if (inFile.isFile() && !inFile.isHidden() &&
            StringUtil.getExtension(inputPath).matches(".bz2")) {
            File outFile = null;
            if (outputPath.length() == 0) {
                outFile = new File(StringUtil.removeExtension(inputPath));
            } else {
                if (StringUtil.getExtension(outputPath).matches(".bz2")) {
                    outputPath = StringUtil.removeExtension(outputPath);
                }
                outFile = new File(outputPath);
                if (outFile.isDirectory()) {
                    throw new IllegalStateException(
                        "If input path is a file, output path must be as well.");
                }

            }

            decompressFile(inFile, outFile);
            if (deleteOldFiles) {
                inFile.delete();
            }
            //Is a file that is decompressed already or hidden.
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
                    BZip2Compress.copyFile(fromFile, toFile);
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

                    decompressFileOrDirectory(inputPath + "//" + fileNames[i],
                                              outputPath + "//" + fileNames[i]);
                    //It is a directory, recurse if requested.
                } else if (fullTree) {
                    (new File(outputPath + "//" + fileNames[i])).mkdir();
                    decompressFileOrDirectory(inputPath + "//" + fileNames[i],
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


    public static void decompressFile(File fromFile, File toFile) throws IOException {
        int bufferSize = 10 * 1024 * 1024;
        System.out.println("Decompressing: " + fromFile.getAbsolutePath() + " to " +
                           toFile.getAbsolutePath());
        CBZip2InputStream input = new CBZip2InputStream(new BufferedInputStream(new
            FileInputStream(fromFile), bufferSize));
        BufferedOutputStream output = new BufferedOutputStream(new FileOutputStream(
            toFile));
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


    public void setParameters(HashMap<String, String> parameters) {
        inputPath = ( (String) parameters.get("inputPath")).trim();
        outputPath = ( (String) parameters.get("outputPath")).trim();
        deleteOldFiles = Boolean.valueOf( (String) parameters.get("deleteOldFiles")).
                         booleanValue();
        fullTree = Boolean.valueOf( (String) parameters.get("fullTree")).
                   booleanValue();
    }


}
