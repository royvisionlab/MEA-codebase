package edu.ucsc.neurobiology.vision.io;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;

import com.martiansoftware.jsap.FlaggedOption;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.UnflaggedOption;

import edu.ucsc.neurobiology.vision.Vision;
import edu.ucsc.neurobiology.vision.anf.Join;
import edu.ucsc.neurobiology.vision.util.StringUtil;
import edu.ucsc.neurobiology.vision.util.VisionJSAP;

/*
 * Class EIMerger
 * 
 * This class is to merge EI files of multiple dataset runs into
 * a big ei file for the concatenated dataset.
 * The concept is inspired from Join, except that join doesn't behave well with EI
 * and just appends at end of file instead of matching IDs appropriately.
 * 
 * This is meant to be used during mVision concatenated analysis, to generate global EIs
 * and allow EI-based duplicate removal on the concatenated dataset
 * 
 * @author Vincent Deo - Stanford University - Nov 9th 2015
 */
public class EIMerger {

    /*
     * main method to allow call from command line Copied from Join's main
     */
    public static void main(String[] args) throws Exception {
        VisionJSAP jsap = new VisionJSAP(Join.class.getName(),
                new com.martiansoftware.jsap.Parameter[] {
                        new UnflaggedOption("in", JSAP.STRING_PARSER,
                                JSAP.REQUIRED, "Input folders."),
                        new UnflaggedOption("out", JSAP.STRING_PARSER,
                                JSAP.REQUIRED, "Output folder."), });
        JSAPResult parsedArgs = jsap.parse(args);

        String inputFolders = (String) parsedArgs.getString("in").trim();
        String[] folders = StringUtil.decomposeString(inputFolders, ";");
        String outputFolder = ((String) parsedArgs.getString("out")).trim();

        System.out.println("Starting EI Merging...");
        eiMerging(folders, outputFolder);
        System.out.println("EI Merging done!");
    }

    public static void eiMerging(String[] inputs, String outputFolder)
            throws IOException {
        new File(outputFolder).mkdirs();
        final int nFiles = inputs.length;

        PhysiologicalImagingFile[] eiFile = new PhysiologicalImagingFile[nFiles];
        int[][] idList = new int[nFiles][];

        for (int n = 0; n < nFiles; n++) {
            String dataSetName = new File(inputs[n]).getName();
            File f = new File(inputs[n] + File.separator + dataSetName + ".ei");
            if (IOUtil.isValidFile(f)) {
                eiFile[n] = new PhysiologicalImagingFile(inputs[n]
                        + File.separator + dataSetName + ".ei");
                idList[n] = eiFile[n].getIDList();
            } else {
                throw new FileNotFoundException();
            }
        }

        int[] mergedIDList = unionArrays(idList);
        Arrays.sort(mergedIDList);

        String dataSetName = new File(outputFolder).getName();
        PhysiologicalImagingFile newEIFile = new PhysiologicalImagingFile(
                outputFolder + File.separator + dataSetName + ".ei",
                eiFile[0].nlPoints, eiFile[0].nrPoints,
                eiFile[0].getArrayID() & 0xFFFF);
        int nPoints = newEIFile.nlPoints + newEIFile.nrPoints + 1;
        int nElectrodes = newEIFile.nElectrodes;
        for (int id = 0; id < mergedIDList.length; id++) {
            float[][][] bufferedImage = new float[2][nElectrodes][nPoints];
            int nSpikes = 0;
            for (int d = 0; d < nFiles; d++) {
                try {
                    float[][][] ei = eiFile[d].getImage(mergedIDList[id]);
                    int locSpikes = eiFile[d].getNSpikes(mergedIDList[id]);
                    nSpikes += locSpikes;
                    for (int e = 0; e < nElectrodes; e++) {
                        for (int t = 0; t < nPoints; t++) {
                            bufferedImage[0][e][t] += locSpikes * ei[0][e][t];
                            bufferedImage[1][e][t] += locSpikes * ei[1][e][t]
                                    * ei[1][e][t];
                        }
                    }
                } catch (Exception e) {
                    System.out.println("Error seeking id " + mergedIDList[id]
                            + " in file " + d);
                }
            }
            // Normalization loop
            for (int e = 0; e < nElectrodes; e++) {
                for (int t = 0; t < nPoints; t++) {
                    bufferedImage[0][e][t] /= (float) nSpikes;
                    bufferedImage[1][e][t] = (float) Math
                            .sqrt(bufferedImage[1][e][t] / (float) nSpikes);
                }
            }
            // Write
            System.out.println();
            newEIFile.appendImage(mergedIDList[id], nSpikes, bufferedImage);
        }
        for (int d = 0; d < eiFile.length; d++) {
            eiFile[d].close();
        }
        newEIFile.close();
    }

    /* Union of multiple arrays */
    public static int[] unionArrays(int[]... arrays) {
        int maxSize = 0;
        int counter = 0;

        for (int[] array : arrays)
            maxSize += array.length;
        int[] accumulator = new int[maxSize];

        for (int[] array : arrays)
            for (int i : array)
                if (!isDuplicated(accumulator, counter, i))
                    accumulator[counter++] = i;

        int[] result = new int[counter];
        for (int i = 0; i < counter; i++)
            result[i] = accumulator[i];

        return result;
    }

    public static boolean isDuplicated(int[] array, int counter, int value) {
        for (int i = 0; i < counter; i++)
            if (array[i] == value)
                return true;
        return false;
    }
}
