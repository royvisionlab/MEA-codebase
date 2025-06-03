package edu.ucsc.neurobiology.vision.anf;

import java.io.*;
import java.util.*;

import com.martiansoftware.jsap.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.calculations.*;
import edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.io.chunk.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * Calculation that joins the results of several automated neuron findings.
 * The neuron, sta, ei and params files are joined into single files which can be viewed by
 * the NeuronViewer.  Globals file is hacked into proper shape too.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class Join
    extends AbstractCalculation {

    String[] folder;
    String outputFolder;


    public void startCalculation() throws Exception {
        new File(outputFolder).mkdirs();

        final int nFiles = folder.length;
        NeuronFile[] neuronFile = new NeuronFile[nFiles];
        STAFile[] staFile = new STAFile[nFiles];
        ParametersFile[] paramsFile = new ParametersFile[nFiles];
        PhysiologicalImagingFile[] eiFile = new PhysiologicalImagingFile[nFiles];

        // Open files from each input folder to be joined
        int nNeurons = 0;
        boolean staFilesExist = true;
        boolean paramsFilesExist = true;
        boolean eiFilesExist = true;
        for (int n = 0; n < nFiles; n++) {
            String dataSetName = new File(folder[n]).getName();
            neuronFile[n] = new NeuronFile(
                folder[n] + File.separator + dataSetName + ".neurons");

            File f = new File(folder[n] + File.separator + dataSetName + ".sta");
            if (IOUtil.isValidFile(f)) {
                staFile[n] = new STAFile(f);
            } else {
                staFilesExist = false;
            }

            f = new File(folder[n] + File.separator + dataSetName + ".params");
            if (IOUtil.isValidFile(f)) {
                paramsFile[n] = new ParametersFile(f);
            } else {
                paramsFilesExist = false;
            }
            
            f = new File(folder[n] + File.separator + dataSetName + ".ei");
            if (IOUtil.isValidFile(f)) {
                eiFile[n] = new PhysiologicalImagingFile(folder[n] + File.separator + dataSetName + ".ei");
            } else {
                eiFilesExist = false;
            }

            nNeurons += neuronFile[n].getNumberOfNeurons();
        }

        // Create output folder, files
        String dataSetName = new File(outputFolder).getName();

        NeuronFile newNeuronFile = new NeuronFile(
            outputFolder + File.separator + dataSetName + ".neurons",
            neuronFile[0].getHeader(),
            2 * nNeurons,
            neuronFile[0].getTTLTimes());

        STAFile newSTAFile = null;
        if (staFilesExist) {
            newSTAFile = new STAFile(
                outputFolder + File.separator + dataSetName + ".sta",
                2 * nNeurons,
                staFile[0].getWidth(),
                staFile[0].getHeight(),
                staFile[0].getSTADepth(),
                staFile[0].getSTAOffset(),
                staFile[0].getStixelWidth(),
                staFile[0].getStixelHeight(),
                staFile[0].getRefreshTime()
                         );
        }
        
        GlobalsFile oldGlobals = new GlobalsFile(folder[0] + File.separator + new File(folder[0]).getName() + ".globals", ChunkFile.READ);
        GlobalsFile newGlobals = new GlobalsFile(outputFolder + File.separator + dataSetName + ".globals", ChunkFile.WRITE);
        if (oldGlobals.imageCalibrationParamsExists()) {
            GlobalsFile.ImageCalibrationParams calParams = oldGlobals.getImageCalibrationParams();
            calParams.arrayNParts = 1;
            calParams.arrayPart = 1;
            newGlobals.setImageCalibrationParams(calParams);
        }
        if (oldGlobals.runTimeMovieParamsExists()) {
            newGlobals.setRunTimeMovieParams(oldGlobals.getRunTimeMovieParams());
        }
        oldGlobals.close();
        newGlobals.close();
 //       edu.ucsc.neurobiology.vision.convert.BZip2Compress.copyFile(
 //       		new File(folder[0] + File.separator + new File(folder[0]).getName() + ".globals"),
 //       		new File(outputFolder + File.separator + dataSetName + ".globals"));

        ParametersFile newParamsFile = null;
        if (paramsFilesExist) {
            newParamsFile = new ParametersFile(
                outputFolder + File.separator + dataSetName + ".params",
                paramsFile[0].getColumnNamesAndTypes(),
                2 * nNeurons);
        }
        
        PhysiologicalImagingFile newEIFile = null;
        if (eiFilesExist) {
            newEIFile = new PhysiologicalImagingFile(
                    outputFolder+File.separator+dataSetName + ".ei",
                    eiFile[0].nlPoints, eiFile[0].nrPoints, eiFile[0].getArrayID() & 0xFFFF); 
        }

        
        // join
        int id = 1;
        for (int n = 0; n < nFiles; n++) {
            int[] idList = neuronFile[n].getIDList();
            
            // If ei file is available, use this to determine what the original electrode numbers were.
            // This info is not stored anywhere else.  This cannot be easily fixed because of file
            // format lock in.
            int[] parentElectrodeNumbers = null;
            if (eiFilesExist) {
                parentElectrodeNumbers = ElectrodeMapFactory.getElectrodeMap(eiFile[n].getArrayID()).getParentElectrodeNumbers();
            } 
  
            for (int i = 0; i < idList.length; i++) {
                System.out.println(id);
                int electrode = neuronFile[n].getElectrode(idList[i]);
                if (eiFilesExist) {
                    electrode = parentElectrodeNumbers[electrode-1];
                }
                newNeuronFile.addNeuron(
                    electrode,
                    id,
                    neuronFile[n].getSpikeTimes(idList[i]),
                    neuronFile[n].getSpikeCount(idList[i])
                    );

                if (staFilesExist) {
                    try {
                        newSTAFile.addSTA(id, staFile[n].getSTA(idList[i]));
                    } catch (IOException ex) {
                        System.err.println("Problem reading STA " + idList[i] +
                                           " from segment " + n);
                        ex.printStackTrace();
                       // throw ex;
                    }
                }

                if (paramsFilesExist) {
                    Object[] row = null;
                    try {
                        row = paramsFile[n].getRow(idList[i]);
                    } catch (Exception ex) {
                        System.err.println("Problem reading Params for neuron " +
                                           idList[i] +
                                           " from segment " + n);
                        throw ex;
                    }
                    row[0] = new Double(id);
                    newParamsFile.addRow(row);
                }
                
                if (eiFilesExist) {  	
                    try {
                        newEIFile.appendImage(id,eiFile[n].getNSpikes(idList[i]), eiFile[n].getExpandedImage(idList[i]));
                    } catch (Exception ex) {
                        
                        System.err.println("Problem reading EI for neuron " + idList[i] + " from segment " + n);
                        ex.printStackTrace();
                        throw ex;
                    }
                }

                id++;
            }
        }

        // close everything
        for (int n = 0; n < nFiles; n++) {
            neuronFile[n].close();
            if (staFilesExist) {
                staFile[n].close();
            }
            if (paramsFilesExist) {
                paramsFile[n].close(false);
            }
            
            if (eiFilesExist) {
                eiFile[n].close();
            }
        }
        newNeuronFile.close();
        if (staFilesExist) {
            newSTAFile.close();
        }
        if (paramsFilesExist) {
            newParamsFile.close(true);
        }
        
        if (eiFilesExist) {
            newEIFile.close();
        }

        Vision.getInstance().getCalculationManager().calculationDone();
    }


    public void setParameters(HashMap<String, String> parameters) {
        String inputFolders = ( (String) parameters.get("Input Folders")).trim();
        folder = StringUtil.decomposeString(inputFolders, ";");

        outputFolder = ( (String) parameters.get("outputFolder")).trim();
    }


    public static void main(String[] args) throws Exception {
        VisionJSAP jsap = new VisionJSAP( 
                Join.class.getName(), 
                new com.martiansoftware.jsap.Parameter[] {
                        new UnflaggedOption("in", JSAP.STRING_PARSER, JSAP.REQUIRED, "Input folders."),
                        new UnflaggedOption("out", JSAP.STRING_PARSER, JSAP.REQUIRED, "Output folder."),	
                        new FlaggedOption( "config", JSAP.STRING_PARSER, "config.xml", JSAP.NOT_REQUIRED, 'c', "config", 
                        "Configuration file" )
                }
        );
        JSAPResult parsedArgs = jsap.parse(args);

        HashMap<String,String> p = new HashMap<String,String>();
        p.clear();
        p.put("Input Folders", parsedArgs.getString("in"));
        p.put("outputFolder", "" + parsedArgs.getString("out"));
        Vision.getInstance().getCalculationManager().runCalculation("Join", p);
    }

}
