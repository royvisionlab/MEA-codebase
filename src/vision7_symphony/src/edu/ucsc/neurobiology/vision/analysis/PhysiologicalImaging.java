package edu.ucsc.neurobiology.vision.analysis;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.calculations.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * This can either calculate a single EI, if neuronIDs is given, or it can calculate all EIs, if neuronIDs < 0.
 * 
 * This subtracts the mean as the first 10 samples on each electrode, so you better put at least 10 left samples...
 * 
 * If the Append option is true, and there is an existing EI file in the neuronFile directory, then this will try to
 * append EIs to that EI file, overwriting EIs for existing Neuron IDs.  If the Append option is false and there is
 * an existing EI file, the file is clobbered.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class PhysiologicalImaging extends AbstractCalculation {
    private static final int MEAN_LSAMPLES = 10;
    
    private String neuronFileName, rawFileName;
    private NeuronFile neurons;
    private int nlPoints, nrPoints;
    private int nSpikesToAverage;
    private int[] neuronIDs;
    private boolean append;
    private boolean recalc;
    
    private RawDataHeader header;
    private int nElectrodes;
    private int[] neuronList;

    private Vision app;
    private RawDataFile rawData;
    private PhysiologicalImagingFile eifile;


    public void startCalculation() throws Exception {
        app = Vision.getInstance();
        app.sendMessage("Preparing...");
        
        rawData = new RawDataFile(new File(rawFileName));
        header = rawData.getHeader();
        nElectrodes = header.getNumberOfElectrodes();
        neurons = new NeuronFile(neuronFileName);
        neuronList = neurons.getIDList();

        int[] idsToCalculate;
        if (neuronIDs.length == 1 && neuronIDs[0] < 0) idsToCalculate = neuronList;
        else                                           idsToCalculate = neuronIDs;        
        
        String fileName = StringUtil.removeExtension(new File(neuronFileName).getAbsolutePath()) + VisionParams.EI_FILE_EXTENSION;
        File file = new File(fileName);
        boolean existingEIFile = file.exists();
        if (!existingEIFile)
            eifile = new PhysiologicalImagingFile(fileName, nlPoints, nrPoints, header.getArrayID());
        else if (append) {
            eifile = new PhysiologicalImagingFile(fileName);
            if (nlPoints != eifile.nlPoints || nrPoints != eifile.nrPoints || header.getArrayID() != eifile.arrayID)
                throw new IllegalArgumentException("Existing EI file is not compatible and cannot be appended to.");
        } else {
            file.delete();
            existingEIFile = false;
            eifile = new PhysiologicalImagingFile(fileName, nlPoints, nrPoints, header.getArrayID());
        }
                
        
        WaveformCalculator average = new WaveformCalculator(nElectrodes, nlPoints, nrPoints);
        app.startProgressBar();
        for (int index = 0; index < idsToCalculate.length; index++) {
            int id = idsToCalculate[index];
            
            if (!neurons.containsID(id)) throw new IllegalArgumentException("Could not find id " + id + " in Neurons File!");
            if (existingEIFile && !recalc && eifile.getIndex(id) != -1) {
                app.sendMessage("Skipping existing Neuron ID " + id);
                continue;
            } 
            app.sendMessage("Processing Neuron " + index + "/" + idsToCalculate.length + ", id " + id);

            int[] times = neurons.getSpikeTimes(id);
            average.reset();
            final int nSpikes = Math.min(nSpikesToAverage, times.length);

            short[][] samples = average.createCompatibleBuffer();
            int nSamples = rawData.getHeader().getNumberOfSamples();
            for (int i = 0; i < nSpikes; i++) {
                double time = times[i];
                final int startSample = (int) Math.round(time) - nlPoints - 1;
                if (startSample > 0 && startSample < nSamples - (samples.length)) {
                    rawData.getData(startSample, samples);
                    average.addSpike(time, samples);
                }
            }
            average.finish();
            final float[][] image = average.getAverage();
            final float[][] error = average.getError();

            // remove the mean
            for (int electrode = 0; electrode < nElectrodes; electrode++) {
                double mean = 0;
                for (int i = 0; i < MEAN_LSAMPLES; i++) mean += image[electrode][i];
                mean /= MEAN_LSAMPLES;
                for (int i = 0; i < image[electrode].length; i++) image[electrode][i] -= mean;
            }
            
            if (existingEIFile) // Already tested for recalc at top of block
                eifile.overwriteOrAppendImage(id, average.getNSpikes(), new float[][][] {image, error});
            else
                eifile.appendImage(id, average.getNSpikes(), new float[][][] {image, error});
            app.setProgress(100 * (index+1) / idsToCalculate.length);
        }
        app.endProgressBar();

        rawData.close();
        eifile.close();

        app.sendMessage("Electrophysiological Imaging done");
        Vision.getInstance().getCalculationManager().calculationDone();
    }
    
    
    public void setParameters(HashMap<String, String> parameters) {
        String datasetFolder = parameters.get("Dataset Folder");
        neuronFileName = datasetFolder + File.separator + new File(datasetFolder).getName() + VisionParams.NEURON_FILE_EXTENSION;

        rawFileName = parameters.get("Raw Data File");
        nlPoints         = Integer.parseInt(parameters.get("Left Samples"));
        nrPoints         = Integer.parseInt(parameters.get("Right Samples"));
        nSpikesToAverage = Integer.parseInt(parameters.get("Spikes To Average"));
        append           = Boolean.parseBoolean(parameters.get("Append"));
        recalc           = Boolean.parseBoolean(parameters.get("Recalculate"));
        
        neuronIDs = parseNeuronIDs(parameters.get("Neuron ID"));
        
        if (nlPoints < MEAN_LSAMPLES)
            System.err.println("EI calculation uses the first " + MEAN_LSAMPLES + 
                    " samples on each electrode to subtract mean; recommend you select > " + 
                    MEAN_LSAMPLES + " for nlPoints; you selected " + nlPoints);
        
        System.err.println("nlPoints " + nlPoints);
        System.err.println("nrPoints " + nrPoints);
        System.err.println("nSpikesToAverage " + nSpikesToAverage);
        System.err.println("Neuron ID " + neuronIDs);
    }
    
    private int[] parseNeuronIDs(String idString) {
        String[] idStrings = idString.split(",");
        int[] ids = new int[idStrings.length];
        int i = 0;
        for (String s : idStrings) ids[i++] = Integer.parseInt(s);
        return ids;
    }
}