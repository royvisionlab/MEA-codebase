package edu.ucsc.neurobiology.vision.analysis;

import java.io.*;
import java.text.*;
import java.util.*;

import javax.swing.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.calculations.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.parameters.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * Forked from Imaging.java, the original Fast EI calculation.  
 * 
 * I am trying to optimize for the situation where we are running EIs on a large 
 * raw data directory but only using portion of the spikes, e.g. 1000 spikes.  This
 * is common, and especially critical for online analysis.  The goals are:
 * 
 * 	1) Clean up the spike selection code.  This was inefficient and about 25% of time
 * was being spent in these methods.  Changing to a setup where the spike stream is 
 * presorted and already has skips incorporated.
 * 
 *  2) Avoid processing samples that are not used.  Ultimately it would be nice not to 
 * even read samples that are not needed.  For now, just try to keep processing of these
 * to a minimum, potentially multithread more of this processing.
 * 
 * In the process, I discovered that the splining being done during EI analysis was pointless,
 * so that has also been removed.
 *
 * @author Peter H. Li, The Salk Institute
 */
public class ImagingNew
    extends AbstractCalculation implements SampleListener {

    private String neuronFileName, rawDataSource;
    private int nlPoints, nrPoints, nPoints;
    private int nElectrodes;
    private long nSpikes;
    private long sampleIndex;
    private LinkedList<short[]> sampleList;
    private Imaging.DataContainer spikeData;
    private NeuronFile.OrderlySkippingSpikeStream spikeIterator;
    private NeuronFile.ExtendedSpike currentSpike;
    private float[] mean;
    private boolean movingBaselineSubtraction;
    private double alpha;

    private long spikeNumber = 1; // DO NOT CHANGE, the algorithm needs a starting value of 1
    private int spikesToUse;
    private MultipleCompressedSampleInputStream sampleInputStream;
    private RawDataHeader header;
    private NeuronFile neuronFile;
    private Vision app;
    private long t1;
    private JPanel diagnosticPanel;
    private java.util.Timer timer;
    private JLabel samplesReadLabel;
    private int numThreads;
    private int[] edgeIDs;
    private EICalculator[] eiCalcs;
    
    private PhysiologicalImagingFile eiFile;


    public ImagingNew() {
        if (Vision.isGUIBased())
            diagnosticPanel = new JPanel(new TableLayout(0, 0));
    }


    public void startCalculation() throws IOException {
        app = Vision.getInstance();
        app.sendMessage("Preparing...");

        // get the raw data file
        sampleInputStream = new MultipleCompressedSampleInputStream(rawDataSource);
        header = sampleInputStream.getHeader();
        nElectrodes = header.getNumberOfElectrodes();

        mean = new float[nElectrodes];

        neuronFile = new NeuronFile(neuronFileName);
        spikeIterator = new NeuronFile.OrderlySkippingSpikeStream(neuronFile, spikesToUse);
        nSpikes = spikeIterator.getNumberOfSpikes();

        try {
            String eiFileName = StringUtil.removeExtension(neuronFileName) + VisionParams.EI_FILE_EXTENSION;
            eiFile = new PhysiologicalImagingFile(eiFileName, nlPoints, nrPoints, header.getArrayID());
        } catch (IOException e) {
            e.printStackTrace();
        }

        // load the first spike
        currentSpike = spikeIterator.next();
        spikeNumber++;

        spikeData = new Imaging.DataContainer();
        spikeData.spike = new short[nPoints + 2][nElectrodes];

        // Flexible threading, could rewrite to use abstracted ParallelUtils for edgeIDs
        int i1 = 0;
        int i2;
        int[] idList = neuronFile.getIDList();
        int istep = idList.length / numThreads;
        edgeIDs = new int[numThreads];
        eiCalcs = new EICalculator[numThreads];
        for (int i = 0; i < numThreads; i++) {
            i2 = istep * (i+1);
            if (i == (numThreads - 1)) i2 = idList.length;
            
            int[] l = new int[i2 - i1];
            System.arraycopy(idList, i1, l, 0, l.length);
            edgeIDs[i] = l[l.length-1];
            eiCalcs[i] = new EICalculator(nlPoints, nrPoints, l, nElectrodes);
            
            i1 = i2;
        }

        // create the sample buffer
        sampleList = new LinkedList<short[]>();
        for (int i = 0; i < nPoints; i++)
            sampleList.add(new short[nElectrodes]);

        app.startProgressBar();
        timer = new java.util.Timer(true);
        timer.scheduleAtFixedRate(timerTask, 0, 100);

        // start
        t1 = System.currentTimeMillis();
        sampleInputStream.addSampleListener(this);
        sampleInputStream.start();
        for (EICalculator e : eiCalcs) e.start();

        app.sendMessage("Calculating EI...");
    }


    TimerTask timerTask = new TimerTask() {
        public void run() {
            int samplesRead = sampleInputStream.getSamplesRead();
            app.setProgress( (int) (100.0 * samplesRead / header.getNumberOfSamples()));

            // show GUI info if appropriate
            if (Vision.isGUIBased()) {
                NumberFormat f = DecimalFormat.getInstance();
                samplesReadLabel.setText(f.format(sampleInputStream.getSamplesRead()));
            }

            // finish the calculation
            if (sampleInputStream.isFinished()) {
                for (EICalculator e : eiCalcs) e.finish();
                while (anyEICalcsAlive()) {
                    try {
                        Thread.sleep(100);
                    } catch (InterruptedException e) {}
                }

                // save the covariance matrices
                try {
//                    String eiFileName = StringUtil.removeExtension(neuronFileName) +
//                                        ".ei";
//                    PhysiologicalImagingFile eiFile = new PhysiologicalImagingFile(
//                        eiFileName, nlPoints, nrPoints,
//                        sampleInputStream.getHeader().getArrayID());

                    for (EICalculator e : eiCalcs) {
                        for (Integer id : e.calc.keySet()) {
                            WaveformCalculator calc = e.calc.get(id);
                            calc.finish();
                            eiFile.appendImage(id, calc.getCount(), calc.getImage());
                        }
                    }
                    eiFile.close();

                } catch (IOException e) {
                    e.printStackTrace();
                }

                long t2 = System.currentTimeMillis();
                app.sendMessage("Done in: " + (t2 - t1) / 1000. + " s.");

                timer.cancel();
                app.endProgressBar();
                calculationDone();
            }
        }
    };

    
    private boolean anyEICalcsAlive() {
        for (EICalculator e : eiCalcs) {
            if (e.isAlive()) return true;
        }
        return false;
    }
    
    
    public void processSample(short[] sample) {
        if (spikeNumber >= nSpikes) return; // make sure we don't get a infinite loop
        
        // Update the mean for moving average filtering
        if (movingBaselineSubtraction) {
            for (int i = 0; i < nElectrodes; i++)
                mean[i] += alpha * (sample[i] - mean[i]);
        }
        
        // Don't bother to further process samples we aren't going to use
        if (currentSpike.time - sampleIndex > nlPoints) {
            sampleIndex++;
            return;
        }
        
        // copy the sample into sampleList
        short[] s = sampleList.removeFirst();
        System.arraycopy(sample, 0, s, 0, nElectrodes);

        // filter the sample
        if (movingBaselineSubtraction) {
            for (int i = 0; i < nElectrodes; i++)
                s[i] -= mean[i];
        }
    
        // Add the sample to the queue
        sampleList.addLast(s);
        
        // process possible spikes
        while (sampleIndex - currentSpike.time == nrPoints) {
            for (int i = 0; i < nPoints; i++)
                spikeData.spike[i] = sampleList.get(i);

            spikeData.neuronID = currentSpike.neuronID;
            spikeData.time = currentSpike.time;

            for (int i = 0; i < edgeIDs.length; i++) {
                if (spikeData.neuronID <= edgeIDs[i]) {
                    eiCalcs[i].add(spikeData);
                    break;
                }
            }

            currentSpike = spikeIterator.next();
            spikeNumber++;
            if (spikeNumber >= nSpikes) break; // make sure we don't get a infinite loop
        }

        sampleIndex++;
    }
    

    public void finishSampleProcessing() {} 

    
    public void setParameters(HashMap<String, String> parameters) {
        String datasetFolder = parameters.get("Dataset Folder");
        neuronFileName = datasetFolder + File.separator + new File(datasetFolder).getName() + ".neurons";
        rawDataSource = parameters.get("Raw Data File");
        nlPoints = Integer.parseInt(parameters.get("Left Samples"));
        nrPoints = Integer.parseInt(parameters.get("Right Samples"));
        this.nPoints = nlPoints + nrPoints + 1;

        spikesToUse = Integer.parseInt(parameters.get("Spikes To Average"));

        numThreads = Integer.parseInt( (String) parameters.get("nThreads"));

        movingBaselineSubtraction = Boolean.parseBoolean(parameters.get("Moving Baseline Subtraction"));
        double timeConstant = Double.parseDouble(parameters.get("Mean Time Constant"));
        alpha = (float) (1.0 / (timeConstant * 20000.0));

        if (Vision.isGUIBased()) {
            // create the diagnostic panel
            diagnosticPanel.add(new JLabel("Samples read:"));
            samplesReadLabel = new JLabel("0", JLabel.RIGHT);
            diagnosticPanel.add(samplesReadLabel);
        }
    }


    @Override
    public JComponent getDiagnosticPanel() {
        return diagnosticPanel;
    }

}
