package edu.ucsc.neurobiology.vision.analysis;

import java.io.File;
import java.util.HashMap;

import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;

import com.nativelibs4java.opencl.*;
import edu.ucsc.neurobiology.vision.Vision;
import edu.ucsc.neurobiology.vision.analysis.SpikeBlockQueueCL.SpikeBlockCL;
import edu.ucsc.neurobiology.vision.calculations.AbstractCalculation;
import edu.ucsc.neurobiology.vision.io.AsyncBufferCL;
import edu.ucsc.neurobiology.vision.io.NeuronFile;
import edu.ucsc.neurobiology.vision.io.RawDataHeader;
import edu.ucsc.neurobiology.vision.io.SampleBlockQueueCL;
import edu.ucsc.neurobiology.vision.parameters.TableLayout;

public class ImagingCL extends AbstractCalculation {
    private String rawDataFileName, neuronFileName;
    private int nBuffers, bufferSize;
    private int nlPoints, nrPoints;
    private int spikesToUse;
    private float alpha;

    private RawDataHeader header;
    
    private CLContext context;
    private CLQueue clQueue;
    private SampleBlockQueueCL sampleQueue;
    private SpikeBlockQueueCL spikeQueue;
    private EICalculatorCL eiCalculator;
    
    private JPanel diagnosticPanel;
    private JLabel samplesReadLabel;

    public ImagingCL() {
        if (Vision.isGUIBased())
            diagnosticPanel = new JPanel(new TableLayout(0,0));
    }
    
    
    @Override
    public void setParameters(HashMap<String, String> parameters) throws Exception {
        context = JavaCL.createBestContext(CLPlatform.DeviceFeature.GPU);
        clQueue = context.createDefaultQueue();
        
        rawDataFileName = parameters.get("Raw Data File");
        nBuffers = Integer.parseInt(parameters.get("Number of Sample Buffers"));
        bufferSize = Integer.parseInt(parameters.get("Sample Buffer Length"));
        
        String datasetFolder = parameters.get("Dataset Folder");
        neuronFileName = datasetFolder + File.separator + new File(datasetFolder).getName() + ".neurons";
        
        nlPoints = Integer.parseInt(parameters.get("Left Samples"));
        nrPoints = Integer.parseInt(parameters.get("Right Samples"));
        spikesToUse = Integer.parseInt(parameters.get("Spikes To Average"));
        
        float timeConstant = (float) Double.parseDouble(parameters.get("Mean Time Constant"));
        alpha = (1f / (timeConstant * 20000f));
        
        if (Vision.isGUIBased()) {
            diagnosticPanel.add(new JLabel("Samples read:"));
            samplesReadLabel = new JLabel("0", JLabel.RIGHT);
            diagnosticPanel.add(samplesReadLabel);
        }
    }

    
    @Override
    public void startCalculation() throws Exception {
        Vision app = Vision.getInstance();

        app.sendMessage("Preparing...");
        NeuronFile neuronFile = new NeuronFile(neuronFileName);
        sampleQueue = SampleBlockQueueCL.create(context, clQueue, rawDataFileName, nBuffers, bufferSize);
        spikeQueue = SpikeBlockQueueCL.create(context, clQueue, nBuffers, neuronFile, spikesToUse, nlPoints, nrPoints, bufferSize);
        int maxSpikesPerCellPerBlock = spikeQueue.getMaxSpikesPerCellPerBlock();
        header = sampleQueue.getHeader();
        eiCalculator = new EICalculatorCL(neuronFile, header, nlPoints, nrPoints, bufferSize, maxSpikesPerCellPerBlock, 
                context, clQueue);
        
        app.sendMessage("Calculating EIs...");
        app.startProgressBar();
        sampleQueue.start();
        spikeQueue.start();

        AsyncBufferCL<Short> sampleBuffer;
        SpikeBlockCL spikeBlock;
        int sampleIndex = 0;
        while (true) {
            sampleBuffer = sampleQueue.take();
            if (sampleBuffer.getBuffer() == null) break;
            spikeBlock = spikeQueue.take();

            eiCalculator.processSamples(sampleIndex, sampleBuffer, spikeBlock);
            sampleBuffer.close();
            spikeBlock.giveBack();
            sampleIndex += bufferSize;

            app.setProgress( (int) (100.0 * sampleIndex / header.getNumberOfSamples()));
            if (Vision.isGUIBased()) samplesReadLabel.setText(Integer.toString(sampleIndex));
        }
        System.out.println("\n" + sampleIndex / header.getSamplingFrequency() + " s of data processed");
        
        eiCalculator.finishCalculation();
        app.endProgressBar();
        calculationDone();
    }
    
    
    @Override
    public JComponent getDiagnosticPanel() {
        return diagnosticPanel;
    }
}
