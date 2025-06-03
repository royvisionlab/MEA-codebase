package edu.ucsc.neurobiology.vision.snf;

import java.io.*;
import java.util.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.analysis.*;
import edu.ucsc.neurobiology.vision.anf.*;
import edu.ucsc.neurobiology.vision.electrodemap.*;
import edu.ucsc.neurobiology.vision.gui.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.math.fitting.*;
import edu.ucsc.neurobiology.vision.neuronviewer.*;
import edu.ucsc.neurobiology.vision.parameters.*;
import edu.ucsc.neurobiology.vision.plot.*;
import static edu.ucsc.neurobiology.vision.snf.SerialNeuronPanel.*;
import edu.ucsc.neurobiology.vision.util.*;
import edu.ucsc.neurobiology.vision.util.SwingWorker;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 * @author Matthew Grivich, The Salk Institute
 */
public class SerialMappingPanel extends JPanel
 implements Closable {

    public final static int nTracesToShow = 250;
    private int TRUNCATE_DATA = 50000;
    double reductionFactor = 1;

    private static VisionEventQueue eventQueue;
    private JScrollPane terminal;
    private static OutputStream outputStream;

    public double acfT1, acfT2;
    private int nlPointsEI, nrPointsEI;
    private int nlPoints, nrPoints, nPoints;
    private int nDimensions;
    private int maxClusters;
    private int minimizationInterval;
    private double minimizationError;
    private double minSpikeFindingSigma;
    private String title;
    private String rawFileName, originalRawFileName, masterModelFileName,
    slaveModelFileName, sigmasFileName, masterNeuronFileName, slaveNeuronFileName, eiFileName, datasetFolder,
    datasetName, masterDatasetFolder;
    private float[] sigmas;
    private int[][] amplitudeHistogram;
    private int nElectrodes;
    private ElectrodeMajorRawDataFile rawDataFile;
    private RawDataHeader header;
    private ElectrodeMap electrodeMap;
    private int nSamples;
    private int maxElectrodePatternSize;
    private int nSpikes;
    private int[] electrodes;
    private int[] extractionIDList;
    private int currentExtractionIndex;
    private int currentElectrode;
    private RawDataFile originalRawDataFile;
    private DoubleHistogram[] neuronAmplitudeHistograms;
    private float[][][][] averageEI;
    JInternalFrame mainFrame;
    DoubleHistogram[] autocorrelations;
    double[] contaminationIndex, nBadSpikes;
    Gaussian1DFunction[] amplitudeHistogramFit;
    String amplitudeHistogramsFileName;
    int[] neuronSpikeCount;
    JMenu viewMenu, optionsMenu;
    PlotPanel[] autocorrelationPanel, amplitudePanel;
    private float[][] projectedSpikes;
    private short[] rawData;
    PrintWriter html;
    boolean standalone = false;
    short[][][] spikes1;
    JMenuBar menuBar;
    float[] spikeAmplitudeList;
    double[] spikeTimeList, exactNeuronSpikeTimes, t0;
    int[] spikeID;
    ClusteringModelFile masterModelFile, slaveModelFile;
    ClusteringModelFile.Model currentModel;
    VisionHeader modelHeader;
    int currentCleaning;
    NeuronFile masterNeuronFile, slaveNeuronFile;
    PhysiologicalImagingFile imagingFile;
    IntegerList neuronsToRemoveList = new IntegerList();

    boolean canExitSNF = true;
    int maxEMSpikes;
    UserInterface userInterface;
    UserInterface[] userInterfaceList;
    JSplitPane splitPane1, splitPane2;
    int newID = 5000;

    // new
    IntegerParameter maxAmplitude;
    BooleanParameter showEI, autoOutline, autoLabels, showAllClusteringPlots,
    showGaussianFit, showRawData;
    EnumeratorParameter amplitudeScale;
    ParametersDialog optionsDialog;
    int nShownPlots = 3;
    boolean[] isIgnored;
    MyGlassPane myGlassPane;

    public final static int MANUAL = 0, AUTOMATIC = 1, RULES_AUTOMATIC = 2;
    private int automationLevel;


    public SerialMappingPanel(String masterDatasetFolder, String datasetFolder, 
            String originalRawFileName, int automationLevel) throws Exception {
        super(new GridLayout(1,1));
        this.masterDatasetFolder = masterDatasetFolder;
        this.datasetFolder = datasetFolder;
        this.originalRawFileName = originalRawFileName;
        this.automationLevel = automationLevel;
        
        datasetName = new File(datasetFolder).getName();
        String masterDatasetName = new File(masterDatasetFolder).getName();

        rawFileName = datasetFolder + File.separator + datasetName + ".rem";
        masterNeuronFileName = 
            masterDatasetFolder + File.separator + masterDatasetName + VisionParams.NEURON_FILE_EXTENSION;
        slaveNeuronFileName =
            datasetFolder + File.separator + datasetName +
            VisionParams.NEURON_FILE_EXTENSION;
        eiFileName =
            datasetFolder + File.separator + datasetName +
            VisionParams.EI_FILE_EXTENSION;
        sigmasFileName =
            masterDatasetFolder + File.separator + masterDatasetName +
            VisionParams.NOISE_FILE_EXTENSION;
        masterModelFileName =
            masterDatasetFolder + File.separator + masterDatasetName + ".model";
        slaveModelFileName = datasetFolder + File.separator + datasetName + ".model";

        title = datasetFolder;


        System.out.println("Opening raw data files...");
        rawDataFile = new ElectrodeMajorRawDataFile(rawFileName);
        originalRawDataFile = new RawDataFile(new File(originalRawFileName));
        header = rawDataFile.getHeader();
        this.nSamples = (int) (header.getNumberOfSamples() * reductionFactor);
        if (nSamples % 2 == 1) {
            nSamples--;
        }

        this.electrodeMap = ElectrodeMapFactory.getElectrodeMap(header.getArrayID());
        this.nElectrodes = electrodeMap.getNumberOfElectrodes();
        System.out.println("--- max Adjacents = " + maxElectrodePatternSize);

        // load the IgnoredElectrodes File
        isIgnored = SerialNeuronPanel.loadIgnoredElectrodesList(
                masterDatasetFolder, nElectrodes);

        System.out.println("Loading sigmas...");
        
        sigmas = SpikeFinding.getSigmas(sigmasFileName, nElectrodes);
        amplitudeHistogramsFileName =
            StringUtil.removeExtension(rawFileName) + ".sah";
        System.out.println("Loading amplitude histograms...");

        amplitudeHistogram = SerialNeuronPanel.loadAmplitudeHistograms(
                amplitudeHistogramsFileName);

        // open the Model file
        File f;
        System.out.println("Opening the Master Model file...");
        f = new File(masterModelFileName);
        if (f.exists() && f.isFile() && f.canRead()) {
            masterModelFile = new ClusteringModelFile(masterModelFileName);
            modelHeader = masterModelFile.getUserHeader();
            nlPoints = modelHeader.nlPoints;
            nrPoints = modelHeader.nrPoints;
            nlPointsEI = modelHeader.nlPointsEI;
            nrPointsEI = modelHeader.nrPointsEI;
            nPoints = nlPoints + nrPoints + 1;
            acfT1 = modelHeader.acfT1;
            acfT2 = modelHeader.acfT2;
            nDimensions = modelHeader.nDimensions;
            maxClusters = modelHeader.maxClusters + 2;
            //minimizationInterval = modelHeader.minimizationInterval;
            minimizationError = modelHeader.minimizationError;
            maxElectrodePatternSize = modelHeader.maxElectrodePatternSize;
            extractionIDList = masterModelFile.getIDList();
            minSpikeFindingSigma = modelHeader.minThreshold;
        } else {
            throw new IOException("The master model file does not exist");
        }

        System.out.println("Opening/Creating the Slave Model file...");
        f = new File(slaveModelFileName);
        if (f.exists() && f.isFile() && f.canRead()) {
            slaveModelFile = new ClusteringModelFile(slaveModelFileName);
        } else {
            slaveModelFile = new ClusteringModelFile(
                    slaveModelFileName, 5000, masterModelFile.getUserHeader());
        }

        
        System.out.println("Finding the maximum number of spikes...");
        amplitudeHistogramFit = new Gaussian1DFunction[maxClusters];
        neuronAmplitudeHistograms = new DoubleHistogram[maxClusters];
        for (int i = 0; i < neuronAmplitudeHistograms.length; i++) {
            neuronAmplitudeHistograms[i] = new DoubleHistogram(
                    "", 0, amplitudeHistogram[0].length, 1);
        }

        maxEMSpikes = Integer.MIN_VALUE;
        int maxSpikesElectrode = -1;
        for (int electrode = 1; electrode < nElectrodes; electrode++) {

            if (electrodeMap.isDisconnected(electrode) || isIgnored[electrode]) {
                continue;
            }

            int threshold = (int) Math.floor(
                    sigmas[electrode] * modelHeader.minThreshold) - 1;
            int n = 0;
            for (int i = threshold; i < amplitudeHistogram[electrode].length; i++) {
                n += amplitudeHistogram[electrode][i];
            }
            if (n > maxEMSpikes) {
                maxEMSpikes = n;
                maxSpikesElectrode = electrode;
            }
        }
        maxEMSpikes *= 1.2;

        maxEMSpikes *= reductionFactor;

        System.out.println("Maximum number of spikes (+20%): " + maxEMSpikes +
                ", on electrode " + maxSpikesElectrode);

        // memory demanding fields
        

        System.out.println("Allocating spike time memory storage...");
        exactNeuronSpikeTimes = new double[maxEMSpikes];
        spikeTimeList = new double[maxEMSpikes];
        spikeID = new int[maxEMSpikes];
        t0 = new double[maxEMSpikes];
        spikeAmplitudeList = new float[maxEMSpikes];

        projectedSpikes = new float[nDimensions][maxEMSpikes];

        System.out.println("Allocating spike waveform memory storage...");
        spikes1 = new short[maxElectrodePatternSize][nPoints + 2][maxEMSpikes];

        // less memory demanding
        System.out.println("Creating various memory structures...");
        autocorrelations = new DoubleHistogram[maxClusters];
        autocorrelationPanel = new PlotPanel[maxClusters];
        amplitudePanel = new PlotPanel[maxClusters];
        averageEI = new float[maxClusters][][][];
        contaminationIndex = new double[maxClusters];
        nBadSpikes = new double[maxClusters];
        neuronSpikeCount = new int[maxClusters];

        // create the neuron file
        System.out.println("Opening/Creating the neuron file...");
        f = new File(slaveNeuronFileName);

        if (f.exists() && f.isFile() && f.canRead()) {

            System.out.println("Open the existing neuron file.");
            slaveNeuronFile = new NeuronFile(slaveNeuronFileName);

        } else {
            // find the TTL signals
            rawData = new short[rawDataFile.getHeader().getNumberOfSamples()];
            rawDataFile.readRawData(0, rawData, rawData.length);

            int nTTL = SpikeFinder.findSpikes(
                    rawData, SerialNeuronPanel.ttlFindingThreshold,
                    spikeTimeList, spikeAmplitudeList, 0, 0);
            System.out.println("nTTL: " + nTTL);
            System.out.println("First TTL: " + spikeTimeList[0]);

            System.out.println("Create a new neuron file.");
            masterNeuronFile = new NeuronFile(masterNeuronFileName);
            VisionHeader h = masterNeuronFile.getHeader();
            System.out.println("version:" + h.version);
            slaveNeuronFile = new NeuronFile(slaveNeuronFileName, h, 5000, spikeTimeList, nTTL);


        }
        rawData = new short[nSamples];

        // open the EI file
        System.out.println("Opening/Creating the EI file...");
        f = new File(eiFileName);
        if (f.exists() && f.isFile() && f.canRead()) {
            System.out.println("The EI file exists. Continue writing in it.");
            imagingFile = new PhysiologicalImagingFile(eiFileName);
        } else {
            System.out.println("Create a new EI file.");
            imagingFile = new PhysiologicalImagingFile(
                    eiFileName, nlPointsEI, nrPointsEI, header.getArrayID());
        }

        // open the HTML report file
        System.out.println("Opening/Creating the report folder...");
        new File(datasetFolder + File.separator + "report").mkdir();
        String htmlFileName = datasetFolder + File.separator + "report" + File.separator +
        datasetName + ".html";

        File reportFolder = new File(datasetFolder + File.separator + "report");
        File htmlFile = new File(datasetFolder + File.separator + "report" +
                File.separator + datasetName + ".html");
        if (IOUtil.isValidFolder(reportFolder) &&
                IOUtil.isValidFile(htmlFile)) {
            System.out.println("Open the report file.");
            html = new PrintWriter(new FileWriter(htmlFileName, true), true);
        } else {
            System.out.println("Create a new report file.");
            html = new PrintWriter(new FileWriter(htmlFileName, true), true);
            html.println(
                    "<html><head><title>SNF Mapping Report</title></head><body>");
            html.println("<table cellspacing=1 cellpadding=1 border=1><tr>");

            html.println("<td> ID");
            html.println("<td> Salk ID");
            html.println("<td> Electrodes");
            html.println("<td> Color");
            html.println("<td> Spikes");
            html.println("<td> Contam");
            html.println("<td> Bad Spikes");
            html.println("<td> Amplitude");
            html.println("<td> Autocorrelation");
            html.println("<td> Amplitude Hist");
            html.println("<td> EI");
            html.println("<td> Clustering Plot");
        }

        System.out.println("Loading calculation status...");
        int _maxAmplitude = 1000;
        boolean _showEI = false;
        boolean _autoOutline = true;
        boolean _autoLabels = true;
        loadStatus();

        // create the options dialog
        maxAmplitude = new IntegerParameter(
                "Max Amplitude", "", "The spike amplitude scale", _maxAmplitude, 0, 2048);
        showEI = new BooleanParameter("Show EIs", "", "", _showEI);
        showRawData = new BooleanParameter("Show Raw Data", "", "", false);
        autoOutline = new BooleanParameter("Autocorrelation Outline", "", "",
                _autoOutline);
        autoLabels = new BooleanParameter("Autocorrelation Labels", "", "", _autoLabels);
        amplitudeScale = new EnumeratorParameter("Amplitude Scale", "", "");
        amplitudeScale.addChoice(1, "In ADC Counts");
        amplitudeScale.addChoice(2, "In Standard Deviations");

        showAllClusteringPlots = new BooleanParameter("Show All Clustering Plots", "", "", false);
        showGaussianFit = new BooleanParameter("Show Amplitude Histogram Fit", "", "", false);

        ParametersTable paramsTable = new ParametersTable();
        paramsTable.addParameter(maxAmplitude);
        paramsTable.addParameter(amplitudeScale);
        paramsTable.addParameter(showEI);
        paramsTable.addParameter(showRawData);
        paramsTable.addParameter(autoOutline);
        paramsTable.addParameter(autoLabels);
        paramsTable.addParameter(showGaussianFit);
        paramsTable.addParameter(showAllClusteringPlots);
        optionsDialog = new ParametersDialog(mainFrame, "Options", paramsTable);
        
        System.out.println("Initialization done");


        // start the neuron finding
        neuronMapping();

        if (!standalone) {
            Vision.getInstance().getCalculationManager().calculationDone();
        }
    }


    private void loadStatus() throws IOException {
        String statusFileName = datasetFolder + File.separator + datasetName + ".snms";
        File f = new File(statusFileName);
        if (f.exists() && f.canRead() && f.isFile()) {
            RandomAccessFile s = new RandomAccessFile(statusFileName, "rw");

            if (s.readInt() != SerialNeuronPanel.LATEST_VERSION_TAG) {
                throw new IOException("bad status file");
            }

            int nNeurons = s.readInt();
            for (int i = 0; i < nNeurons; i++) {
                neuronsToRemoveList.add(s.readInt());
            }

            currentExtractionIndex = s.readInt();
            currentCleaning = s.readInt();

            s.close();
        } else {
            currentExtractionIndex = 0;
            currentCleaning = 0;
        }

        System.out.println("Neurons To Remove: " + neuronsToRemoveList.size());
    }


    private void saveStatus() {
        String statusFileName =
            datasetFolder + File.separator + datasetName + ".snms";
        try {
            RandomAccessFile s = new RandomAccessFile(statusFileName, "rw");
            s.setLength(0);

            s.writeInt(SerialNeuronPanel.LATEST_VERSION_TAG);
            s.writeInt(neuronsToRemoveList.size());
            for (int i = 0; i < neuronsToRemoveList.size(); i++) {
                s.writeInt(neuronsToRemoveList.get(i));
            }
            s.writeInt(currentExtractionIndex);
            s.writeInt(currentCleaning);

            s.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    private void copySpikes(int electrodeIndex) {
        for (int sIndex = 0; sIndex < nSpikes; sIndex++) {
            // roundoff error can happen
            int time = (int) Math.round(spikeTimeList[sIndex]);

//			System.arraycopy(
//			rawData, time - nlPoints - 1,
//			spikes[sIndex][electrodeIndex], 0,
//			nPoints + 2);

            for (int i = 0; i < nPoints + 2; i++) {
                spikes1[electrodeIndex][i][sIndex] = rawData[time - nlPoints - 1 + i];
            }
        }
    }


    private void loadRawData(boolean doSpikeFinding) throws IOException {
        double t = 0;

        for (int electrodeIndex = 0; electrodeIndex < electrodes.length; electrodeIndex++) {
            double t1 = System.currentTimeMillis();
            rawDataFile.readRawData(electrodes[electrodeIndex], rawData, nSamples);
            t += (System.currentTimeMillis() - t1);

            if (electrodeIndex == 0 && doSpikeFinding) {
                nSpikes = SpikeFinder.findSpikes(
                        rawData, currentModel.threshold,
                        spikeTimeList, spikeAmplitudeList, 2 * nlPointsEI, 2 * nrPointsEI);
                System.out.println(
                        "Spikes: " + nSpikes + ", Threshold: " +
                        StringUtil.format(currentModel.threshold, 1));
            }
            copySpikes(electrodeIndex);
        }

        t /= 1000;
        double mBytesRead = nSamples * 1.5 * electrodes.length / 1024.0 / 1024.0;
        System.out.println(
                "Data reading speed: " + StringUtil.format(mBytesRead / t, 1) +
        " MBytes/sec.");
    }


    class ProjectingThread
    extends Thread {

        private int n1, n2;
        private UniformSpline spline;
        private boolean isWorking = true;
        double[] allignedSpike;


        public ProjectingThread(int n1, int n2) {
            this.n1 = n1;
            this.n2 = n2;
            spline = new UniformSpline(nPoints + 2);
            allignedSpike = new double[nPoints * maxElectrodePatternSize];
        }


        public void run() {
            final int eigenvectorLength = currentModel.eigenvectors[0].length;
            double[] spike = new double[nPoints + 2];

            for (int sIndex = n1; sIndex < n2; sIndex++) {
                int time = (int) Math.round(spikeTimeList[sIndex]); // trunaction error can happen
                for (int electrode = 0; electrode < electrodes.length; electrode++) {
                    if (time > nSamples - nrPoints - 1) {
                        continue;
                    }

//					spline.reSpline(spikes[sIndex][electrode]);

                    for (int i = 0; i < nPoints + 2; i++) {
                        spike[i] = spikes1[electrode][i][sIndex];
                    }
                    spline.reSpline(spike);

                    final int baseIndex = electrode * nPoints;
                    for (int i = 0; i < nPoints; i++) {
                        allignedSpike[baseIndex +
                                      i] = spline.getValueAt(t0[sIndex] + (i - nlPoints) - 1);
                    }
                }

                for (int d = 0; d < nDimensions; d++) {
                    projectedSpikes[d][sIndex] = 0;
                    for (int i = 0; i < eigenvectorLength; i++) {
                        projectedSpikes[d][sIndex] +=
                            allignedSpike[i] * currentModel.eigenvectors[d][i];
                    }
                }
            }

            isWorking = false;
        }
    }


    /**
     * This method is threaded but blocks until the projectiona are done.
     */
    private void projectSpikes() {
        System.out.print("Projecting spikes: ");
        double t1 = System.currentTimeMillis();

        ProjectingThread covThread1 = new ProjectingThread(0, nSpikes / 2);
        ProjectingThread covThread2 = new ProjectingThread(nSpikes / 2, nSpikes);
        covThread1.start();
        covThread2.start();

        while (covThread1.isWorking || covThread2.isWorking) {
            try {
                Thread.sleep(100);
            } catch (InterruptedException e) {
            }
        }

        double t2 = System.currentTimeMillis();
        System.out.println( (t2 - t1) / 1000 + " sec.");
    }


    private void map(int neuronIndex, int neuronID) throws IOException {
        int N = getNeuronSpikeTimes(neuronIndex);

        slaveNeuronFile.addExactNeuron(electrodes[0], neuronID, exactNeuronSpikeTimes, N);
        neuronsToRemoveList.add(neuronID);

        //write a report
        System.out.println("Neuron " + neuronID + " was extracted.");
        String ID = StringUtil.format(neuronID, 0, 3);
        String s = datasetFolder + File.separator + "report" + File.separator;

        html.println("<tr>");
        html.println("<td>" + neuronID);
        int[] idList = slaveNeuronFile.getIDList();
        int n = 0;
        for (int i = 0; i < idList.length; i++) {
            if (slaveNeuronFile.getElectrode(idList[i]) == currentElectrode) {
                n++;
            }
        }
        html.println("<td>" + StringUtil.format(currentElectrode, 0, 3) +
                (char) ('A' + n - 1));
        html.println("<td>" + StringUtil.toString(electrodes));

        html.println("<td>" + PlotUtil.getColorName(neuronIndex));
        html.println("<td>" + N);
        html.println("<td>" + StringUtil.format(contaminationIndex[neuronIndex], 3));
        html.println("<td>" + nBadSpikes[neuronIndex]);
        html.println("<td>" + "not calculated"); // amplitude
        html.println("<td><a href=\"" + ID + "-auto.png\">Autocorrelation</a>");
        html.println("<td><a href=\"" + ID + "-amp.png\">Amplitude Hist</a>");
        html.println("<td><a href=\"" + ID + "-ei.png\">EI</a>");
        html.println("<td>");
        html.println("<a href=\"" + ID + "-clust12.png\">1-2</a>");
        html.println("<a href=\"" + ID + "-clust13.png\">1-3</a>");
        html.println("<a href=\"" + ID + "-clust23.png\">2-3</a>");

        // save the plots
//		userInterface.clusteringPlots[0].autoscale();
        userInterface.clusteringPlots[0].saveAsPNG(s + ID + "-clust12.png" /*, imageSize,
                                                    imageSize*/);
//		userInterface.clusteringPlots[1].autoscale();
        userInterface.clusteringPlots[1].saveAsPNG(s + ID + "-clust13.png" /*, imageSize,
                                                    imageSize*/);
//		userInterface.clusteringPlots[2].autoscale();
        userInterface.clusteringPlots[2].saveAsPNG(s + ID + "-clust23.png" /*, imageSize,
                                                    imageSize*/);
        autocorrelationPanel[neuronIndex].autoscale();
        autocorrelationPanel[neuronIndex].saveAsPNG(s + ID + "-auto.png" /*, imageSize,
                                                    imageSize / 2*/);
        amplitudePanel[neuronIndex].autoscale();
        amplitudePanel[neuronIndex].saveAsPNG(s + ID + "-amp.png" /*, imageSize,
                          imageSize / 2*/);

//		STAViewer.makeImagingPanel(averageEI[neuron], electrodeMap).saveAsPNG(
//		new File(s + ID + "-ei.png"), imageSize, imageSize / 2);
    }


    private void nextElectrode(boolean block) {
        
        SwingWorker worker = new SwingWorker() {
            public Object construct() {
                if (userInterface != null) {
                    userInterface.lock();
                }

                // reached the end of mapping ?
                        if (currentExtractionIndex == extractionIDList.length) {
                            canExitSNF = false;
                            // do the last cleaning
                            saveStatus();
                            try {
                                SerialNeuronPanel.cleaning(
                                        neuronsToRemoveList, slaveNeuronFile, imagingFile, nlPointsEI,
                                        nrPointsEI, minimizationInterval, minimizationError,
                                        electrodeMap, rawData, rawDataFile, sigmas,
                                        minSpikeFindingSigma, spikeTimeList, spikeAmplitudeList,
                                        amplitudeHistogram, null, isIgnored, nSamples);
                                SerialNeuronPanel.saveAmplitudeHistograms(
                                        amplitudeHistogram,
                                        amplitudeHistogramsFileName);
                                canExitSNF = true;
                            } catch (IOException e) {
                                e.printStackTrace();
                            }

                            if (automationLevel == MANUAL) {
                                JOptionPane.showConfirmDialog(
                                        mainFrame,
                                        "The Mapping process is finised. The program will now exit",
                                        "Mapping Done", JOptionPane.OK_OPTION);
                            }

                            canExitSNF = true;
                            endCalculation();
                            return null;
                        } 

                        try {
                            currentModel = masterModelFile.getNeuronExtraction(extractionIDList[currentExtractionIndex]);
                            
                            // a cleaning is now required
                            if (currentCleaning != currentModel.cleaningLevel) {
                                canExitSNF = false;
                                saveStatus();
                                SerialNeuronPanel.cleaning(
                                        neuronsToRemoveList, slaveNeuronFile, imagingFile, nlPointsEI,
                                        nrPointsEI, minimizationInterval, minimizationError,
                                        electrodeMap, rawData, rawDataFile, sigmas,
                                        minSpikeFindingSigma, spikeTimeList, spikeAmplitudeList,
                                        amplitudeHistogram, null, isIgnored, nSamples);
                                SerialNeuronPanel.saveAmplitudeHistograms(
                                        amplitudeHistogram,
                                        amplitudeHistogramsFileName);
                                currentCleaning++;
                                neuronsToRemoveList.clear();
                                canExitSNF = true;
                            }
                        } catch (IOException e) {
                            e.printStackTrace();
                        }

                        System.out.println("\nExtraction " +
                                extractionIDList[currentExtractionIndex]);

                        electrodes = currentModel.electrodes;
                        currentElectrode = electrodes[0];
                        try {
                            loadRawData(true);
                        } catch (IOException e) {
                            e.printStackTrace();
                        }

                        try {
                            calculateExactSpikeTimes(
                                    nlPoints, nrPoints, nSpikes, nSamples,
                                    spikeTimeList, spikeAmplitudeList, t0,
                                    spikes1, minimizationInterval, minimizationError);
                        } catch (CannotEvaluateException e) {
                            e.printStackTrace();
                        }

                        projectSpikes();

                        return null;
            }


            public void finished() {
                setUI(currentModel);
                userInterface.dataUpdated();
                userInterface.updateNeuronsInfo();
                userInterface.updateView(true, true);

                String s = "Manual Mapper - " + title + ", El: " + currentElectrode +
                ", Spikes: " + nSpikes + ", Neurons to map: ";

                for (int i = 0; i < currentModel.neuronIndex.length; i++) {
                    s += PlotUtil.getColorName(currentModel.neuronIndex[i]) + ",";

                }
                mainFrame.setTitle(s);

                if (automationLevel == RULES_AUTOMATIC &&
                        userInterface instanceof EMUserInterface) {
                    int[] oldSpikeCount = new int[currentModel.neuronIndex.length];
                    int[] newSpikeCount = new int[currentModel.neuronIndex.length];
                    double[] oldContamination = new double[currentModel.neuronIndex.
                                                           length];
                    double[] newContamination = new double[currentModel.neuronIndex.
                                                           length];
                    for (int j = 0; j < currentModel.neuronIndex.length; j++) {
                        oldSpikeCount[j] = neuronSpikeCount[currentModel.neuronIndex[j]];
                        oldContamination[j] = contaminationIndex[currentModel.neuronIndex[
                                                                                          j]];
                    }

                    System.out.println("Direct Mapped Clusters:");
                    for (int i = 0; i < currentModel.neuronIndex.length; i++) {
                        System.out.println("Cluster: " +
                                PlotUtil.getColorName(currentModel.
                                        neuronIndex[i]));
                        System.out.println("Spikes: " +
                                neuronSpikeCount[currentModel.neuronIndex[i]]);
                        System.out.println("c: " +
                                StringUtil.format(contaminationIndex[
                                                                     currentModel.
                                                                     neuronIndex[i]], 2) +
                                                                     "(" +
                                                                     (int) nBadSpikes[currentModel.neuronIndex[i]] +
                        ")");
                    }

                    boolean reset = false;
                    //Do refitting
                    try {
                        ( (EMUserInterface) userInterface).refit();
                    } catch (FitFailedException e) {
                        System.out.println("EM: Could not fit beacuse: " +
                                e.getMessage());
                        reset = true;
                    }

                    if (reset == false) {
                        System.out.println("Refit Clusters:");
                        for (int i = 0; i < currentModel.neuronIndex.length; i++) {

                            System.out.println("Cluster: " +
                                    PlotUtil.getColorName(currentModel.
                                            neuronIndex[i]));
                            System.out.println("Spikes: " +
                                    neuronSpikeCount[currentModel.neuronIndex[
                                                                              i]]);
                            System.out.println("c: " +
                                    StringUtil.format(contaminationIndex[
                                                                         currentModel.
                                                                         neuronIndex[i]], 2) +
                                                                         "(" +
                                                                         (int) nBadSpikes[currentModel.neuronIndex[
                                                                                                                   i]] +
                            ")");
                        }

                        for (int j = 0; j < currentModel.neuronIndex.length; j++) {
                            newSpikeCount[j] = neuronSpikeCount[currentModel.neuronIndex[
                                                                                         j]];
                            newContamination[j] = contaminationIndex[currentModel.
                                                                     neuronIndex[
                                                                                 j]];
                        }
                        //if any neuron is bad, reset them all.
                        for (int j = 0; j < currentModel.neuronIndex.length; j++) {
                            double contaminationDifference = oldContamination[j] -
                            newContamination[j];
                            //if new contamination is higher, reset
                            if (contaminationDifference < -.01) {
                                reset = true;
                                //if new contamination is lower, do not reset
                            } else if (contaminationDifference > .01) {
                                ;
                            }
                            //if contamination is approximately equal, look at spike rate
                            else if (newSpikeCount[j] < oldSpikeCount[j]) {
                                reset = true;
                            }
                        }
                    }

                    if (reset) {
                        ( (EMUserInterface) userInterface).dataUpdated();
                        ( (EMUserInterface) userInterface).updateNeuronsInfo();

                    }
                    userInterface.updateView(true, true);
                    System.out.println("Kept Clusters:");
                    for (int i = 0; i < currentModel.neuronIndex.length; i++) {

                        System.out.println("Cluster: " +
                                PlotUtil.getColorName(currentModel.
                                        neuronIndex[i]));
                        System.out.println("Spikes: " +
                                neuronSpikeCount[currentModel.neuronIndex[i]]);
                        System.out.println("c: " +
                                StringUtil.format(contaminationIndex[
                                                                     currentModel.
                                                                     neuronIndex[i]], 2) +
                                                                     "(" +
                                                                     (int) nBadSpikes[currentModel.neuronIndex[i]] +
                        ")");

                    }

                }
                userInterface.unlock();
            }

        };

        if (block) {
            //  DO NOT MAKE NEW THREADS
            worker.construct();
            worker.finished();
        } else {
            worker.start();
        }
    }


    public static void calculateExactSpikeTimes(
            int nlPoints, int nrPoints, int nSpikes, int nSamples,
            double[] spikeTimeList, float[] spikeAmplitudeList, double[] t0,
            short[][][] spikes, int minimizationInterval, double minimizationError
    ) throws CannotEvaluateException {

        System.out.print("Calculating exact spike times: ");
        double t1 = System.currentTimeMillis();

        int nPoints = nlPoints + nrPoints + 1;
        double[] spike = new double[nPoints + 2];
        UniformSpline spline = new UniformSpline(nPoints + 2);

        for (int sIndex = 0; sIndex < nSpikes; sIndex++) {
            int time = (int) spikeTimeList[sIndex]; // no roundoff error can happen
            if (time - nlPoints < 0 || time + nrPoints >= nSamples) {
                t0[sIndex] = time;
                spikeTimeList[sIndex] = time;
                spikeAmplitudeList[sIndex] = time;
                continue;
            }

            // respline the spike on the main electrode
            for (int i = 0; i < nPoints + 2; i++) {
                spike[i] = spikes[0][i][sIndex];
            }
            spline.reSpline(spike);

            t0[sIndex] = FunctionMinimization.brentParabolic(
                    spline, nlPoints + 1 - minimizationInterval,
                    nlPoints + 1 + minimizationInterval, minimizationError);
            spikeTimeList[sIndex] = (double) (time - nlPoints - 1 + t0[sIndex]);
            spikeAmplitudeList[sIndex] = - (float) spline.getValueAt(t0[sIndex]);
            if (Math.abs(spikeTimeList[sIndex] - time) >= 1) {
                spikeTimeList[sIndex] = time;
            }
        }

        double t2 = System.currentTimeMillis();
        System.out.println( + (t2 - t1) / 1000 + " sec.");
    }


    private void neuronMapping() {
        // create the terminal
/*		JTextArea outputTextArea = new JTextArea();
        terminal = new JScrollPane(outputTextArea);
        System.setOut(new UIPrintStream(outputStream, outputTextArea));
        System.setErr(System.err);
*/
        eventQueue = new VisionEventQueue();
        // create user interfaces
        userInterfaceList = new UserInterface[3];
        userInterfaceList[0] = new EMUserInterface(1);
        userInterfaceList[1] = new EMUserInterface(2);
        userInterfaceList[2] = new BoxUserInterface();
        
        splitPane1 = new JSplitPane(JSplitPane.VERTICAL_SPLIT);
        splitPane2 = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
        this.add(splitPane2);
        setUI(0);
    //	splitPane1.setDividerLocation(800);
    //	splitPane2.setDividerLocation(400);


        
        mainFrame = Vision.getInstance().createFrame(this, null, makeMenuBar(), "Manual Neuron Finder");
//		mainFrame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);



        myGlassPane = new MyGlassPane(mainFrame);
        mainFrame.setGlassPane(myGlassPane);




        
        

        if (automationLevel == MANUAL) {
            nextElectrode(false);
        } else {
            canExitSNF = false;
            int[] extractionIDList = masterModelFile.getIDList();
            nextElectrode(true);
            userInterface.lock();
            for (int i = currentExtractionIndex; i < extractionIDList.length; i++) {
                userInterface.extractNeurons();
                currentExtractionIndex++;
                nextElectrode(true);
            }
            userInterface.unlock();
            canExitSNF = true;
        }

        // LOOPS FOREVER
    //	final SynchronizationObject syncObject = new SynchronizationObject();
    //	syncObject.setWorking();
    //	syncObject.waitUntilDone();
    }


    private void setUI(int i) {
        userInterface = userInterfaceList[i];
//		eventQueue.setKeyListener(userInterface);

        splitPane1.setLeftComponent(userInterface.rightPanel);
//		splitPane1.setRightComponent(terminal);
        splitPane2.setRightComponent(splitPane1);
        splitPane2.setLeftComponent(userInterface.leftPanel);

        splitPane1.setDividerLocation(800);
        splitPane2.setDividerLocation(400);
    }


    private void setUI(ClusteringModelFile.Model model) {
        if (model instanceof ClusteringModelFile.EMModel) {
            ClusteringModelFile.EMModel m = (ClusteringModelFile.EMModel) model;
            if (m.nGaussians == m.nClusters) {
                setUI(0);
            } else {
                setUI(1);
            }
        } else {
            setUI(2);
        }
    }
    
    private void endCalculation() {
   //     if (canClose()) {
            Vision.getInstance().removeFrame(mainFrame);
    //    }
        
    }


    public boolean canClose() {
        if (canExitSNF) {
            System.out.println("Saving amplitude histograms...");
            try {
                saveStatus();
                SerialNeuronPanel.saveAmplitudeHistograms(
                        amplitudeHistogram,
                        amplitudeHistogramsFileName);
                html.close();
                slaveNeuronFile.close();
                imagingFile.close();
                rawDataFile.close();
                originalRawDataFile.close();
                masterModelFile.close();
                slaveModelFile.close();
            } catch (IOException e) {
                e.printStackTrace();
            }

            Vision.getInstance().getCalculationManager().calculationDone();
            return true;
        } else {
            System.out.println(
                    "SNM cannot exit now, an unstoppable process is going on.");
            return false;
        }
    }


    private int getNeuronSpikeTimes(int neuron) {
        int n = 0;
        for (int sIndex = 0; sIndex < nSpikes; sIndex++) {
            if (spikeID[sIndex] == neuron) {
                exactNeuronSpikeTimes[n] = spikeTimeList[sIndex];
                n++;
            }
        }
        return n;
    }

    private JMenu[] makeMenuBar() {
        viewMenu = new JMenu("SNM: View");
        viewMenu.setMnemonic('V');
        
        optionsMenu = makeOptionsMenu();
        
        return new JMenu[] {viewMenu, optionsMenu};
    }

//	private JMenuBar makeMenuBar() {
//		final JMenuBar bar = new JMenuBar();
//		JMenuItem item;
//
//		viewMenu = new JMenu("View");
//		bar.add(viewMenu);
//
//		optionsMenu = makeOptionsMenu();
//		bar.add(optionsMenu);
//
//		return bar;
//	}


//	private void calculateEI(int neuron) {
//	int n = getNeuronSpikeTimes(neuron);
//	System.out.println(n);
//	averageEI[neuron] = WaveformCalculator.calculateEI(
//	exactNeuronSpikeTimes, n, 100, nlPointsEI, nrPointsEI, originalRawDataFile);
//	}

    public double[] getNeuronSpikeTimes1(int neuron) {
        DoubleList list = new DoubleList();
        for (int sIndex = 0; sIndex < nSpikes; sIndex++) {
            if (spikeID[sIndex] == neuron) {
                list.add(spikeTimeList[sIndex]);
            }
        }
        return list.toArray();
    }


    private void calculateEI(int neuron) throws IOException {
        double[] t = getNeuronSpikeTimes1(neuron);
        for (int i = 0; i < t.length; i++) {
            t[i] -= 50;
        }
        averageEI[neuron] = WaveformCalculator1.calculateEI(
                t, t.length, 100, nlPointsEI, nrPointsEI, originalRawDataFile);
        for (int e = 1; e < averageEI[neuron][0].length; e++) {
            double a = 0;
            for (int i = 0; i < 10; i++) {
                a += averageEI[neuron][0][e][i];
            }
            a /= 10;
            for (int i = 0; i < averageEI[neuron][0][e].length; i++) {
                averageEI[neuron][0][e][i] -= a;
            }
        }
    }


    private JMenu makeOptionsMenu() {
        JMenu menu = new JMenu("SNM: Options");
        menu.setMnemonic('O');
        JMenuItem item;

        item = new JMenuItem("Show Spike Time Dependence");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                double dt = DoubleDialog.showDoubleInputDialog(mainFrame, "Binning",
                        "Type the spike rate histogram binning in seconds", 1, 5e-5, 1e10);

                JFrame f = new JFrame("Spike Rate Time Dependence");
                f.getContentPane().setLayout(new GridLayout(0, 1));

                int[] ttl = null;
                try {
                    ttl = slaveNeuronFile.getTTLTimes();
                    System.out.println("nTTL " + ttl.length);
                } catch (IOException e) {
                    e.printStackTrace();
                }

                int nNeurons = userInterface.getNeuronsCount();
                for (int n = 0; n < nNeurons; n++) {
                    PlotPanel p = new PlotPanel();

                    // add the spike rate hist
                    DoubleHistogram h = new DoubleHistogram("", 0, nSamples / 20000.0, dt);
                    double[] t = getNeuronSpikeTimes1(n);
                    for (int i = 0; i < t.length; i++) {
                        h.fill(t[i] / 20000, 1);
                    }
                    h.scale(1.0 / dt);
                    p.addData(h, new HistogramStyle());

//					// add the TTLs
//ScatterPlot sp = new ScatterPlot();
//for (int sIndex = 0; sIndex < ttl.length; sIndex++) {
//					sp.add(ttl[sIndex] / 20000.0, 0);
//					}
//					ScatterPlotStyle style = new ScatterPlotStyle();
//					style.setSymbolType(SymbolType.VERTICAL_LINE);
//					style.setSymbolColor(Color.red);
//					style.setSymbolSize(50);
//					p.addData(sp, style);

                    // show the plots
                    p.setLabels("", PlotUtil.getColorName(n));
                    p.autoscale();
                    f.add(p);
                }

                f.setBounds(0, 0, 800, 600);
                f.setVisible(true);
            }
        });
        menu.add(item);

        item = new JMenuItem("Cross-Correlation");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ex) {
                JFrame f = new JFrame("Cross-correlation plots");
                f.getContentPane().setLayout(new GridLayout(0, 3));

                int n = userInterface.getNeuronsCount();
                for (int i = 0; i < n; i++) {
                    double[] t1 = getNeuronSpikeTimes1(i);
                    if (t1.length == 0) {
                        continue;
                    }
                    for (int j = i + 1; j < n; j++) {
                        double[] t2 = getNeuronSpikeTimes1(j);
                        if (t2.length == 0) {
                            continue;
                        }
                        DoubleHistogram h = CrossCorrelationCalculator.
                        getCrossCorrelationHistogram(
                                t1, t2, 1, 400);
                        h.scale(1.0 / Math.min(t1.length, t2.length));

                        // find the cross-correlation coefficient.
                        double[] cc = h.toArray();
                        double c = Double.NEGATIVE_INFINITY;
                        for (int k = 0; k < cc.length - 2 * coincidenceTime; k++) {
                            double v = 0;
                            for (int m = 0; m < 2 * coincidenceTime; m++) {
                                v += cc[k + m];
                            }
                            if (v > c) {
                                c = v;
                            }
                        }

                        PlotPanel p = new PlotPanel();
                        p.addToLegend(PlotUtil.getColorName(i) + " - " +
                                PlotUtil.getColorName(j) + " : " +
                                StringUtil.format(c * 100, 1) + "%");
                        p.addData(h, new HistogramStyle());
                        p.autoscale();
                        f.add(p);
                    }
                }

                f.setBounds(0, 0, 800, 600);
                f.setVisible(true);
            }
        });
        menu.add(item);

        item = new JMenuItem("Show Eigenvectors");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ex) {
                JDialog f = new JDialog((JFrame) null, "Eigenvectors", false);
                f.getContentPane().setLayout(new GridLayout(0, 1));

                ScatterPlotStyle style = new ScatterPlotStyle();
                style.setConnectingPoints(true);
                style.setConnectionPeriod(nPoints + 0);
                for (int d = 0; d < nDimensions; d++) {
                    ScatterPlot sp = new ScatterPlot();
                    double[] eigenvector = currentModel.eigenvectors[d];
                    for (int i = 0; i < eigenvector.length; i++) {
                        sp.add(i, eigenvector[i]);
                    }

                    PlotPanel p = new PlotPanel();
                    p.setAxisVisible(false);
                    p.setGridVisible(true);
                    p.addToLegend("" + (d + 1));
                    p.addData(sp, style);
                    p.setLabels("Time (samples)", "PC" + (d + 1));
                    p.autoscale();
                    p.replotAllData();
                    f.add(p);
                }

                f.setSize(600, 600);
                f.setLocationRelativeTo(mainFrame);
                f.setVisible(true);
            }
        });
        menu.add(item);

        item = new JMenuItem("PCs as a function of time");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ex) {
                JDialog f = new JDialog((JFrame) null, "PCs as a function of time", true);
                final int nNeurons = userInterface.getNeuronsCount();

                JPanel panel = new JPanel(new GridLayout(1, 0));
                panel.add(new JLabel("Dimension"));
                final JSpinner spinner = new JSpinner(
                        new SpinnerNumberModel(2, 1, nDimensions, 1));
                panel.add(spinner);
                f.add(panel, BorderLayout.NORTH);
                final PlotPanel p = new PlotPanel();
                f.add(p, BorderLayout.CENTER);

                spinner.addChangeListener(new ChangeListener() {
                    public void stateChanged(ChangeEvent e) {
                        p.removeAllData();
                        int d = ( (Integer) spinner.getValue()).intValue() - 1;

                        if (nNeurons == 0) {
                            ScatterPlot sp = new ScatterPlot();
                            for (int i = 0; i < nSpikes; i++) {
                                sp.add(spikeTimeList[i] / 20000.0, projectedSpikes[d][i]);
                            }

                            p.addData(sp, new ScatterPlotStyle());
                        } else {
                            for (int neuron = 0; neuron < nNeurons; neuron++) {
                                ScatterPlot sp = new ScatterPlot();
                                for (int i = 0; i < nSpikes; i++) {
                                    if (spikeID[i] == neuron) {
                                        sp.add(spikeTimeList[i] / 20000.0,
                                                projectedSpikes[d][i]);
                                    }
                                }

                                p.addData(sp, userInterface.clusteringPlotStyles[neuron]);
                            }
                        }

                        p.setLabels("Time (s)", "Principal Component " + (d + 1));
                        p.autoscale();
                        p.replotAllData();
                    }
                });

                spinner.setValue(new Integer(1));

                f.setSize(800, 400);
                f.setLocationRelativeTo(mainFrame);
                f.setVisible(true);
            }
        });
        menu.add(item);

        menu.add(new JSeparator());

        item = new JMenuItem("Options");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                optionsDialog.setVisible(true);
                if (showAllClusteringPlots.isSelected()) {
                    nShownPlots = dim.length;
                } else {
                    nShownPlots = 3;
                }
                userInterface.updateNeuronsInfo();
                userInterface.updateView(true, true);
            }
        });
        menu.add(item);

        return menu;
    }


    ActionListener viewListener = new ActionListener() {
        public void actionPerformed(ActionEvent e) {
            JComponent source = (JComponent) e.getSource();
            int neuron = Integer.parseInt(source.getName());
            try { //BEU
                calculateEI(neuron);
            } catch (IOException ex) {
                ex.printStackTrace();
            }

            new PhysiologicalImagePanel(
                    averageEI[neuron], sigmas, minSpikeFindingSigma, electrodeMap,
                    currentElectrode).showInAWindow(PlotUtil.getColorName(neuron));

//			PlotUtilities.showData(
//			"" + PlotUtilities.getColorName(neuron),
//			VisionUtilities.makeImagingPanel(averageEI[neuron], electrodeMap, -1));
        }
    };


    public abstract class UserInterface
    implements KeyListener {

        protected JPanel leftPanel, rightPanel, controlPanel, clusteringPanel;
        protected JButton resetBox, nextBox;
        protected PlotPanel[] clusteringPlots;
        protected JCheckBox grid = new JCheckBox("Grid", false);
        protected JCheckBox axes = new JCheckBox("Axes", false);
        public ScatterPlotStyle[] clusteringPlotStyles;
        PlotPanel[] rawDataPanels = new PlotPanel[maxClusters];


        private void addPlots(int plotIndex) {
            // add the plots
            for (int j = 0; j < maxClusters; j++) {
                ScatterPlotStyle style = new ScatterPlotStyle(
                        SymbolType.FILLED_SQUARE, 1,
                        PlotUtil.getColor(j), false, Color.black, 1);
                clusteringPlots[plotIndex].addData(
                        getClusteringPlot(j, dim[plotIndex][0], dim[plotIndex][1]), style);
            }
        }


        public UserInterface() {
            clusteringPlotStyles = new ScatterPlotStyle[maxClusters];
            for (int j = 0; j < maxClusters; j++) {
                clusteringPlotStyles[j] = new ScatterPlotStyle(
                        SymbolType.FILLED_SQUARE, 1,
                        PlotUtil.getColor(j), false, Color.black, 1);
            }

            // set up the clustering plots
            clusteringPlots = new PlotPanel[dim.length];
            for (int i = 0; i < dim.length; i++) {
                clusteringPlots[i] = makePlotPanel(i);
                clusteringPlots[i].setOverscaleFactor(1.2);
                clusteringPlots[i].setName("" + i);
                clusteringPlots[i].setLegendVisible(true);
                clusteringPlots[i].setAxisVisible(false);
                clusteringPlots[i].setLabels("Principal Feature " + (dim[i][0] + 1),
                        "Principal Feature " + (dim[i][1] + 1));
                // add the selection listeners
                clusteringPlots[i].addSelectionAction(new SelectionAction("Show Raw Data") {
                    public void selectionPerformed(JComponent source, Selection selection) {
                        canExitSNF = false;

                        PlotPanel p = (PlotPanel) source;
                        int neuron = p.getSelectedOption();
                        if (neuron >= getNeuronsCount()) {
                            neuron = -1;
                        }
                        showWaveformOverlap(Integer.parseInt(source.getName()), neuron,
                                selection.getSelection());
                        canExitSNF = true;
                    }
                });
            } // for i

            controlPanel = new JPanel();

            controlPanel.add(grid);
            grid.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    for (int i = 0; i < clusteringPlots.length; i++) {
                        clusteringPlots[i].setGridVisible(grid.isSelected());
                    }
                }
            });
            controlPanel.add(axes);
            axes.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    for (int i = 0; i < clusteringPlots.length; i++) {
                        clusteringPlots[i].setAxisVisible(axes.isSelected());
                    }
                }
            });

            resetBox = new JButton("Reset");
            resetBox.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    lock();

                    dataUpdated();
                    updateNeuronsInfo();
                    updateView(true, true);

                    unlock();
                }
            });
            controlPanel.add(resetBox);

            nextBox = new JButton("next");
            nextBox.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    canExitSNF = false;
                    lock();

                    extractNeurons();
                    currentExtractionIndex++;
                    nextElectrode(false);

                    canExitSNF = true;
                    unlock();
                }
            });
            controlPanel.add(nextBox);

            leftPanel = new JPanel(new BorderLayout());
            clusteringPanel = new JPanel(new GridLayout(0, 1));
            leftPanel.add(clusteringPanel, BorderLayout.CENTER);

            leftPanel.add(controlPanel, BorderLayout.SOUTH);
            rightPanel = new JPanel(new GridLayout(0, 2));
        }


        public void extractNeurons() {
            for (int i = 0; i < currentModel.neuronIndex.length; i++) {
                try {
                    map(currentModel.neuronIndex[i], currentModel.neuronID[i]);
                } catch (IOException ex) {
                    ex.printStackTrace();
                }
            }

            ClusteringModelFile.Model model = userInterface.getModel();
            model.extractionID = currentModel.extractionID;
            model.neuronIndex = currentModel.neuronIndex;
            model.neuronID = currentModel.neuronID;
            try {
                slaveModelFile.addExtraction(model);
            } catch (IOException ex) {
                ex.printStackTrace();
            }
            saveStatus();
        }


        public PlotPanel getWaveformOverlap(int neuronIndex) {
            BinaryRandom random = new BinaryRandom(nTracesToShow /
                    (double) neuronSpikeCount[neuronIndex]);
            ScatterPlot sp = new ScatterPlot("");

            for (int sIndex = 0; sIndex < nSpikes; sIndex++) {
                if (spikeID[sIndex] == neuronIndex && random.random()) {
                    for (int e = 0; e < electrodes.length; e++) {
                        for (int i = 0; i < nPoints + 2; i++) {
                            sp.add(e * (nPoints + 2) + i, spikes1[e][i][sIndex]);
                        }
                    }
                }
            }

            ScatterPlotStyle style1 = new ScatterPlotStyle();
            PlotPanel p = new PlotPanel();
            ScatterPlotStyle style = new ScatterPlotStyle();
            style.setConnectingPoints(true);
            style.setConnectionPeriod(nPoints + 2);
            style.setSymbolType(SymbolType.NONE);
            p.addData(sp, style);
            p.setAxisVisible(false);
//			p.setLabels("Time (samples)", "Amplitude (ADC)");
//p.getXAxis().setSecondaryLabels(1 / 20.0, 0, "(ms)");
p.autoscale();

return p;
        }


        public void showWaveformOverlap(int plotIndex, final int neuronIndex,
                SelectionRegion r) {
            int n = 0;
            double t = 0;

            IntegerList spikeIndexList = new IntegerList();
            for (int sIndex = 0; sIndex < nSpikes; sIndex++) {
                double x = projectedSpikes[dim[plotIndex][0]][sIndex];
                double y = projectedSpikes[dim[plotIndex][1]][sIndex];
                if (r.contains(x, y)) {
                    if (neuronIndex == -1) {
                        spikeIndexList.add(sIndex);
                    } else {
                        if (spikeID[sIndex] == neuronIndex) {
                            spikeIndexList.add(sIndex);
                        }
                    }
                }
            }

            BinaryRandom random = new BinaryRandom(nTracesToShow /
                    (double) spikeIndexList.size());
            ScatterPlot sp = new ScatterPlot("");

            for (int j = 0; j < spikeIndexList.size(); j++) {
                int sIndex = spikeIndexList.get(j);
                if (random.random()) {
                    for (int e = 0; e < electrodes.length; e++) {
                        for (int i = 0; i < nPoints + 2; i++) {
                            sp.add(e * (nPoints + 2) + i, spikes1[e][i][sIndex]);
                        }
                    }
                    if (n == 0) {
                        t = spikeTimeList[sIndex];
                    }
                    n++;
                }
            }

            if (n == 0) {
                System.out.println("Sorry, No traces to show");
                return;
            }

            ScatterPlotStyle style1 = new ScatterPlotStyle();
            PlotPanel p = new PlotPanel();
            ScatterPlotStyle style = new ScatterPlotStyle();
            style.setConnectingPoints(true);
            style.setConnectionPeriod(nPoints + 2);
            style.setSymbolType(SymbolType.NONE);
            p.addData(sp, style);
            p.setLabels("Time (samples)", "Amplitude (ADC)");
            p.getXAxis().setSecondaryLabels(1 / 20.0, 0, "(ms)");
            p.autoscale();
            if (n == 1) {
                PlotUtil.showData(
                        "Time: " + StringUtil.format(t, 1) + ", electrodes: " +
                        StringUtil.toString(electrodes), p);
            } else {
                PlotUtil.showData(
                        n + " traces on " + StringUtil.toString(electrodes), p);
            }
        }


        public abstract PlotPanel makePlotPanel(int i);


        public final void updateView(boolean complete, boolean autoscale) {
            int nGaussians = getNeuronsCount();

            viewMenu.removeAll();
            for (int neuron = 0; neuron < nGaussians; neuron++) {
                JMenuItem item = new JMenuItem("Show EI for " +
                        PlotUtil.getColorName(neuron));
                item.setName("" + neuron);
                item.addActionListener(viewListener);
                viewMenu.add(item);
            }

            if (complete) {
                // clear the view
                rightPanel.removeAll();
                clusteringPanel.removeAll();
                rightPanel.repaint();
                leftPanel.repaint();

                int n = 2;
                if (showEI.isSelected()) {
                    n++;
                }
                if (showRawData.isSelected()) {
                    n++;
                }
                rightPanel.setLayout(new GridLayout(0, n));
                Arrays.fill(rawDataPanels, null);

                for (int neuron = 0; neuron < nGaussians; neuron++) {
                    // add the autocorrelation plot
                    autocorrelationPanel[neuron] = new PlotPanel();
                    autocorrelationPanel[neuron].setLabels("Time between spikes (ms)",
                    "# of spike pairs");
                    autocorrelationPanel[neuron].setAxisVisible(autoLabels.getValue());
                    HistogramStyle autoStyle = new HistogramStyle();
                    autoStyle.setFillingTowers(true);
                    autoStyle.setFillingColor(PlotUtil.getColor(neuron));
                    if (autoOutline.getValue()) {
                        autoStyle.setOutlineColor(Color.black);
                    } else {
                        autoStyle.setOutlineColor(PlotUtil.getColor(neuron));
                    }
                    autocorrelationPanel[neuron].addData(autocorrelations[neuron],
                            autoStyle);
                    ArrayList<String> additionalLegend = new ArrayList<String>();
                    additionalLegend.add("Cluster: " + PlotUtil.getColorName(neuron));
                    additionalLegend.add("Spikes: " + neuronSpikeCount[neuron]);
                    additionalLegend.add(
                            "c: " + StringUtil.format(contaminationIndex[neuron], 2) +
                            "(" + (int) nBadSpikes[neuron] + ")");
                    autocorrelationPanel[neuron].setAdditionalLegend(additionalLegend);
                    autocorrelationPanel[neuron].setLegendVisible(true);
                    autocorrelationPanel[neuron].autoscale();
                    rightPanel.add(autocorrelationPanel[neuron]);

                    // add the EI
                    if (showEI.getValue()) {
                        try { //BEU
                            calculateEI(neuron);
                        } catch (IOException ex) {
                            ex.printStackTrace();
                        }
                        JComponent eiPanel = new PhysiologicalImagePanel(
                                averageEI[neuron], sigmas, minSpikeFindingSigma, electrodeMap,
                                currentElectrode);
                        rightPanel.add(eiPanel);
                    }

                    // add the Raw Data
                    if (showRawData.getValue()) {
                        rawDataPanels[neuron] = getWaveformOverlap(neuron);
                        rightPanel.add(rawDataPanels[neuron]);
                    }

                    // add the amplitude histogram
                    amplitudePanel[neuron] = new PlotPanel();
                    amplitudePanel[neuron].addData(
                            neuronAmplitudeHistograms[neuron], new HistogramStyle());

                    if (showGaussianFit.isSelected()) {
                        ArrayList<String> legend = new ArrayList<String>();
                        amplitudePanel[neuron].addData(
                                amplitudeHistogramFit[neuron],
                                new FunctionStyle("Amplitude Hist"));
                        legend.add("Mean: " + StringUtil.format(
                                amplitudeHistogramFit[neuron].getMean(), 1));
                        legend.add("RMS: " + StringUtil.format(
                                amplitudeHistogramFit[neuron].getSigma(), 1));
                        legend.add("Chi2: " + StringUtil.format(
                                amplitudeHistogramFit[neuron].getChiSquared(), 2));
                        amplitudePanel[neuron].setAdditionalLegend(legend);
                    }

                    amplitudePanel[neuron].autoscale();
                    amplitudePanel[neuron].setXRange(0, 1024);
                    if ( (int) amplitudeScale.getValue() == 1) {
                        amplitudePanel[neuron].setXRange(0, maxAmplitude.getValue());
                        amplitudePanel[neuron].setLabels("ADC Counts", "");
                    } else {
                        amplitudePanel[neuron].setXRange(0,
                                maxAmplitude.getValue() / sigmas[currentElectrode]);
                        amplitudePanel[neuron].setLabels("Standard Deviations", "");
                    }
                    rightPanel.add(amplitudePanel[neuron]);
                }

                // apply the same scale to all raw data plots
                if (showRawData.isSelected()) {
                    double yMin = Double.POSITIVE_INFINITY;
                    double yMax = Double.NEGATIVE_INFINITY;
                    for (int i = 0; i < nGaussians; i++) {
                        double[] r = rawDataPanels[i].getRange();
                        if (r[2] < yMin) {
                            yMin = r[2];
                        }
                        if (r[3] > yMax) {
                            yMax = r[3];
                        }
                    }
                    for (int i = 0; i < nGaussians; i++) {
                        rawDataPanels[i].setYRange(yMin, yMax);
                    }
                }

            } //if all

            if (showAllClusteringPlots.isSelected()) {
                clusteringPanel.setLayout(new GridLayout(0, 2));
            } else {
                clusteringPanel.setLayout(new GridLayout(0, 1));
            }

            for (int i = 0; i < clusteringPlots.length; i++) {
                clusteringPlots[i].clearOptions();
                for (int neuron = 0; neuron < nGaussians; neuron++) {
                    clusteringPlots[i].addOption(PlotUtil.getColorName(neuron));
                }
                clusteringPlots[i].addOption("All");
                clusteringPlots[i].setSelectedOption(nGaussians);
            }

            for (int i = 0; i < nShownPlots; i++) {
                clusteringPanel.add(clusteringPlots[i]);

                clusteringPlots[i].removeAllData();
                addPlots(i);

                PlotData[] contours = getClusterContours(i);
                PlotStyle style = getContoursStyle();
                for (int j = 0; j < contours.length; j++) {
                    clusteringPlots[i].addData(contours[j], style);
                }

                ArrayList<String> additionalLegend = new ArrayList<String>();
                additionalLegend.add("View: " + (dim[i][0] + 1) + " - " + (dim[i][1] + 1));
                clusteringPlots[i].setAdditionalLegend(additionalLegend);

                if (autoscale) {
                    clusteringPlots[i].autoscale();
                }
            }

            mainFrame.validate();
            mainFrame.repaint();
        }


        public abstract PlotData[] getClusterContours(int plotIndex);


        public abstract PlotStyle getContoursStyle();


        public abstract ClusteringModelFile.Model getModel();


        public abstract int getNeuronsCount();


        public abstract void sortSpikes();


        public final void updateNeuronsInfo() {
            int nNeurons = getNeuronsCount();

            sortSpikes();

            // extract the neuron times and compute the autocorrelations
            for (int neuron = 0; neuron < nNeurons; neuron++) {
                averageEI[neuron] = null;

                neuronAmplitudeHistograms[neuron].clear();
                int N = getNeuronSpikeTimes(neuron);
                neuronSpikeCount[neuron] = N;
                if ( (int) amplitudeScale.getValue() == 1) {
                    neuronAmplitudeHistograms[neuron].redefine(0, 2048);
                    for (int sIndex = 0; sIndex < nSpikes; sIndex++) {
                        if (spikeID[sIndex] == neuron) {
                            neuronAmplitudeHistograms[neuron].fill(
                                    spikeAmplitudeList[sIndex], 1);
                        }
                    }
                } else {
                    neuronAmplitudeHistograms[neuron].redefine(0,
                            2048 / sigmas[currentElectrode]);
                    for (int sIndex = 0; sIndex < nSpikes; sIndex++) {
                        if (spikeID[sIndex] == neuron) {
                            neuronAmplitudeHistograms[neuron].fill(
                                    spikeAmplitudeList[sIndex] / sigmas[currentElectrode], 1);
                        }
                    }
                }

                autocorrelations[neuron] = AutocorrelationCalculator.calculateISI(
                        exactNeuronSpikeTimes, N, 50, 0.05);
                nBadSpikes[neuron] = autocorrelations[neuron].count(acfT1, acfT2);
                double T = nSamples / VisionParams.SAMPLES_PER_MILLISECOND;
                double dT = acfT2 - acfT1;
                contaminationIndex[neuron] = nBadSpikes[neuron] * T / (dT * N * N);

                if (showGaussianFit.isSelected()) {
                    try {
                        amplitudeHistogramFit[neuron] =
                            neuronAmplitudeHistograms[neuron].fitToGaussian(10);
                    } catch (FitFailedException e) {
                        amplitudeHistogramFit[neuron] = new Gaussian1DFunction(0, 1, 1, 1);
                    }
                }
            }
        }


        public abstract void dataUpdated();


        public abstract ScatterPlotData getClusteringPlot(
                final int neuronID, final int d1, final int d2);


        public abstract void lock();


        public abstract void unlock();


        public void keyPressed(KeyEvent keyEvent) {
            switch (keyEvent.getKeyCode()) {
            case KeyEvent.VK_F5:

//				skipElectrode();
break;
            }
        }


        public void keyReleased(KeyEvent e) {
        }


        public void keyTyped(KeyEvent e) {
        }

    } // UserInterface


    public class EMUserInterface
    extends UserInterface {

        // EM related
        public final int binsPerDimension = 20;
        public static final int maxEMIterations = 100;

        JButton refitButton;
        private double[] min, max, binInterval;
        private ExpectationMaximization expectationMaximization;
        ClusteringModelFile.EMModel model;


        public EMUserInterface(int type) {
            super();

            min = new double[nDimensions];
            max = new double[nDimensions];
            binInterval = new double[nDimensions];

            switch (type) {
            case 1:
                expectationMaximization = new ExpectationMaximization1(
                        nDimensions, maxClusters, maxEMSpikes);
                break;

            case 2:
                expectationMaximization = new ExpectationMaximization2(
                        nDimensions, maxClusters, maxEMSpikes);
                break;
            }


            // clicks are used to add gaussians
            /*
                         for (int i = 0; i < nShownPlots; i++) {
                clusteringPlots[i].addClickListener(new ClickListener() {
                    public void clickPerformed(
                        Component source, MouseEvent e, double x, double y) {
                        if (e.getButton() != MouseEvent.BUTTON1) {
                            return;
                        }
                        int i = Integer.parseInt(source.getName());
                        // remove a gaussian
                        if ( (e.getModifiersEx() & MouseEvent.CTRL_DOWN_MASK) ==
                            MouseEvent.CTRL_DOWN_MASK) {
                            expectationMaximization.removeClosestGaussian(
                                x, y, dim[i][0], dim[i][1]);
                        } else {
                            // add a gaussian
                            double[] mean = new double[nDimensions];
                            double[] sigma = new double[nDimensions];
                            for (int d = 0; d < nDimensions; d++) {
                                mean[d] = (min[d] + max[d]) * 0.5;
                                sigma[d] = binInterval[d] * binInterval[d];
                            }
                            mean[dim[i][0]] = x;
                            mean[dim[i][1]] = y;
                            expectationMaximization.addGaussian(mean, sigma);
                        }
                        updateView(false);
                        refitButton.setEnabled(
                            expectationMaximization.getClustersCount() != 0);
                        nextBox.setEnabled(false);
                    }
                });
                         } // for i
             */

            refitButton = new JButton("Fit");
            refitButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    SwingWorker worker = new SwingWorker() {
                        public Object construct() {
                            try {
                                refit();
                            } catch (FitFailedException e) {
                                System.out.println("EM: Could not fit beacuse: " +
                                        e.getMessage());
                            }

                            return null;
                        }


                        public void finished() {
                            updateNeuronsInfo();
                            updateView(true, true);

                            unlock();
                        }
                    };
                    worker.start();

                }
            });
            controlPanel.add(refitButton);
        }


        public void refit() throws FitFailedException {

//			try {
//			expectationMaximization.setState(true, false, true);
//			expectationMaximization.fit(1e-6, 10,
//			SerialNeuronFinding.EMUserInterface.maxEMIterations);

            expectationMaximization.setState(true, true, true);
            expectationMaximization.fit(1e-6, 10,
                    SerialNeuronPanel.EMUserInterface.maxEMIterations);
//			} catch (FitFailedException e) {
//	System.out.println("EM: Could not fit beacuse: " +
//			e.getMessage());
//			}

            updateNeuronsInfo();

        }


        public ClusteringModelFile.Model getModel() {
            ClusteringModelFile.EMModel m = new ClusteringModelFile.EMModel();

            m.cleaningLevel = currentModel.cleaningLevel;
            m.electrodes = currentModel.electrodes;
            m.threshold = currentModel.threshold;
            m.nDimensions = currentModel.nDimensions;
            m.eigenvectors = currentModel.eigenvectors;

            m.nClusters = expectationMaximization.getClustersCount();
            m.nGaussians = expectationMaximization.getGaussiansCount();
            m.probability = new double[m.nGaussians];
            m.means = new double[m.nGaussians][];
            m.covariances = new double[m.nGaussians][];
            for (int i = 0; i < m.nGaussians; i++) {
                m.probability[i] = expectationMaximization.getGaussianProbability(i);
                m.means[i] = expectationMaximization.getMeans(i);
                m.covariances[i] = expectationMaximization.getCovariances(i);
            }

            return m;
        }


        public PlotPanel makePlotPanel(int i) {
            return new PlotPanel();
        }


        public PlotData[] getClusterContours(int i) {
            return expectationMaximization.getEllipses(dim[i][0], dim[i][1]);
        }


        public PlotStyle getContoursStyle() {
            return new FunctionStyle("Contour", Color.red, 2);
        }


        public void sortSpikes() {
            // assign an neuron ID for every spike
            for (int sIndex = 0; sIndex < nSpikes; sIndex++) {
                spikeID[sIndex] = expectationMaximization.getCluster(sIndex);
            }
        }


        public int getNeuronsCount() {
            return expectationMaximization.getClustersCount();
        }


        public void dataUpdated() {
            model = (ClusteringModelFile.EMModel) currentModel;

            for (int d = 0; d < nDimensions; d++) {
                min[d] = MathUtil.min(projectedSpikes[d]);
                max[d] = MathUtil.max(projectedSpikes[d]);
                binInterval[d] = (max[d] - min[d]) / binsPerDimension;
            }

            expectationMaximization.reset(projectedSpikes, TRUNCATE_DATA, nSpikes, model);
//			try {
//	expectationMaximization.setState(false, false, true);
//expectationMaximization.fit(
//			1e-6, 10, SerialNeuronFinding.EMUserInterface.maxEMIterations);
//			} catch (FitFailedException e) {
//			System.out.println("EM: Could not fit beacuse: " + e.getMessage());
//			System.out.println("The weights are not adjusted");
//			}
        }


        public ScatterPlotData getClusteringPlot(
                final int neuronID, final int d1, final int d2) {

            return new ScatterPlotData() {
                int j = 0;

                public void getDataPoint(int i, double[] x) {
                    if (expectationMaximization.isValid()) {
                        if (i == 0) {
                            j = 0;
                        }

                        if (j >= nSpikes) {
                            return;
                        }

                        while (spikeID[j] != neuronID) {
                            j++;
                            if (j >= nSpikes) {
                                return;
                            }
                        }

                        x[0] = projectedSpikes[d1][j];
                        x[1] = projectedSpikes[d2][j];

                        j++;
                    } else {
                        x[0] = projectedSpikes[d1][i];
                        x[1] = projectedSpikes[d2][i];
                    }
                }


                public int getPointCount() {
                    if (expectationMaximization.isValid()) {
                        return neuronSpikeCount[neuronID];
                    } else {
                        if (neuronID == 0) {
                            return nSpikes;
                        } else {
                            return 0;
                        }
                    }
                }


                public boolean hasXErrors() {
                    return false;
                }


                public boolean hasYErrors() {
                    return false;
                }


                public String getDescription() {
                    return "SP " + neuronSpikeCount[neuronID] + ", " + neuronID;
                }
            };
        }


        public void lock() {
            myGlassPane.setVisible(true);

            viewMenu.setEnabled(false);
            optionsMenu.setEnabled(false);
            refitButton.setEnabled(false);
            nextBox.setEnabled(false);
            resetBox.setEnabled(false);
        }


        public void unlock() {
            myGlassPane.setVisible(false);

            viewMenu.setEnabled(true);
            optionsMenu.setEnabled(true);
            refitButton.setEnabled(true);
            nextBox.setEnabled(true);
            resetBox.setEnabled(true);
        }


        public void keyPressed(KeyEvent keyEvent) {
            super.keyPressed(keyEvent);

            switch (keyEvent.getKeyCode()) {
            case KeyEvent.VK_F1:
                refitButton.doClick();
                break;

            case KeyEvent.VK_F2:
                nextBox.doClick();
                break;
            }
        }
    }


    public class BoxUserInterface
    extends UserInterface implements ChangeListener {


        public BoxUserInterface() {
            super();

            // add the selection listeners to add a clusters
            for (int i = 0; i < clusteringPlots.length; i++) {
                clusteringPlots[i].addPassiveSelectionChangeListener(this);
            }
        }


        public void stateChanged(ChangeEvent e) {
            updateNeuronsInfo();
            updateView(true, false);

            unlock();
        }


        public ClusteringModelFile.Model getModel() {
            ClusteringModelFile.ManualModel m = new ClusteringModelFile.ManualModel();

            m.cleaningLevel = currentModel.cleaningLevel;
            m.electrodes = currentModel.electrodes;
            m.threshold = currentModel.threshold;
            m.nDimensions = currentModel.nDimensions;
            m.eigenvectors = currentModel.eigenvectors;
            m.dimension1 = ( (ClusteringModelFile.ManualModel) currentModel).dimension1;
            m.dimension2 = ( (ClusteringModelFile.ManualModel) currentModel).dimension2;

            m.nClusters = getNeuronsCount();
            m.selections = new ArrayList<double[]>();
            for (int n = 0, j = 0; n < clusteringPlots.length; n++) {
                for (int i = 0; i < clusteringPlots[n].getPassiveSelectionsCount(); i++) {
                    m.selections.add(clusteringPlots[n].getPassiveSelection(i).toArray());
                }
            }

            return m;
        }


        public PlotPanel makePlotPanel(int i) {
            if (i == 0) {
                return new PlotPanel("", true, true, false, true);
            } else {
                return new PlotPanel("", true, false, false, false);
            }
        }


        public PlotData[] getClusterContours(int plotIndex) {
            return new ScatterPlot[0];
        }


        public PlotStyle getContoursStyle() {
            ScatterPlotStyle style = new ScatterPlotStyle();
            style.setSymbolType(SymbolType.NONE);
            style.setConnectingPoints(true);
            style.setConnectionLineThickness(1);
            return style;
        }


        /*
                public void sortSpikes() {
                    // assign an neuron ID for every spike
                    Arrays.fill(spikeID, 0);
                    for (int j = 0; j < boxPanel.getPassiveSelectionsCount(); j++) {
             SelectionRegion r = boxPanel.getPassiveSelection(j).getSelection();
                        for (int sIndex = 0; sIndex < nSpikes; sIndex++) {
                            double x = projectedSpikes[0][sIndex];
                            double y = projectedSpikes[1][sIndex];
                            if (r.contains(x, y)) {
                                spikeID[sIndex] = j + 1;
                            }
                        }
                    }
                }
         */
        public void sortSpikes() {
            // assign an neuron ID for every spike
            Arrays.fill(spikeID, 0);
            for (int n = 0, neuronIndex = 1; n < clusteringPlots.length; n++) {
                int d1 = dim[n][0];
                int d2 = dim[n][1];
                for (int j = 0; j < clusteringPlots[n].getPassiveSelectionsCount(); j++) {
                    SelectionRegion r = clusteringPlots[n].getPassiveSelection(j).
                    getSelection();
                    for (int sIndex = 0; sIndex < nSpikes; sIndex++) {
                        double x = projectedSpikes[d1][sIndex];
                        double y = projectedSpikes[d2][sIndex];
                        if (r.contains(x, y)) {
                            spikeID[sIndex] = neuronIndex;
                        }
                    }
                    neuronIndex++;
                }
            }
        }


//		public int getNeuronsCount() {
//		return boxPanel.getPassiveSelectionsCount() + 1;
//		}
        public int getNeuronsCount() {
            int N = 0;
            for (int n = 0; n < clusteringPlots.length; n++) {
                N += clusteringPlots[n].getPassiveSelectionsCount();
            }
            return N + 1;
        }


        public void dataUpdated() {
            ClusteringModelFile.ManualModel model =
                (ClusteringModelFile.ManualModel) currentModel;
            for (int n = 0; n < clusteringPlots.length; n++) {
                clusteringPlots[n].resetPassiveSelections();

            }

            for (int i = 0; i < model.selections.size(); i++) {
                int d1 = model.dimension1[i];
                int d2 = model.dimension2[i];
                PlotPanel boxPanel = clusteringPlots[getPlotByDimesnions(d1, d2)];

                Selection s = Selection.makeSelection(boxPanel, boxPanel.getXAxis(),
                        boxPanel.getYAxis(), (double[]) model.selections.get(i));
                boxPanel.addPassiveSelection(s);
            }
        }


        private int getPlotByDimesnions(int d1, int d2) {
            for (int i = 0; i < dim.length; i++) {
                if (dim[i][0] == d1 && dim[i][1] == d2) {
                    return i;
                }
            }
            return -1;
        }


        public ScatterPlotData getClusteringPlot(
                final int neuronID, final int d1, final int d2) {

            return new ScatterPlotData() {
                int j = 0;

                public void getDataPoint(int i, double[] x) {
                    if (getNeuronsCount() != 0) {
                        if (i == 0) {
                            j = 0;
                        }

                        if (j >= nSpikes) {
                            return;
                        }

                        while (spikeID[j] != neuronID) {
                            j++;
                            if (j >= nSpikes) {
                                return;
                            }
                        }

                        x[0] = projectedSpikes[d1][j];
                        x[1] = projectedSpikes[d2][j];

                        j++;
                    } else {
                        x[0] = projectedSpikes[d1][i];
                        x[1] = projectedSpikes[d2][i];
                    }
                }


                public int getPointCount() {
                    if (getNeuronsCount() != 0) {
                        return neuronSpikeCount[neuronID];
                    } else {
                        if (neuronID == 0) {
                            return nSpikes;
                        } else {
                            return 0;
                        }
                    }
                }


                public boolean hasXErrors() {
                    return false;
                }


                public boolean hasYErrors() {
                    return false;
                }


                public String getDescription() {
                    return "SP " + neuronSpikeCount[neuronID] + ", " + neuronID;
                }
            };
        }


        public void lock() {
            myGlassPane.setVisible(true);

            viewMenu.setEnabled(false);
            optionsMenu.setEnabled(false);
            nextBox.setEnabled(false);
            resetBox.setEnabled(false);
        }


        public void unlock() {
            myGlassPane.setVisible(false);

            viewMenu.setEnabled(true);
            optionsMenu.setEnabled(true);
            nextBox.setEnabled(true);
            resetBox.setEnabled(true);
        }


        public void keyPressed(KeyEvent keyEvent) {
            super.keyPressed(keyEvent);

            switch (keyEvent.getKeyCode()) {
            case KeyEvent.VK_F2:
                nextBox.doClick();
                break;
            }
        }

    }




}
