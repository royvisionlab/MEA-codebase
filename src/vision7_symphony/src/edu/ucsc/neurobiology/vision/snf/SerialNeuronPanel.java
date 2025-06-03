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
import static edu.ucsc.neurobiology.vision.util.VisionParams.*;
import edu.ucsc.neurobiology.vision.util.*;
import edu.ucsc.neurobiology.vision.util.SwingWorker;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 * @author Matthew Grivich, The Salk Institute
 */
public class SerialNeuronPanel
extends JPanel implements Closable {

    enum CalculateMode {
        SWITCH_ELECTRODE, RELOAD_DATA, NEW_SPIKE_FINDING, RECLUSTER, ISOLATE_CLUSTER,
    }


    // old versions: none
    public final static int LATEST_VERSION_TAG = 0xABCDEF00;

    private static VisionEventQueue eventQueue;
//	private static OutputStream outputStream;
//	private JScrollPane terminal;

    public final static int nTracesToShow = 250;
    private int TRUNCATE_DATA = 50000;
    public ElectrodeUsage electrodeUsage;
    private final int maxElectrodePatternSize = 14;
    private final static double allowedCleaningError = 3;
    public final static double ttlFindingThreshold = 100;

    public static final int[][] dim = { {0, 1}
    , {0, 2}
    , {1, 2}
    , {0, 3}
    , {1, 3}
    , {2, 3}
    , {0, 4}
    , {1, 4}
    , {2, 4}
    , {3, 4}
    };
    int nShownPlots = 3;

    private final int maxClusters = 16;
    public static final int imageSize = 400;
    public static final int nEISpikes = 1000;
    public double minSpikeFindingSigma = -1;
    private final int nlPoints = 10, nrPoints = 15;
    private final int nPoints = nlPoints + nrPoints + 1;
    private final static int nlPointsEI = 1 * 20, nrPointsEI = 3 * 20;
    public final static int ignoreAtBeginning = Math.max(2 * nlPointsEI, 100);
    private final int nDimensions = 5;
    private final int minimizationInterval = 1;
    private final double minimizationError = 1e-6;
    public static final double coincidenceTime = 10;
    public static final double maxCorrelation = 0.25;

    // common params
    private String rawFileName, originalRawFileName, datasetFolder, sigmasFileName,
    datasetName;
    private float[] sigmas;
    private int[][] amplitudeHistogram;
    private int nElectrodes;
    private ElectrodeMajorRawDataFile rawDataFile;
    private RawDataHeader header;
    private ElectrodeMap electrodeMap;
    private int nSamples;
    private int nSpikes;
    private int[] electrodes;
    private int currentElectrode;
    private RawDataFile originalRawDataFile;
    private DoubleHistogram[] neuronAmplitudeHistograms;
    private float[][][][] averageEI;
    int _maxSpikes;
    DoubleHistogram[] autocorrelations;
    double[] contaminationIndex, nBadSpikes;
    Gaussian1DFunction[] amplitudeHistogramFit;
    NeuronFile neuronFile;
    boolean[] skipElectrode, isIgnored, skipElectrodeForever;
    String amplitudeHistogramsFileName;
    int[] neuronSpikeCount;
    PlotPanel[] autocorrelationPanel, amplitudePanel;
    private float[][] projectedSpikes;
    private short[] rawData;
    private double floatingSigma;
    PrintWriter html;
    PCA pca;
    PhysiologicalImagingFile imagingFile;
    IntegerList neuronsToRemoveList = new IntegerList();
    short[][][] spikes1;
    JMenu viewMenu;
    JMenu[] menuList;
    float[] spikeAmplitudeList;
    double[] spikeTimeList; // the list of all spike times
    double[] exactNeuronSpikeTimes, t0;
    int[] spikeID;
    ClusteringModelFile modelFile;
    boolean canExitSNF;

    // the number of times the cleaning was already done
    // 0 - means the cleaning was never yet done
    int cleaningLevel;

    JInternalFrame mainFrame;
    UserInterface userInterface;
    JSplitPane splitPane1, splitPane2;
    UserInterface[] userInterfaceList;
    int currentExtractionID;
    ElectrodeDialog elecrodeDialog;
    ParametersDialog optionsDialog;
    IntegerParameter maxAmplitude, nSamplesParameter /*, truncateDataParameter*/;
    BooleanParameter showEI, autoOutline, autoLabels, showAllClusteringPlots,
    showGaussianFit, showRawData;
    EnumeratorParameter amplitudeScale;

    int _maxAmplitude;
    boolean _showEI;
    boolean _autoOutline, _autoLabels;

    int nSamplesToAnalyze;
    MyGlassPane myGlassPane;


    public static final boolean[] loadIgnoredElectrodesList(
            String datasetFolder, int nElectrodes) {

        boolean[] isIgnored = new boolean[nElectrodes];
        String ignoreFileName = datasetFolder + File.separator + "ignore-electrodes.txt";
        File f = new File(ignoreFileName);

        if (f.exists() && f.isFile() && f.canRead()) {
            try {
                LineNumberReader r = new LineNumberReader(new FileReader(ignoreFileName));
                System.out.print("Ignored electrodes: ");
                while (true) {
                    String line = r.readLine();
                    if (line == null) {
                        break;
                    }
                    line = line.trim();
                    if (line.length() == 0) {
                        break;
                    }

                    int electrode = Integer.parseInt(line);
                    isIgnored[electrode] = true;
                    System.out.print(electrode + ", ");
                }
                System.out.println();
                r.close();
            } catch (IOException e) {
                Vision.reportFatalException(
                        "The ignore-electrodes.txt file cannot be opened", e);
            }
        }

        return isIgnored;
    }


    public SerialNeuronPanel(String datasetFolder, String originalRawFileName) {
        super(new GridLayout(1, 1));

        this.datasetFolder = datasetFolder;
        this.originalRawFileName = originalRawFileName;
        String datasetName = new File(datasetFolder).getName();
        rawFileName = datasetFolder + File.separator + datasetName + ".rem";
        sigmasFileName = datasetFolder + File.separator + datasetName + ".noise";
        /*
                // redirect the output to a folder-specific file
                try {
                    outputStream = new FileOutputStream(
                        datasetFolder + File.separator + "output.txt", true);
                    PrintStream output = new PrintStream(outputStream, true);
                    System.setErr(output);
                    System.setOut(output);

         System.out.println("\n======== SNF Started ======= " + new Date() + "\n");
                } catch (IOException e) {
                    VisionUtilities.reportFatalException(
                        "The output.txt file cannot be opened, created or written", e);
                }
         */


        System.out.println("Opening raw data files...");

        try {
            rawDataFile = new ElectrodeMajorRawDataFile(rawFileName);
        } catch (IOException e) {
            Vision.reportFatalException("The .rem file cannot be opened.", e);
        }
        try {
            originalRawDataFile = new RawDataFile(new File(originalRawFileName));
        } catch (IOException e) {
            Vision.reportFatalException("The .bin file cannot be opened.", e);
        }
        header = rawDataFile.getHeader();
        this.nSamples = header.getNumberOfSamples();
        this.nSamplesToAnalyze = nSamples;

        System.out.println("Creating electrode adjacency tables...");
        this.electrodeMap = ElectrodeMapFactory.getElectrodeMap(header.getArrayID());
        this.nElectrodes = electrodeMap.getNumberOfElectrodes();
        skipElectrode = new boolean[nElectrodes];
        skipElectrodeForever = new boolean[nElectrodes];

        System.out.println("Loading sigmas...");
        try {
            sigmas = SpikeFinding.getSigmas(sigmasFileName, nElectrodes);
        } catch (IOException e) {
            Vision.reportFatalException("The .noise file cannot be oppened", e);
        }
        amplitudeHistogramsFileName =
            StringUtil.removeExtension(rawFileName) + ".sah";
        System.out.println("Loading amplitude histograms...");
        try {
            amplitudeHistogram = loadAmplitudeHistograms(amplitudeHistogramsFileName);
        } catch (IOException e) {
            Vision.reportFatalException(
                    "The amplitude histograms (.sah) file cannot be oppened", e);
        }
        elecrodeDialog = new ElectrodeDialog();

        System.out.println("Finding the maximum number of spikes...");
        amplitudeHistogramFit = new Gaussian1DFunction[maxClusters];
        neuronAmplitudeHistograms = new DoubleHistogram[maxClusters];
        for (int i = 0; i < neuronAmplitudeHistograms.length; i++) {
            neuronAmplitudeHistograms[i] = new DoubleHistogram(
                    "", 0, amplitudeHistogram[0].length, 1);
        }

        // memory demanding fields
        System.out.println("Allocating raw data memory storage...");
        rawData = new short[nSamples];

        int response = createOrOpenFiles();

        // load the IgnoredElectrodes File
        isIgnored = loadIgnoredElectrodesList(datasetFolder, nElectrodes);

        _maxSpikes = Integer.MIN_VALUE;
        int maxSpikesElectrode = -1;
        for (int electrode = 1; electrode < nElectrodes; electrode++) {
            if (electrodeMap.isDisconnected(electrode) || isIgnored[electrode]) {
                continue;
            }

            int threshold = (int) Math.floor(sigmas[electrode] * minSpikeFindingSigma) -
            1;
            int n = 0;
            for (int i = threshold; i < amplitudeHistogram[electrode].length; i++) {
                n += amplitudeHistogram[electrode][i];
            }
            if (n > _maxSpikes) {
                _maxSpikes = n;
                maxSpikesElectrode = electrode;
            }
        }
        _maxSpikes *= 1.2;
        if (_maxSpikes < 400000) {
            _maxSpikes = 400000;
        }
        System.out.println("Maximum number of spikes (+20%): " + _maxSpikes +
                ", on electrode " + maxSpikesElectrode);

        System.out.println("Allocating memory for a max of " + _maxSpikes +
                " spikes on electrode " + maxSpikesElectrode + "...");
        exactNeuronSpikeTimes = new double[_maxSpikes];
        spikeTimeList = new double[_maxSpikes];
        spikeAmplitudeList = new float[_maxSpikes];
        spikeID = new int[_maxSpikes];
        t0 = new double[_maxSpikes];
        projectedSpikes = new float[nDimensions][_maxSpikes];
        spikes1 = new short[maxElectrodePatternSize][nPoints + 2][_maxSpikes];

        // memory demanding
        System.out.println("Creating various memory structures...");
        autocorrelations = new DoubleHistogram[maxClusters];
        autocorrelationPanel = new PlotPanel[maxClusters];
        amplitudePanel = new PlotPanel[maxClusters];
        averageEI = new float[maxClusters][][][];
        contaminationIndex = new double[maxClusters];
        nBadSpikes = new double[maxClusters];
        neuronSpikeCount = new int[maxClusters];

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

//		truncateDataParameter = new IntegerParameter(
//		"Spikes used for COV and EM", "", "", TRUNCATE_DATA, 1, Integer.MAX_VALUE);
        nSamplesParameter = new IntegerParameter(
                "Samples to analyze", "", "", nSamples, 1, nSamples);

        System.out.println("Initialization done");

        // start the neuron finding
        manualNeuronFinding();
    }


    private void initializeStatus() {
        floatingSigma = minSpikeFindingSigma;
        for (int i = 0; i < nElectrodes; i++) {
            skipElectrode[i] = false;
        }
        cleaningLevel = 0;
//		currentExtractionID = 1;
        electrodeUsage = ElectrodeUsage.SEVEN_ELECTRODES;
        _maxAmplitude = 1000;
        _showEI = true;
        _autoOutline = true;
        _autoLabels = true;
    }


    private int createOrOpenFiles() {
        datasetName = new File(datasetFolder).getName();

        File neuronFileName = new File(datasetFolder + File.separator + datasetName +
                VisionParams.NEURON_FILE_EXTENSION);
        File reportFolder = new File(datasetFolder + File.separator + "report");
        File htmlFileName = new File(datasetFolder + File.separator + "report" +
                File.separator + datasetName + ".html");
        File eiFileName = new File(datasetFolder + File.separator + datasetName +
                VisionParams.EI_FILE_EXTENSION);
        File modelFileName = new File(datasetFolder + File.separator + datasetName +
                ".model");

        System.out.println("Loading calculation status...");
        File statusFile = new File(datasetFolder + File.separator + datasetName + ".snfs");
        if (IOUtil.isValidFile(statusFile)) {
            try {
                loadStatus(statusFile);
            } catch (IOException e) {
                JOptionPane.showMessageDialog(
                        mainFrame,
                        "The status file is not compatible with this version of SNF." +
                        "\nThis version of SNF cannot be used with this dataset.",
                        "Incompatible status file", JOptionPane.ERROR_MESSAGE);
                return -1;
            }
        } else {
            initializeStatus();
        }

        // loading or creating all the other files
        if (IOUtil.isValidFile(neuronFileName) &&
                IOUtil.isValidFile(htmlFileName) &&
                IOUtil.isValidFile(eiFileName) &&
                IOUtil.isValidFile(modelFileName) &&
                IOUtil.isValidFolder(reportFolder)) {
            System.out.println("All output files exist, open them.");

            System.out.println("Opening/Creating the neuron file...");
            try {
                neuronFile = new NeuronFile(neuronFileName.getAbsolutePath());
            } catch (IOException e) {
                Vision.reportFatalException(
                        "The neurons (.neurons) file cannot be oppened", e);
            }

            System.out.println("Opening/Creating the report folder...");
            reportFolder.mkdir();
            try {
                html = new PrintWriter(new FileWriter(htmlFileName, true), true);
            } catch (IOException e) {
                Vision.reportFatalException(
                        "The report (.html) file cannot be oppened", e);
            }

            System.out.println("Opening/Creating the EI file...");
            try {
                imagingFile = new PhysiologicalImagingFile(eiFileName.getAbsolutePath());
            } catch (IOException e) {
                Vision.reportFatalException(
                        "The electrode immaging (.ei) file cannot be oppened", e);
            }

            System.out.println("Opening/Creating the Model file...");
            try {
                modelFile = new ClusteringModelFile(modelFileName.getAbsolutePath());
            } catch (IOException e) {
                Vision.reportFatalException(
                        "The model file (.model) file cannot be oppened", e);
            }

            System.out.println("Loading calculation status...");
        } else if (!IOUtil.isValidFile(neuronFileName) &&
                !IOUtil.isValidFile(htmlFileName) &&
                !IOUtil.isValidFile(eiFileName) &&
                !IOUtil.isValidFile(modelFileName) &&
                !IOUtil.isValidFolder(reportFolder)) {

            System.out.println("None of the output files exist, create them.");

            // Since you are starting SNF for the first time you can choose the value

            minSpikeFindingSigma = DoubleDialog.showDoubleInputDialog(
                    mainFrame, "Threshold",
                    "Input the Threshold Value In Standard Deviations:", 4, 1, 1e10);

            System.out.println("Opening/Creating the neuron file...");
            try {
                rawDataFile.readRawData(0, rawData,
                        rawDataFile.getHeader().getNumberOfSamples());
            } catch (IOException e) {
                Vision.reportFatalException("The .rem file cannot be read", e);
            }
            spikeTimeList = new double[100000];
            spikeAmplitudeList = new float[100000];
            int nTTL = SpikeFinder.findSpikes(
                    rawData, ttlFindingThreshold, spikeTimeList, spikeAmplitudeList, 0, 0);
//			System.out.println("nTTL: " + nTTL);
//			for (int i = 0; i < nTTL; i++) {
//			System.out.println(spikeTimeList[i]);
//			}
//			System.out.println("end TTL");

            // create the Header
            VisionHeader visionHeader = new VisionHeader();
            visionHeader.version = NeuronFile.DOUBLE_VERSION;
            visionHeader.arrayID = header.getArrayID();
            visionHeader.nSamples = header.getNumberOfSamples();
            visionHeader.samplingFrequency = header.getSamplingFrequency();
            visionHeader.meanTimeConstant = -1;
            visionHeader.threshold = minSpikeFindingSigma;
            visionHeader.electrodeUsage = null;
            visionHeader.minCovarianceSpikes = -1;
            visionHeader.maxCovarianceSpikes = TRUNCATE_DATA;
            visionHeader.nlPoints = nlPoints;
            visionHeader.nrPoints = nrPoints;
            visionHeader.nDimensions = nDimensions;

            visionHeader.binsPerDimension = EMUserInterface.binsPerDimension;
            visionHeader.clusteringSignificance = EMUserInterface.clusteringSignificance;
            visionHeader.minClusters = EMUserInterface.minClusters;
            visionHeader.maxClusters = maxClusters;
            visionHeader.nEMSpikes = TRUNCATE_DATA;
            visionHeader.minNeuronSpikes = -1;
            visionHeader.acfT1 = ACFT1;
            visionHeader.acfT2 = ACFT2;
            visionHeader.maxContamination = -1;
            visionHeader.coincidenceTime = coincidenceTime;
            visionHeader.maxCorrelation = maxCorrelation;

            try {
                neuronFile = new NeuronFile(neuronFileName.getAbsolutePath(),
                        visionHeader, 5000,
                        spikeTimeList, nTTL);
            } catch (IOException e) {
                Vision.reportFatalException(
                        "The neurons (.neurons) file cannot be created", e);
            }

            System.out.println("Openning/Creating the report folder...");
            reportFolder.mkdir();
            try {
                html = new PrintWriter(new FileWriter(htmlFileName), true);
                html.println(
                        "<html><head><title>Serial Neuron Finding Report</title></head><body>");
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
            } catch (IOException e) {
                Vision.reportFatalException(
                        "The report (.html) file cannot be created or written", e);
            }

            System.out.println("Openning/Creating the EI file...");
            try {
                imagingFile = new PhysiologicalImagingFile(
                        eiFileName.getAbsolutePath(), nlPointsEI, nrPointsEI,
                        header.getArrayID());
            } catch (IOException e) {
                Vision.reportFatalException(
                        "The electrode immaging (.ei) file cannot be created", e);
            }

            System.out.println("Openning/Creating the Model file...");
            VisionHeader header = new VisionHeader(
                    visionHeader);
            header.nlPointsEI = nlPointsEI;
            header.nrPointsEI = nrPointsEI;
            //header.minimizationInterval = minimizationInterval;
            header.minimizationError = minimizationError;
            header.maxElectrodePatternSize = maxElectrodePatternSize;
            header.minThreshold = minSpikeFindingSigma;
            try {
                modelFile = new ClusteringModelFile(
                        modelFileName.getAbsolutePath(), 5000, header);
            } catch (IOException e) {
                Vision.reportFatalException(
                        "The model file (.model) file cannot be created", e);
            }
        } else {
            JOptionPane.showMessageDialog(
                    mainFrame, "Some of the SNF output files exist, but others don't. " +
                    "\n" +
                    "\nNeuron File - " + IOUtil.isValidFile(neuronFileName) +
                    "\nReport Folder - " + IOUtil.isValidFolder(reportFolder) +
                    "\nHTML Report File - " + IOUtil.isValidFile(htmlFileName) +
                    "\nEI File - " + IOUtil.isValidFile(eiFileName) +
                    "\nModel File - " + IOUtil.isValidFile(modelFileName) +
                    "\nStatus File - " + IOUtil.isValidFile(statusFile) +
                    "\n" +
                    "\nDelete the existing files or put into the folder the missing files, then restart SNF." +
                    "\nSNF will now exit.",
                    "Invalid dataset folder",
                    JOptionPane.ERROR_MESSAGE);
            return -1;
        }

        try {
            minSpikeFindingSigma = modelFile.getUserHeader().minThreshold;
        } catch (IOException e) {
            Vision.reportFatalException(
                    "The model file (.model) file cannot be read", e);
        }
        currentExtractionID = modelFile.getNextExtractionID();
        System.out.println("currentExtractionID = " + currentExtractionID);

        return 0;
    }


    private void showEI(int neuron) {
        if (averageEI[neuron] == null) {
            calculateEI(neuron);
        }

        new PhysiologicalImagePanel(
                averageEI[neuron], sigmas, minSpikeFindingSigma, electrodeMap,
                currentElectrode).showInAWindow(PlotUtil.getColorName(neuron));
    }


    public int getNeuronSpikeTimes(int neuron) {
        int n = 0;
        for (int sIndex = 0; sIndex < nSpikes; sIndex++) {
            if (spikeID[sIndex] == neuron) {
                exactNeuronSpikeTimes[n] = spikeTimeList[sIndex];
                n++;
            }
        }
        return n;
    }


    public double[] getNeuronSpikeTimes1(int neuron) {
        DoubleList list = new DoubleList();
        for (int sIndex = 0; sIndex < nSpikes; sIndex++) {
            if (spikeID[sIndex] == neuron) {
                list.add(spikeTimeList[sIndex]);
            }
        }
        return list.toArray();
    }


    public abstract class UserInterface
    implements KeyListener {

        protected JPanel leftPanel, rightPanel, controlPanel, clusteringPanel;
        protected JButton extractClusterBox;
        protected JComboBox fitPatternBox;
        protected PlotPanel[] clusteringPlots;
        public final ScatterPlotStyle[] clusteringPlotStyles;
        protected JCheckBox grid = new JCheckBox("Grid", false);
        protected JCheckBox axes = new JCheckBox("Axes", false);
        PlotPanel[] rawDataPanels = new PlotPanel[maxClusters];


        public void setExtractState(boolean state) {
            extractClusterBox.setEnabled(state);

            if (nSamplesToAnalyze != nSamples) {
                extractClusterBox.setEnabled(false);
            }
        }


        public UserInterface() {
            // set up the clustering plots
            clusteringPlotStyles = new ScatterPlotStyle[maxClusters];
            for (int j = 0; j < maxClusters; j++) {
                clusteringPlotStyles[j] = new ScatterPlotStyle(
                        SymbolType.FILLED_SQUARE, 1,
                        PlotUtil.getColor(j), false, Color.black, 1);
            }

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

            extractClusterBox = new JButton("Extract");
            extractClusterBox.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    // show the check-box dialog
                    JPanel p = new JPanel(new GridLayout(0, 1));
                    for (int i = 0; i < userInterface.getNeuronsCount(); i++) {
                        JCheckBox box = new JCheckBox(PlotUtil.getColorName(i), false);
                        p.add(box);
                    }
                    p.add(new JSeparator());
                    final JCheckBox skipForeverBox = new JCheckBox("Skip Forever?", false);
                    p.add(skipForeverBox);
                    int selection = JOptionPane.showOptionDialog(
                            mainFrame, p, "Choose Which Clusters",
                            JOptionPane.OK_CANCEL_OPTION,
                            JOptionPane.QUESTION_MESSAGE, null, null, null);
                    if (selection != JOptionPane.OK_OPTION) {
                        return;
                    }

                    // get the user's response
                    final int nNeurons = userInterface.getNeuronsCount();
                    final boolean[] selected = new boolean[nNeurons];
                    for (int i = 0; i < nNeurons; i++) {
                        selected[i] = ( (JCheckBox) p.getComponent(i)).isSelected();
                    }
                    if (MathUtil.countValues(true, selected) == 0) {
                        System.out.println("You have to select at least one cluster");
                        return;
                    }

                    // do the removal
                    SwingWorker extractWorker = new SwingWorker() {
                        public Object construct() {
                            lock();
                            canExitSNF = false;
                            final double t1 = System.currentTimeMillis();

                            try {
                                extractNeurons(selected);
                                if (skipForeverBox.isSelected()) {
                                    skipElectrodeForever[currentElectrode] = true;
                                }
                            } catch (IOException e) {
                                Vision.reportFatalException(
                                        "Error during neuron extraction", e);
                            }

                            try {
                                saveStatus();
                            } catch (IOException e) {
                                Vision.reportException(
                                        "The status file could not be saved", e);
                            }

                            final double t2 = System.currentTimeMillis();
                            System.out.println("Removal Took " + (t2 - t1) / 1000.0 +
                            " sec. ");
                            canExitSNF = true;
                            return null;
                        }


                        public void finished() {
                            nextElectrode(CalculateMode.SWITCH_ELECTRODE);
                        }
                    };
                    extractWorker.start();
                }
            });
            controlPanel.add(extractClusterBox);

            fitPatternBox = new JComboBox();
            fitPatternBox.addItemListener(new ItemListener() {
                public void itemStateChanged(ItemEvent e) {
                    if (e.getStateChange() != ItemEvent.SELECTED) {
                        return;
                    }
                    final int neuron = ( (JComboBox) e.getSource()).getSelectedIndex() -
                    1;
                    if (neuron < 0) {
                        return;
                    }

                    lock();

                    SwingWorker worker = new SwingWorker() {
                        public Object construct() {
                            System.out.println("Switching pattern to match " + neuron);
                            if (averageEI[neuron] == null) {
                                calculateEI(neuron);
                            }

                            int[] newElectrodes = new int[maxElectrodePatternSize];
                            newElectrodes[0] = currentElectrode;
                            double[] amp = new double[nElectrodes];
                            for (int electrode = 1; electrode < nElectrodes; electrode++) {
                                amp[electrode] =
                                    -MathUtil.min(averageEI[neuron][0][electrode]);
                            }

                            for (int i = 1; i < newElectrodes.length; ) {
                                int electrode = MathUtil.maxIndex(amp);
                                amp[electrode] = 0;
                                if (!electrodeMap.isDisconnected(electrode) &&
                                        (electrode != newElectrodes[0])) {
                                    newElectrodes[i++] = electrode;
                                }
                            }

                            if (MathUtil.areEqual(newElectrodes, electrodes)) {
                                System.out.println("There is no better pattern.");
                                unlock();
                            } else {
                                System.out.print("The best pattern is: ");
                                IOUtil.printArray(newElectrodes);
                                electrodes = newElectrodes;

                                nextElectrode(CalculateMode.RELOAD_DATA);
                            }

                            return null;
                        }


                        public void finished() {}
                    };
                    worker.start();
                };
            });
            controlPanel.add(fitPatternBox);

            leftPanel = new JPanel(new BorderLayout());
            clusteringPanel = new JPanel(new GridLayout(0, 1));
            leftPanel.add(clusteringPanel, BorderLayout.CENTER);

            leftPanel.add(controlPanel, BorderLayout.SOUTH);
            rightPanel = new JPanel();
        }


        private void addPlots(int plotIndex) {
            // add the plots
            for (int j = 0; j < maxClusters; j++) {
                clusteringPlots[plotIndex].addData(
                        getClusteringPlot(j, dim[plotIndex][0], dim[plotIndex][1]),
                        clusteringPlotStyles[j]);
            }
        }


        public void repaintData() {
            for (int i = 0; i < nShownPlots; i++) {
                clusteringPlots[i].replotAllData();
            }
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


        // 0 - partial
        // 1 - all
        // 2 - no autoscale
        boolean DEBUG = false;
        public final void updateView(boolean complete, boolean autoscale) {
            if (DEBUG) {
                System.out.println("updateView - in");
            }
            int nGaussians = getNeuronsCount();

            fitPatternBox.removeAllItems();
            fitPatternBox.addItem("Match");
            viewMenu.removeAll();
            for (int neuron = 0; neuron < nGaussians; neuron++) {
                fitPatternBox.addItem("" + PlotUtil.getColorName(neuron));
                JMenuItem item = new JMenuItem("Show EI for " +
                        PlotUtil.getColorName(neuron));
                item.setName("" + neuron);
                item.addActionListener(viewListener);
                viewMenu.add(item);
            }

            if (complete) {
                // clear the views
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
                    ArrayList additionalLegend = new ArrayList();
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
                        if (neuronSpikeCount[neuron] > 2) {
                            calculateEI(neuron);
                            JComponent eiPanel = new PhysiologicalImagePanel(
                                    averageEI[neuron], sigmas, minSpikeFindingSigma,
                                    electrodeMap,
                                    currentElectrode);
                            rightPanel.add(eiPanel);
                        } else {
                            rightPanel.add(new JLabel("    Not enough spikes to compute."));
                        }
                    }

                    // add the Raw Data
                    if (showRawData.getValue()) {
                        rawDataPanels[neuron] = getWaveformOverlap(neuron);
                        rightPanel.add(rawDataPanels[neuron]);
                    }

                    // add the amplitude histogram
                    amplitudePanel[neuron] = new PlotPanel();
                    amplitudePanel[neuron].addBackgroundText(KeyUtil.getKeyText(neuron));
                    amplitudePanel[neuron].addData(
                            neuronAmplitudeHistograms[neuron], new HistogramStyle());
                    if (showGaussianFit.isSelected()) {
                        amplitudePanel[neuron].addData(
                                amplitudeHistogramFit[neuron],
                                new FunctionStyle("Amplitude Hist"));
                        ArrayList legend = new ArrayList();
                        legend.add("Mean: " + StringUtil.format(
                                amplitudeHistogramFit[neuron].getMean(), 1));
                        legend.add("RMS: " + StringUtil.format(
                                amplitudeHistogramFit[neuron].getSigma(), 1));
                        legend.add("Chi2: " + StringUtil.format(
                                amplitudeHistogramFit[neuron].getChiSquared(), 2));
                        amplitudePanel[neuron].setAdditionalLegend(legend);
                    }
                    amplitudePanel[neuron].autoscale();
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

                ArrayList additionalLegend = new ArrayList();
                additionalLegend.add("View: " + (dim[i][0] + 1) + " - " + (dim[i][1] + 1));
                clusteringPlots[i].setAdditionalLegend(additionalLegend);

                if (autoscale) {
                    clusteringPlots[i].autoscale();
                }

                clusteringPlots[i].replotAllData();
            }

            if (DEBUG) {
                System.out.println("updateView - out");
            }
        }


        public abstract PlotData[] getClusterContours(int plotIndex);


        public abstract PlotStyle getContoursStyle();


        public abstract ClusteringModelFile.Model getModel();


        public abstract int getNeuronsCount();


        public abstract void sortSpikes();


        public final void updateNeuronsInfo() {
            if (DEBUG) {
                System.out.println("updateNeuronsInfo - in");
            }

            final int nNeurons = getNeuronsCount();
            sortSpikes();

            // extract the neuron times and compute the autocorrelations
            for (int neuron = 0; neuron < nNeurons; neuron++) {
                averageEI[neuron] = null;

                // calculate the spike amplitude histograms
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

                // calculate the autocorrelation
                autocorrelations[neuron] = AutocorrelationCalculator.calculateISI(
                        exactNeuronSpikeTimes, N, 50, 0.05);
                nBadSpikes[neuron] = autocorrelations[neuron].count(ACFT1, ACFT2);
                double T = nSamples / VisionParams.SAMPLES_PER_MILLISECOND;
                double dT = ACFT2 - ACFT1;
                contaminationIndex[neuron] = nBadSpikes[neuron] * T / (dT * N * N);

                // calculate the gaussian fit
                if (showGaussianFit.isSelected()) {
                    try {
                        amplitudeHistogramFit[neuron] =
                            neuronAmplitudeHistograms[neuron].fitToGaussian(10);
                    } catch (FitFailedException e) {
                        amplitudeHistogramFit[neuron] = new Gaussian1DFunction(0, 1, 1, 1);
                    }
                }

                // calculate the cross-correlation with previous neurons
                // FIXME
//				getCorrelation(exactNeuronSpikeTimes)
            }

            if (DEBUG) {
                System.out.println("updateNeuronsInfo - out");
            }
        }


        public abstract void dataUpdated();


        public abstract ScatterPlotData getClusteringPlot(
                final int neuronID, final int d1, final int d2);


        public void lock() {
            myGlassPane.setVisible(true);

            for (int i = 0; i < menuList.length; i++) {
                menuList[i].setEnabled(false);
            }
        }


        public void unlock() {
            for (int i = 0; i < menuList.length; i++) {
                menuList[i].setEnabled(true);
            }

            myGlassPane.setVisible(false);
        }


        public void keyPressed(KeyEvent e) {
//			switch (e.getKeyCode()) {
//case KeyEvent.VK_F5:
//	skipElectrode();
//break;
//}

// CTRL + Key is used to show the EI
            if ( (e.getModifiersEx() & e.ALT_DOWN_MASK) == e.ALT_DOWN_MASK) {
                int nNeurons = userInterface.getNeuronsCount();
                int code = e.getKeyCode();
                for (int i = 0; i < nNeurons; i++) {
                    if (code == KeyUtil.getKeyCode(i)) {
                        showEI(i);
                    }
                }
            }

            // SHIFT + Key is used to highlight a cluster
            if ( (e.getModifiersEx() & e.SHIFT_DOWN_MASK) == e.SHIFT_DOWN_MASK) {
                int nNeurons = userInterface.getNeuronsCount();
                int code = e.getKeyCode();
                for (int i = 0; i < nNeurons; i++) {
                    if (code == KeyUtil.getKeyCode(i)) {
                        clusteringPlotStyles[i].setSymbolSize(2);
                        repaintData();
                    }
                }
            } else {
                for (int i = 0; i < clusteringPlotStyles.length; i++) {
                    clusteringPlotStyles[i].setSymbolSize(1);
                }
                repaintData();
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
        public final static int binsPerDimension = 20;
        public final static double clusteringSignificance = 3;
        public final static int maxSpikeLoss = 0;
        public final static int minClusters = 2;
        public static final int maxEMIterations = 100;

        public JButton refitButton;
        private double[] min, max, binInterval;
        private ExpectationMaximization expectationMaximization;
        ArrayList guessMeans = new ArrayList();
        ArrayList guessSigmas = new ArrayList();


        public EMUserInterface(int type, float[][] _pJX) {
            super();

            switch (type) {
            case 1:
                expectationMaximization = new ExpectationMaximization1(
                        nDimensions, maxClusters, _maxSpikes, _pJX);
                break;

            case 2:
                expectationMaximization = new ExpectationMaximization2(
                        nDimensions, maxClusters, _maxSpikes, _pJX);
                break;
            }


            min = new double[nDimensions];
            max = new double[nDimensions];
            binInterval = new double[nDimensions];

            // set up the clustering plots
            for (int i = 0; i < nShownPlots; i++) {
                // clicks are used to add or remove gaussians
                clusteringPlots[i].addClickListener(new ScaledMouseListener() {
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
//								mean[d] = (min[d] + max[d]) * 0.5;
                                mean[d] = 0;
                                sigma[d] = binInterval[d] * binInterval[d];
                            }

                            mean[dim[i][0]] = x;
                            mean[dim[i][1]] = y;

                            expectationMaximization.addGaussian(mean, sigma);
                        }

                        updateView(false, false);

//						refitButton.setEnabled(
//						expectationMaximization.getClustersCount() != 0);
                        setExtractState(false);
                        fitPatternBox.setEnabled(false);
                    }

                    public void enteredPerformed(Component source,
                            MouseEvent event, double x, double y) {	
                    }

                    public void exitedPerformed(Component source,
                            MouseEvent event, double x, double y) {	
                    }

                    public void pressPerformed(Component source,
                            MouseEvent event, double x, double y) {
                    }

                    public void releasePerformed(Component source,
                            MouseEvent event, double x, double y) {
                    }
                    

                    
                    
                });
            } // for i

            refitButton = new JButton("Fit");
            refitButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    SwingWorker worker = new SwingWorker() {
                        public Object construct() {
                            if (DEBUG) {
                                System.out.println("refitButton - started");
                            }

                            lock();
                            doEMFit();

                            updateNeuronsInfo();
                            return null;
                        }


                        public void finished() {
                            if (DEBUG) {
                                System.out.println("refitButton - finished");
                            }

                            updateView(true, true);
                            unlock();
                        }
                    };
                    worker.start();
                }
            });
            controlPanel.add(refitButton);
        }


        public PlotPanel makePlotPanel(int i) {
            return new PlotPanel();
        }


        public ClusteringModelFile.Model getModel() {
            ClusteringModelFile.EMModel m = new ClusteringModelFile.EMModel();

            double[][] eigenvectors = new double[nDimensions][];
            for (int i = 0; i < nDimensions; i++) {
                eigenvectors[i] = pca.getEigenVector(i);
            }

            m.cleaningLevel = cleaningLevel;
            m.electrodes = electrodes;
            m.threshold = sigmas[currentElectrode] * floatingSigma;
            m.nDimensions = expectationMaximization.getDimensionsCount();
            m.eigenvectors = eigenvectors;
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


        private void doEMFit() {
            if (expectationMaximization.getClustersCount() == 0) {
                for (int i = 0; i < guessMeans.size(); i++) {
                    double[] m = (double[]) guessMeans.get(i);
                    double[] c = (double[]) guessSigmas.get(i);
                    expectationMaximization.addGaussian(m, c);
                }
            }

            try {
                expectationMaximization.fit(1e-6, 5, maxEMIterations);
            } catch (FitFailedException e) {
                System.out.println("EM fit failed because: " + e.getMessage() +
                        " Try fitting again.");
            }
        }


        public void dataUpdated() {
            expectationMaximization.reset(projectedSpikes, TRUNCATE_DATA, nSpikes);
            guessMeans.clear();
            guessSigmas.clear();

            PCANFClustering.densityClustering(projectedSpikes, nSpikes, binsPerDimension,
                    maxSpikeLoss, clusteringSignificance,
                    guessMeans,
                    guessSigmas, min, max, binInterval, false,
                    1, maxClusters);

            setExtractState(false);
            fitPatternBox.setEnabled(false);
            refitButton.setEnabled(true);
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
            super.lock();

            refitButton.setEnabled(false);
            setExtractState(false);
            fitPatternBox.setEnabled(false);
        }


        public void unlock() {
            super.unlock();

            refitButton.setEnabled(true);
            setExtractState(true);
            fitPatternBox.setEnabled(true);
        }


        public void keyPressed(KeyEvent keyEvent) {
            super.keyPressed(keyEvent);

            switch (keyEvent.getKeyCode()) {
            case KeyEvent.VK_F1:
                refitButton.doClick();
                break;

            case KeyEvent.VK_F2:
                extractClusterBox.doClick();
                break;
            }
        }
    }


    public class ManualUserInterface
    extends UserInterface implements ChangeListener {


        public ManualUserInterface() {
            super();

            for (int i = 0; i < clusteringPlots.length; i++) {
                clusteringPlots[i].addPassiveSelectionChangeListener(this);
            }
        }


        public void stateChanged(ChangeEvent e) {
            unlock();
            updateNeuronsInfo();
            updateView(true, false);
        }


        public PlotPanel makePlotPanel(int i) {
            return new PlotPanel("", true, true, true, true);
        }


        public ClusteringModelFile.Model getModel() {
            ClusteringModelFile.ManualModel m = new ClusteringModelFile.ManualModel();

            double[][] eigenvectors = new double[nDimensions][];
            for (int i = 0; i < nDimensions; i++) {
                eigenvectors[i] = pca.getEigenVector(i);
            }

            m.cleaningLevel = cleaningLevel;
            m.electrodes = electrodes;
            m.threshold = sigmas[currentElectrode] * floatingSigma;
            m.nDimensions = nDimensions;
            m.eigenvectors = eigenvectors;

            m.nClusters = getNeuronsCount();
            m.dimension1 = new int[m.nClusters];
            m.dimension2 = new int[m.nClusters];
            m.selections = new ArrayList();
            for (int n = 0, j = 0; n < clusteringPlots.length; n++) {
                for (int i = 0; i < clusteringPlots[n].getPassiveSelectionsCount(); i++) {
                    m.dimension1[j] = dim[n][0];
                    m.dimension2[j] = dim[n][1];
                    m.selections.add(clusteringPlots[n].getPassiveSelection(i).toArray());
                    j++;
                }
            }

            return m;
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


        public int getNeuronsCount() {
            int N = 0;
            for (int n = 0; n < clusteringPlots.length; n++) {
                N += clusteringPlots[n].getPassiveSelectionsCount();
            }
            return N + 1;
        }


        public void dataUpdated() {
            for (int n = 0; n < clusteringPlots.length; n++) {
                clusteringPlots[n].resetPassiveSelections();
            }

            updateNeuronsInfo();

            setExtractState(false);
            fitPatternBox.setEnabled(false);
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
            super.lock();

            setExtractState(false);
            fitPatternBox.setEnabled(false);
        }


        public void unlock() {
            super.unlock();

            setExtractState(true);
            fitPatternBox.setEnabled(true);
        }


        public void keyPressed(KeyEvent keyEvent) {
            super.keyPressed(keyEvent);

            switch (keyEvent.getKeyCode()) {
            case KeyEvent.VK_F2:
                extractClusterBox.doClick();
                break;
            }
        }

    }


    private void loadStatus(File f) throws IOException {
        RandomAccessFile s = new RandomAccessFile(f, "rw");

        int version = s.readInt();
        cleaningLevel = s.readInt();

        for (int i = 0; i < nElectrodes; i++) {
            skipElectrode[i] = s.readBoolean();
        }
        for (int i = 0; i < nElectrodes; i++) {
            skipElectrodeForever[i] = s.readBoolean();
        }

        int nNeurons = s.readInt();
        for (int i = 0; i < nNeurons; i++) {
            neuronsToRemoveList.add(s.readInt());
        }

        floatingSigma = s.readDouble();
        electrodeUsage = ElectrodeUsage.values()[s.readInt()];
        _maxAmplitude = s.readInt();
        _showEI = s.readBoolean();
        _autoOutline = s.readBoolean();
        _autoLabels = s.readBoolean();

        // do not change ANYTHING above

        s.close();

        // do the checks
        if (version != LATEST_VERSION_TAG || cleaningLevel < 0 || nNeurons < 0 ||
                floatingSigma < 0 || electrodeUsage == null || _maxAmplitude < 0) {
            throw new IOException("Bad status file.");
        }
    }


    private void saveStatus() throws IOException {
        String statusFileName = datasetFolder + File.separator + datasetName + ".snfs";
        RandomAccessFile s = new RandomAccessFile(statusFileName, "rw");
        s.setLength(0);

        s.writeInt(LATEST_VERSION_TAG);
        s.writeInt(cleaningLevel);

        for (int i = 0; i < nElectrodes; i++) {
            s.writeBoolean(skipElectrode[i]);
        }
        for (int i = 0; i < nElectrodes; i++) {
            s.writeBoolean(skipElectrodeForever[i]);
        }

        s.writeInt(neuronsToRemoveList.size());
        for (int i = 0; i < neuronsToRemoveList.size(); i++) {
            s.writeInt(neuronsToRemoveList.get(i));
        }

        s.writeDouble(floatingSigma);
        s.writeInt(electrodeUsage.ordinal());
        s.writeInt(maxAmplitude.getValue());
        s.writeBoolean(showEI.getValue());
        s.writeBoolean(autoOutline.getValue());
        s.writeBoolean(autoLabels.getValue());

        // do not change ANYTHING above

        s.close();
    }


    private void copySpikes(int electrodeIndex) {
        for (int sIndex = 0; sIndex < nSpikes; sIndex++) {
            // roundoff error can happen
            int time = (int) Math.round(spikeTimeList[sIndex]);

            for (int i = 0; i < nPoints + 2; i++) {
                spikes1[electrodeIndex][i][sIndex] = rawData[time - nlPoints - 1 + i];
            }
        }
    }


    private void loadRawData(boolean doSpikeFinding) throws IOException {
        double t = 0;

        for (int electrodeIndex = 0; electrodeIndex < electrodes.length; electrodeIndex++) {
            double t1 = System.currentTimeMillis();

            // MOD
//			rawDataFile.readRawData(electrodes[electrodeIndex], rawData);
            rawDataFile.readRawData(electrodes[electrodeIndex], rawData,
                    nSamplesToAnalyze);

            t += (System.currentTimeMillis() - t1);

            if (electrodeIndex == 0 && doSpikeFinding) {
                nSpikes = SpikeFinder.findSpikes(
                        rawData, sigmas[currentElectrode] * floatingSigma,
                        spikeTimeList, spikeAmplitudeList, ignoreAtBeginning, 2 * nrPointsEI);
                System.out.println(
                        "Spikes: " + nSpikes + ", Threshold: " +
                        StringUtil.format(floatingSigma, 2) + " : " +
                        StringUtil.format(sigmas[currentElectrode] * floatingSigma, 1));
            }
            copySpikes(electrodeIndex);
        }

        t /= 1000;
        double mBytesRead = nSamples * 1.5 * electrodes.length / 1024.0 / 1024.0;
        System.out.println(
                "Data reading speed: " + StringUtil.format(mBytesRead / t, 1) +
        " MBytes/sec.");
    }


    class CovarianceThread
    extends Thread {

        private int n1, n2;
        private CovarianceMatrix covariance;
        private short[] s;
        private UniformSpline spline;
        private boolean isWorking = true;
        private double[] allignedSpike;
        private BinaryRandom r;


        public CovarianceThread(int n1, int n2, int nSpikesToUse) {
            this.n1 = n1;
            this.n2 = n2;
            covariance = new CovarianceMatrix(electrodes.length * nPoints);
            spline = new UniformSpline(nPoints + 2);
            allignedSpike = new double[nPoints * maxElectrodePatternSize];
            r = new BinaryRandom( (double) Math.min(nSpikesToUse, n2 - n1) /
                    (double) nSpikes);
        }


        public void run() {
            double[] spike = new double[nPoints + 2];

            for (int sIndex = n1; sIndex < n2; sIndex++) {
                if (r.random()) {
                    int time = (int) Math.round(spikeTimeList[sIndex]); // roundoff error can happen
                    for (int electrode = 0; electrode < electrodes.length; electrode++) {
                        if (time > nSamples - nrPoints) {
                            continue;
                        }

//						spline.reSpline(spikes[sIndex][electrode]);
                        for (int i = 0; i < nPoints + 2; i++) {
                            spike[i] = spikes1[electrode][i][sIndex];
                        }
                        spline.reSpline(spike);

                        final int baseIndex = electrode * nPoints;
                        for (int i = 0; i < nPoints; i++) {
                            allignedSpike[baseIndex + i] =
                                spline.getValueAt(t0[sIndex] + (i - nlPoints) - 1);
                        }
                    }
                    covariance.addData(allignedSpike);
                }
            }
            isWorking = false;
        }
    }


    private float[] calculateCovarianceMatrix() {
        System.out.print("Calculating covariances: ");
        double t1 = System.currentTimeMillis();

        CovarianceThread covThread1 = new CovarianceThread(
                0, nSpikes / 2, TRUNCATE_DATA / 2);
        CovarianceThread covThread2 = new CovarianceThread(
                nSpikes / 2, nSpikes, TRUNCATE_DATA / 2);
        covThread1.start();
        covThread2.start();
        while (covThread1.isWorking || covThread2.isWorking) {
            try {
                Thread.sleep(100);
            } catch (InterruptedException e) {
            }
        }

        double t2 = System.currentTimeMillis();
        System.out.println("took " + (t2 - t1) / 1000);

        covThread1.covariance.combine(covThread2.covariance);
        return covThread1.covariance.getCovariance();
    }


    /*
        public static void calculateExactSpikeTimes(
            int nlPoints, int nrPoints, int nSpikes, int nSamples,
            double[] spikeTimeList, float[] spikeAmplitudeList, double[] t0,
            short[][][] spikes, int minimizationInterval, double minimizationError
            ) {
            System.out.print("Calculating exact spike times: ");
            double t1 = System.currentTimeMillis();
            int nPoints = nlPoints + nrPoints + 1;
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
                spline.reSpline(spikes[sIndex][0]);
                t0[sIndex] = FunctionMinimization.brentParabolic(
                    spline, nlPoints + 1 - minimizationInterval,
                    nlPoints + 1 + minimizationInterval, minimizationError);
                spikeTimeList[sIndex] = (double) (time - nlPoints - 1 + t0[sIndex]);
                if (time > 5320 && time < 5330) {
                    System.out.println("==> " + time + " : " + spikeTimeList[sIndex]);
                }
                spikeAmplitudeList[sIndex] = - (float) spline.getValueAt(t0[sIndex]);
                if (Math.abs(spikeTimeList[sIndex] - time) >= 1) {
                    //System.out.println("problem -- " + (spikeTimeList[sIndex] - time));
                    spikeTimeList[sIndex] = time;
                }
            }
            double t2 = System.currentTimeMillis();
            System.out.println( + (t2 - t1) / 1000 + " sec.");
        }
     */

    class ProjectingThread
    extends Thread {

        private int n1, n2;
        private float[] uniformX;
        private UniformSpline spline;
        private boolean isWorking = true;
        double[] allignedSpike;


        public ProjectingThread(int n1, int n2) {
            this.n1 = n1;
            this.n2 = n2;
            uniformX = new float[nPoints];
            for (int i = 0; i < uniformX.length; i++) {
                uniformX[i] = i;
            }
            spline = new UniformSpline(nPoints + 2);
            allignedSpike = new double[nPoints * maxElectrodePatternSize];
        }


        public void run() {
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
                    projectedSpikes[d][sIndex] =
                        (float) pca.project(allignedSpike, d);
                }
            }

            isWorking = false;
        }
    }


    private void projectSpikes() {
        System.out.print("Projecting spikes: ");
        double t1 = System.currentTimeMillis();

        ProjectingThread prjThread1 = new ProjectingThread(0, nSpikes / 2);
        ProjectingThread prjThread2 = new ProjectingThread(nSpikes / 2, nSpikes);
        prjThread1.start();
        prjThread2.start();

        while (prjThread1.isWorking || prjThread2.isWorking) {
            try {
                Thread.sleep(100);
            } catch (InterruptedException e) {
            }
        }

        double t2 = System.currentTimeMillis();
        System.out.println( (t2 - t1) / 1000 + " sec.");
    }


    private void nextElectrode(final CalculateMode mode) {
        new SwingWorker() {
            public Object construct() {
                canExitSNF = false;
                userInterface.lock();

                switch (mode) {
                case SWITCH_ELECTRODE:

                    // MOD
                    nSamplesToAnalyze = nSamplesParameter.getValue();

                    // reset sigma
                    floatingSigma = minSpikeFindingSigma;

                    currentElectrode = elecrodeDialog.geNextElectrode();
                    if (currentElectrode == -1) {
                        return null;
                    }
                    electrodes = electrodeMap.getAdjacentsTo(currentElectrode,
                            electrodeUsage.ordinal());
                    System.out.print("\nElectrode Pattern: ");
                    IOUtil.printArray(electrodes);

                    try {
                        loadRawData(true);
                    } catch (IOException e) {
                        Vision.reportFatalException(
                                "The .rem file could not be read",
                                e);
                    }
                    break
                    ;

                case RELOAD_DATA:
                case NEW_SPIKE_FINDING:
                    System.out.print("\nElectrode Pattern: ");
                    IOUtil.printArray(electrodes);
                    try {
                        loadRawData(true);
                    } catch (IOException e) {
                        Vision.reportFatalException(
                                "The .rem file could not be read",
                                e);
                    }
                    break
                    ;

//					case ISOLATE_CLUSTER:
//	break;

                case RECLUSTER:
                    break;
                }

                try {
                    SerialMappingPanel.calculateExactSpikeTimes(
                            nlPoints, nrPoints, nSpikes, nSamples,
                            spikeTimeList, spikeAmplitudeList, t0,
                            spikes1, minimizationInterval, minimizationError);
                } catch (CannotEvaluateException e) {
                    e.printStackTrace();
                }

                float[] covariance = calculateCovarianceMatrix();
                pca = new PCA(covariance);
                try {
                    pca.doPCA();
                } catch (TooManyIterationsException e) {
                    Vision.reportException("PCA analysis failed", e);
                }
                pca.printPercentageEigenvalues(nDimensions);
                projectSpikes();

                userInterface.dataUpdated();

                return "ok";
            }


            public void finished() {
                if (this.getValue() != null) {
                    userInterface.updateView(true, true);
                    setMenuBarEnabled(true);

                    mainFrame.setTitle(
                            "Manual Classifier - " + datasetFolder +
                            ", El: " + currentElectrode + ", Spikes: " + nSpikes);

                    canExitSNF = true;
                    userInterface.unlock();
                }
            }
        }

        .start();
    }


    private void setMenuBarEnabled(boolean enabled) {
        for (int i = 0; i < menuList.length; i++) {
            menuList[i].setEnabled(enabled);
        }
    }


    public static void showData(short[] rawData, int n1, int n2, String text) {
        ScatterPlot h = new ScatterPlot("");
        for (int i = n1; i < n2; i++) {
            h.add(i, rawData[i]);
        }
        ScatterPlotStyle style = new ScatterPlotStyle();
        style.setConnectingPoints(true);
        style.setSymbolType(SymbolType.NONE);
        PlotUtil.showData(text, h, style);
    }


    private void cleaning() {
        canExitSNF = false;
        userInterface.lock();

        // write a HTML line for cleaning
        html.println("<tr>");
        html.println("<td>Cleaning");

        System.out.println(
                "\nRemoving the " + neuronsToRemoveList.size() + " waiting neurons ...");

        // run the cleaning in a new thread
        SwingWorker sw = new SwingWorker() {
            public Object construct() {
                try {
                    saveStatus();
                } catch (IOException e) {
                    Vision.reportException("The status file could not be saved",
                            e);
                }

                try {
                    cleaning(
                            neuronsToRemoveList,
                            neuronFile,
                            imagingFile,
                            nlPointsEI,
                            nrPointsEI,
                            minimizationInterval,
                            minimizationError,
                            electrodeMap,
                            rawData,
                            rawDataFile,
                            sigmas,
                            minSpikeFindingSigma,

                            spikeTimeList,
                            spikeAmplitudeList,
                            amplitudeHistogram,
                            skipElectrode,
                            isIgnored,
                            nSamples
                    );
                } catch (IOException e) {
                    Vision.reportFatalException(
                            "The .rem file could not be accessed during cleanning", e);
                }

                try {
                    saveAmplitudeHistograms(amplitudeHistogram,
                            amplitudeHistogramsFileName);
                } catch (IOException e) {
                    Vision.reportException(
                            "The amplitude histograms could not be saved", e);
                }

                try {
                    saveStatus();
                } catch (IOException e) {
                    Vision.reportException("The status file could not be saved",
                            e);
                }

                return null;
            }


            public void finished() {
                cleaningLevel++;
                canExitSNF = true;
                userInterface.unlock();
                nextElectrode(CalculateMode.SWITCH_ELECTRODE);
            }
        };
        sw.start();
    }


    public static void cleaning(
            IntegerList neuronsToRemoveList,
            NeuronFile neuronFile,
            PhysiologicalImagingFile imagingFile,
            final int nlPointsEI,
            final int nrPointsEI,
            final int minimizationInterval,
            final double minimizationError,
            ElectrodeMap electrodeMap,
            final short[] rawData,
            ElectrodeMajorRawDataFile rawDataFile,
            float[] sigmas,
            double minSpikeFindingSigma,

            double[] spikeTimeList,
            float[] spikeAmplitudeList,
            int[][] amplitudeHistogram,
            boolean[] skipElectrode,
            boolean[] isIgnored,
            int nSamples

    ) throws IOException {

        double t1 = System.currentTimeMillis();
        final int nNeurons = neuronsToRemoveList.size();
        // load the spike times
        double[][] times = new double[nNeurons][];
        int[][] intTimes = new int[nNeurons][];
        for (int i = 0; i < nNeurons; i++) {
            times[i] = neuronFile.getExactSpikeTimes(neuronsToRemoveList.get(i));
            intTimes[i] = new int[times[i].length];
            for (int j = 0; j < intTimes[i].length; j++) {
                intTimes[i][j] = (int) Math.round(times[i][j]);
            }
        }

        // prepare structures
        int nPointsEI = nlPointsEI + nrPointsEI + 1;
        final UniformSpline spline = new UniformSpline(nPointsEI);
        int nElectrodes = electrodeMap.getNumberOfElectrodes();
        float[][][][] eiImages = new float[nNeurons][2][nElectrodes][];

        // calculate the EI's exact spike time on the main electrode
        double[] t0EI = new double[nNeurons];
        for (int neuron = 0; neuron < nNeurons; neuron++) {
            int mainElectrode = neuronFile.getElectrode(neuronsToRemoveList.get(
                    neuron));
            System.out.println("main " + mainElectrode);
            rawDataFile.readRawData(mainElectrode, rawData, nSamples);
            ///fuuk
            float[][] img = WaveformCalculator.calculateEI(times[neuron],
                    times[neuron].length, nEISpikes, nlPointsEI, nrPointsEI, rawData);
            spline.reSpline(img[0]);
            try {
                t0EI[neuron] = FunctionMinimization.brentParabolic(spline,
                        nlPointsEI - minimizationInterval, nlPointsEI + minimizationInterval,
                        minimizationError);
            } catch (CannotEvaluateException e) {
                e.printStackTrace();
            }
            System.out.println("t0 " + t0EI[neuron]);
        }

        // do the removal, electrode by electrode
//		int nSamples = rawDataFile.getHeader().getNumberOfSamples();

        // do the cleaning
        for (int electrode = 0; electrode < nElectrodes; electrode++) {
            // for disconnected electrodes or ignored electrodes
            if (electrodeMap.isDisconnected(electrode) || isIgnored[electrode]) {
                for (int neuron = 0; neuron < nNeurons; neuron++) {
                    eiImages[neuron][0][electrode] = new float[nPointsEI];
                    eiImages[neuron][1][electrode] = new float[nPointsEI];
                }
                continue;
            }

            // load the raw data
            rawDataFile.readRawData(electrode, rawData, nSamples);

            final double threshold = sigmas[electrode] * minSpikeFindingSigma;
            System.out.println("Electrode " + electrode);
            boolean cleaningPerformed = false;

            for (int neuron = 0; neuron < nNeurons; neuron++) {
                // calculate the EI
                float[][] img = WaveformCalculator.calculateEI(times[neuron],
                        times[neuron].length, nEISpikes, nlPointsEI, nrPointsEI, rawData);
                eiImages[neuron][0][electrode] = img[0];
                eiImages[neuron][1][electrode] = img[1];
                spline.reSpline(img[0]);

                // evaluate whether the cleaning is required for this neuron
                final int tEI = MathUtil.minIndex(img[0]);
                int nBad = 0;
                for (int i = 0; i < times[neuron].length; i++) {
                    int time = (int) Math.round(times[neuron][i] - (t0EI[neuron] - tEI));
                    if ( -rawData[time - 1] > threshold ||
                            -rawData[time + 0] > threshold ||
                            -rawData[time + 1] > threshold) {
                        nBad++;
                    }
                }

                if (100.0 * nBad / times[neuron].length > allowedCleaningError) {
                    // cleaning is required, do it
                    nBad = 0;
                    for (int sIndex = 0; sIndex < times[neuron].length; sIndex++) {
                        final double tDouble = times[neuron][sIndex];
                        final int tInt = intTimes[neuron][sIndex];
                        final int startSample = tInt - nlPointsEI + 1;
                        final int endSample = tInt + nrPointsEI - 1;
                        final int t2 = (int) tDouble + 1;
                        double maxCorrelation = Double.NEGATIVE_INFINITY, correlation = 0;
                        double ttt, n, bestT = tDouble;
                        int t, i;

                        for (n = (int) tDouble; n <= t2; n += 0.1) {
                            correlation = 0;
                            ttt = t0EI[neuron] - n;
                            for (t = startSample + 10, i = 0; t < startSample + 40; t++,
                            i++) {
                                correlation += rawData[t] * spline.getValueAt(ttt + t);
                            }
                            if (correlation > maxCorrelation) {
                                maxCorrelation = correlation;
                                bestT = n;
                            }
                        }

                        ttt = t0EI[neuron] - bestT;
                        for (t = startSample, i = 0; t < endSample; t++, i++) {
                            rawData[t] -= spline.getValueAt(ttt + t);
                        }

                        // test the goodness of the removal
                        if ( -rawData[tInt - 1] > threshold ||
                                -rawData[tInt] > threshold ||
                                -rawData[tInt + 1] > threshold) {
                            nBad++;
                        }
                    } // for sIndex

                    System.out.println(
                            neuronsToRemoveList.get(neuron) + ": " +
                            StringUtil.format(100.0 * nBad / times[neuron].length, 1) +
                            "% error");

                    cleaningPerformed = true;
                }

            } // for neuron

            // recalculate amplitude histogram and save the data if required
            if (cleaningPerformed) {
                if (amplitudeHistogram != null && spikeTimeList != null &&
                        spikeAmplitudeList != null) {
                    // update the raw data amplitude histogram
                    Arrays.fill(amplitudeHistogram[electrode], 0);
                    int n = SpikeFinder.findSpikes(
                            rawData, sigmas[electrode] * minSpikeFindingSigma,
                            spikeTimeList, spikeAmplitudeList, ignoreAtBeginning,
                            2 * nrPointsEI);
                    for (int i = 0; i < n; i++) {
                        int amp = (int) spikeAmplitudeList[i];
                        if (amp > 0 && amp < amplitudeHistogram[electrode].length) {
                            amplitudeHistogram[electrode][amp]++;
                        }
                    }
                }

                rawDataFile.writeRawData(electrode, rawData, nSamples);
            }
        } // for electrode

        // save the EIs
        for (int neuron = 0; neuron < nNeurons; neuron++) {
            int id = neuronsToRemoveList.get(neuron);
            imagingFile.appendImage(id, 0, eiImages[neuron]);
        }

        // reset the structures
        neuronsToRemoveList.clear();
        if (skipElectrode != null) {
            Arrays.fill(skipElectrode, false);
        }

        double t2 = System.currentTimeMillis();
        System.out.println("Cleaning took: " + (t2 - t1) / 1000 + " sec.");
    }


    private void extractNeurons(boolean[] extract) throws IOException {
        canExitSNF = false;

        IntegerList neuronIndex = new IntegerList();
        IntegerList neuronID = new IntegerList();

        loop:for (int neuron = 0; neuron < extract.length; neuron++) {
            if (!extract[neuron]) {
                continue;
            }

            final int N = getNeuronSpikeTimes(neuron);
            double[] tList1 = new double[N];
            System.arraycopy(exactNeuronSpikeTimes, 0, tList1, 0, N);

            // check to see if this neuron is different from all the ones waiting for removal
            double maxP = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < neuronsToRemoveList.size(); i++) {
                int id = neuronsToRemoveList.get(i);
                double[] tList2 = neuronFile.getExactSpikeTimes(id);

                DoubleHistogram ccH = CrossCorrelationCalculator.
                getCrossCorrelationHistogram(
                        tList1, tList2, 1, 400);
                double[] cc = ccH.toArray();
                double maxV = Double.NEGATIVE_INFINITY;
                for (int k = 0; k < cc.length - 2 * coincidenceTime; k++) {
                    double v = 0;
                    for (int m = 0; m < 2 * coincidenceTime; m++) {
                        v += cc[k + m];
                    }
                    if (v > maxV) {
                        maxV = v;
                    }
                }
                double p = maxV / Math.min(tList1.length, tList2.length);
                if (p > maxP) {
                    maxP = p;
                }
                if (p > maxCorrelation) {
                    PlotUtil.showData("" + p, ccH, new HistogramStyle());
                    String message = "Sorry, this neuron seems to match neuron " + id +
                    " with a correlation of " + p + "." +
                    "\nThis neuron has been marked for removal allready. " +
                    "\nChoose another one.";

                    System.out.println(message);
                    JOptionPane.showMessageDialog(
                            mainFrame, message, "Duplicate Detected",
                            JOptionPane.ERROR_MESSAGE);

                    continue loop;
                }
            }
            System.out.println("Max correlation : " + StringUtil.format(maxP, 4));

            // calculate the amplitudes on the electrode array
            if (averageEI[neuron] == null) {
                calculateEI(neuron);
            }
            int nCleanedElectrodes = 0;
            double[] amplitudeList = new double[nElectrodes];
            double[] sigmaList = new double[nElectrodes];
            for (int electrode = 1; electrode < nElectrodes; electrode++) {
                if (electrodeMap.isDisconnected(electrode)) {
                    continue;
                }
                final int minIndex = MathUtil.minIndex(averageEI[neuron][0][
                                                                            electrode]);
                amplitudeList[electrode] = -averageEI[neuron][0][electrode][minIndex];
                sigmaList[electrode] = averageEI[neuron][1][electrode][minIndex];
            }

            final int bestElectrode = MathUtil.maxIndex(amplitudeList);
            if (amplitudeList[currentElectrode] < amplitudeList[bestElectrode]) {
                String message =
                    "This neuron may be best extracted from another electrode with higher amplitude: " +
                    "\nThis Electrode ( " + currentElectrode + ") : " +
                    amplitudeList[currentElectrode] +
                    "\nBest Electrode ( " + bestElectrode + ") : " +
                    amplitudeList[bestElectrode] +
                    "\nDo you really want to extract the neuron?";

                int r = JOptionPane.showOptionDialog(
                        mainFrame, message, "Possible Inefficiency",
                        JOptionPane.YES_NO_OPTION,
                        JOptionPane.QUESTION_MESSAGE, null, null, null);

                if (r == JOptionPane.NO_OPTION) {
                    continue;
                }
            }

            // skip the electrodes with high amplitude
            for (int electrode = 1; electrode < nElectrodes; electrode++) {
                if (electrodeMap.isDisconnected(electrode)) {
                    continue;
                }

                if (amplitudeList[electrode] > sigmas[electrode] * minSpikeFindingSigma) {
                    skipElectrode[electrode] = true;
                    nCleanedElectrodes++;
                }
            }

            // save the neuron
            final int id = neuronFile.addExactNeuron(
                    currentElectrode, exactNeuronSpikeTimes, N);
            neuronIndex.add(neuron);
            neuronID.add(id);
            neuronsToRemoveList.add(id);

            //write a report
            System.out.println("Neuron " + id + " added for removal. " +
                    neuronsToRemoveList.size() + " neurons waiting removal.");

            String ID = StringUtil.format(id, 0, 3);
            String s = datasetFolder + File.separator + "report" + File.separator;

            html.println("<tr>");
            html.println("<td>" + id);
            int[] idList = neuronFile.getIDList();
            int n = 0;
            for (int i = 0; i < idList.length; i++) {
                if (neuronFile.getElectrode(idList[i]) == currentElectrode) {
                    n++;
                }
            }
            html.println("<td>" + StringUtil.format(currentElectrode, 0, 3) +
                    (char) ('A' + n - 1));
            html.println("<td>" + StringUtil.toString(electrodes));
            html.println("<td>" + PlotUtil.getColorName(neuron));
            html.println("<td>" + N);
            html.println("<td>" + StringUtil.format(contaminationIndex[neuron], 3));
            html.println("<td>" + nBadSpikes[neuron]);
            html.println("<td>" +
                    StringUtil.format(amplitudeList[currentElectrode], 0));
            html.println("<td><a href=\"" + ID + "-auto.png\">Autocorrelation</a>");
            html.println("<td><a href=\"" + ID + "-amp.png\">Amplitude Hist</a>");
            html.println("<td><a href=\"" + ID + "-ei.png\">EI</a>");
            html.println("<td>");
            html.println("<a href=\"" + ID + "-clust12.png\">1-2</a>");
            html.println("<a href=\"" + ID + "-clust13.png\">1-3</a>");
            html.println("<a href=\"" + ID + "-clust23.png\">2-3</a>");

            // save the plots
//			userInterface.clusteringPlots[0].autoscale();
            userInterface.clusteringPlots[0].saveAsPNG(s + ID + "-clust12.png"
            /*, imageSize, imageSize*/);
//			userInterface.clusteringPlots[1].autoscale();
            userInterface.clusteringPlots[1].saveAsPNG(s + ID +
                    "-clust13.png"
                    /*, imageSize,
                  imageSize*/);
//			userInterface.clusteringPlots[2].autoscale();
            userInterface.clusteringPlots[2].saveAsPNG(s + ID +
                    "-clust23.png"
                    /*, imageSize,
                  imageSize*/);
            //          autocorrelationPanel[neuron].autoscale();
            autocorrelationPanel[neuron].saveAsPNG(s + ID + "-auto.png" /*, imageSize,
                                                                    imageSize / 2*/);
            //         amplitudePanel[neuron].autoscale();
            amplitudePanel[neuron].saveAsPNG(
                    s + ID + "-amp.png" /*, imageSize, imageSize / 2*/);

            PhysiologicalImagePanel p = new PhysiologicalImagePanel(
                    averageEI[neuron], sigmas, minSpikeFindingSigma, electrodeMap,
                    currentElectrode);
            GraphicsIO.saveAsPNG(p, new File(s + ID + "-ei.png"), imageSize,
                    imageSize / 2);

            canExitSNF = true;
        }

        if (neuronID.size() != 0) {
            ClusteringModelFile.Model model = userInterface.getModel();
            model.extractionID = currentExtractionID;
            currentExtractionID++;
            model.neuronIndex = neuronIndex.toArray();
            model.neuronID = neuronID.toArray();
            modelFile.addExtraction(model);
        }
    }


    /**
     * Returns true if the application can exit and false if it cannot exit because
     * there is unsaved data
     */
    public boolean canClose() {
        if (canExitSNF) {
            System.out.println("Saving amplitude histograms...");
            try {
                saveAmplitudeHistograms(amplitudeHistogram,
                        amplitudeHistogramsFileName);
            } catch (IOException e) {
                Vision.reportException(
                        "The amplitude histograms could not be saved", e);
            }

            try {
                saveStatus();
            } catch (IOException e) {
                Vision.reportException("The status file could not be saved", e);
            }

            try {
                html.close();
                neuronFile.close();
                imagingFile.close();
                rawDataFile.close();
                originalRawDataFile.close();
                modelFile.close();
            } catch (IOException e) {
                Vision.reportException(
                        "One of the opened files could not be closed", e);
            }
            Vision.getInstance().getCalculationManager().calculationDone();
            return true;
        } else {
            System.out.println(
                    "SNF cannot exit now, an unstoppable process is going on.");
            return false;
        }
    }


    private void manualNeuronFinding() {
        /*
                 // create the terminal
                 JTextArea outputTextArea = new JTextArea();
                 terminal = new JScrollPane(outputTextArea);
                 System.setOut(new UIPrintStream(outputStream, outputTextArea));
                 System.setErr(System.out);
         */

        eventQueue = new VisionEventQueue();
        //    Toolkit.getDefaultToolkit().getSystemEventQueue().push(eventQueue);

        // create user interfaces
        userInterfaceList = new UserInterface[3];
        float[][] pJX = new float[maxClusters * 2][_maxSpikes];
        userInterfaceList[0] = new EMUserInterface(2, pJX);
        userInterfaceList[1] = new EMUserInterface(1, pJX);
        userInterfaceList[2] = new ManualUserInterface();

        // make the parameters dialog
        ParametersTable paramsTable = new ParametersTable();
        paramsTable.addParameter(maxAmplitude);
        paramsTable.addParameter(amplitudeScale);
        paramsTable.addParameter(showEI);
        paramsTable.addParameter(showRawData);
        paramsTable.addParameter(autoOutline);
        paramsTable.addParameter(autoLabels);
        paramsTable.addParameter(showAllClusteringPlots);
        paramsTable.addParameter(showGaussianFit);
        paramsTable.addParameter(nSamplesParameter);
//		paramsTable.addParameter(truncateDataParameter);
        optionsDialog = new ParametersDialog(mainFrame, "Options", paramsTable);


//		mainFrame = Vision.getInstance().getMainFrame();
        this.menuList = makeMenuBar();

        splitPane1 = new JSplitPane(JSplitPane.VERTICAL_SPLIT);
        splitPane2 = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
        this.add(splitPane2);
        setUI(0);

//		mainFrame.add(splitPane2);
//		mainFrame.setBounds(200, 0, 1280 - 200, 1000);
//		mainFrame.setVisible(true);
        Vision app = Vision.getInstance();

        mainFrame = app.createFrame(
                this, null, menuList, "Serial Neuron Finder");

        myGlassPane = new MyGlassPane(mainFrame);
        mainFrame.setGlassPane(myGlassPane);


        nextElectrode(CalculateMode.SWITCH_ELECTRODE);

//		final SynchronizationObject syncObject = new SynchronizationObject();
//		syncObject.setWorking();
//		syncObject.waitUntilDone();
    }


    private void setUI(int i) {
        userInterface = userInterfaceList[i];
        eventQueue.setKeyListener(userInterface);

        splitPane1.setLeftComponent(userInterface.rightPanel);
//		splitPane1.setRightComponent(terminal);
        splitPane2.setRightComponent(splitPane1);
        splitPane2.setLeftComponent(userInterface.leftPanel);

        splitPane1.setDividerLocation(800);
        splitPane2.setDividerLocation(400);
    }


    private void endCalculation() {
        if (canClose()) {
            Vision.getInstance().removeFrame(mainFrame);
        }
    }


    public static void saveAmplitudeHistograms(int[][] hist, String fileName) throws
    IOException {

        DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new
                FileOutputStream(fileName), 1024));
        dos.writeInt(hist.length);
        dos.writeInt(hist[0].length);
        for (int electrode = 0; electrode < hist.length; electrode++) {
            for (int i = 0; i < hist[electrode].length; i++) {
                dos.writeInt(hist[electrode][i]);
            }
        }
        dos.close();
    }


    public void skipElectrode() {
        if (canExitSNF) {
            canExitSNF = false;
            skipElectrode[currentElectrode] = true;
            nextElectrode(CalculateMode.SWITCH_ELECTRODE);
        }
    }


    public void skipElectrodeForever() {
        if (canExitSNF) {
            canExitSNF = false;
            skipElectrodeForever[currentElectrode] = true;
            nextElectrode(CalculateMode.SWITCH_ELECTRODE);
        }
    }


    public void leaveElectrode() {
        if (canExitSNF) {
            canExitSNF = false;
            nextElectrode(CalculateMode.SWITCH_ELECTRODE);
        }
    }


    private JMenu makeActionsMenu() {
        JMenu menu = new JMenu("SNF: Actions");
        menu.setMnemonic('A');
        JMenuItem item;

        item = new JMenuItem("Analyze full data");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                nSamplesToAnalyze = nSamples;
                nextElectrode(CalculateMode.RELOAD_DATA);
            }
        });
        menu.add(item);

        menu.add(new JSeparator());

        item = new JMenuItem("Save Status");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                try {
                    saveStatus();
                } catch (IOException e) {
                    Vision.reportException("The status file could not be saved",
                            e);
                }
            }
        });
        menu.add(item);

        menu.add(new JSeparator());

        item = new JMenuItem("Skip This Electrode");
        item.setAccelerator(KeyStroke.getKeyStroke(
                KeyEvent.VK_S, ActionEvent.CTRL_MASK /* + ActionEvent.SHIFT_MASK*/));
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                skipElectrode();
            }
        });
        menu.add(item);

        item = new JMenuItem("Skip This Electrode Pattern");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                if (canExitSNF) {
                    canExitSNF = false;
                    for (int i = 0; i < electrodes.length; i++) {
                        skipElectrode[electrodes[i]] = true;
                    }
                    nextElectrode(CalculateMode.SWITCH_ELECTRODE);
                }
            }
        });
        menu.add(item);

        item = new JMenuItem("Skip This Electrode Forever");
        item.setAccelerator(KeyStroke.getKeyStroke(
                KeyEvent.VK_F, ActionEvent.CTRL_MASK /* + ActionEvent.SHIFT_MASK*/));
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                skipElectrodeForever();
            }
        });
        menu.add(item);

        item = new JMenuItem("Skip This Electrode Pattern Forever");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                if (canExitSNF) {
                    canExitSNF = false;
                    for (int i = 0; i < electrodes.length; i++) {
                        skipElectrodeForever[electrodes[i]] = true;
                    }
                    nextElectrode(CalculateMode.SWITCH_ELECTRODE);
                }
            }
        });
        menu.add(item);

        item = new JMenuItem("Leave This Electrode...");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                leaveElectrode();
            }
        });
        menu.add(item);

        menu.add(new JSeparator());

        item = new JMenuItem("Do Neuron Cleaning");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ex) {
                String message =
                    "There still are unprocessed electrodes, they may contain good neurons.\n" +
                    " Are you sure you want to start the removal now?";
                int r = JOptionPane.showOptionDialog(
                        mainFrame, message, "Removal Needed", JOptionPane.YES_NO_OPTION,
                        JOptionPane.QUESTION_MESSAGE, null, null, null);

                if (r == JOptionPane.YES_OPTION) {
                    setMenuBarEnabled(false);
                    cleaning();
                }
            }
        });
        menu.add(item);
        /*
                item = new JMenuItem("Subtract EIs");
                item.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent ex) {
                        JPanel p = new JPanel(new GridLayout(0, 1));
                        JComboBox box1 = new JComboBox();
                        JComboBox box2 = new JComboBox();
                        for (int i = 0; i < userInterface.getNeuronsCount(); i++) {
                            box1.addItem(new Integer(i));
                            box2.addItem(new Integer(i));
                        }
                        p.add(box1);
                        p.add(box2);
                        int result = JOptionPane.showOptionDialog(
             mainFrame, p, "Choose clusters", JOptionPane.OK_CANCEL_OPTION,
                            JOptionPane.QUESTION_MESSAGE, null, null, null);
                        if (result == JOptionPane.OK_OPTION) {
                            int n1 = ( (Integer) box1.getSelectedItem()).intValue();
                            int n2 = ( (Integer) box2.getSelectedItem()).intValue();
                            calculateEI(n1);
                            calculateEI(n2);
                            float[][] ei1 = averageEI[n1][0];
                            float[][] ei2 = averageEI[n2][0];
                            float[][] ei = new float[nElectrodes][ei1[0].length];
             for (int electrode = 0; electrode < nElectrodes; electrode++) {
                                for (int i = 0; i < ei[0].length; i++) {
             ei[electrode][i] = ei1[electrode][i] - ei2[electrode][i];
                                }
                            }
                menu.add(item);
         */

        menu.add(new JSeparator());
        item = new JMenuItem("Exit");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                endCalculation();
            }
        });
        menu.add(item);

        return menu;
    }


    JMenu[] makeMenuBar() {
        viewMenu = new JMenu("SNF: View");
        viewMenu.setMnemonic('V');

        return new JMenu[] {makeActionsMenu(), viewMenu, makeOptionsMenu(),
                makeThresholdMenu()};
    }


    private JMenu makeOptionsMenu() {
        JMenu menu = new JMenu("SNF: Options");
        menu.setMnemonic('O');

        ActionListener uiListener = new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                int n = Integer.parseInt( ( (JComponent) e.getSource()).getName());
                setUI(n);
                userInterface.dataUpdated();

                mainFrame.validate();
                mainFrame.repaint();

                userInterface.updateNeuronsInfo();
                userInterface.updateView(true, true);
            }
        };

        ActionListener patternListener = new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                int n = Integer.parseInt( ( (JComponent) e.getSource()).getName());
                electrodeUsage = ElectrodeUsage.values()[n];
                electrodes = electrodeMap.getAdjacentsTo(currentElectrode,
                        electrodeUsage.ordinal());
                nextElectrode(CalculateMode.RELOAD_DATA);
            }
        };

        String[] s = {
                "Two Gaussians/Cluster Expectation Maximization",
                "One Gaussian/Cluster Expectation Maximization",
                "Manual Boxing"};
        ButtonGroup uiGroup = new ButtonGroup();
        for (int i = 0; i < s.length; i++) {
            JRadioButtonMenuItem button = new JRadioButtonMenuItem(s[i], i == 0);
            button.addActionListener(uiListener);
            button.setName("" + i);
            uiGroup.add(button);
            menu.add(button);
        }
        menu.add(new JSeparator());

        String[] pattern = {"Single electrode pattern",
        "Electrode + 6 neighbors pattern"};
        ButtonGroup patternGroup = new ButtonGroup();
        for (int i = 0; i < pattern.length; i++) {
            JRadioButtonMenuItem button = new JRadioButtonMenuItem(pattern[i],
                    ElectrodeUsage.values()[i] == electrodeUsage);
            button.addActionListener(patternListener);
            button.setName("" + i);
            patternGroup.add(button);
            menu.add(button);
        }
        menu.add(new JSeparator());

        JMenuItem item = new JMenuItem("Show Neuron Statistics");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                PlotPanel p = new PlotPanel();
                p.setAxisVisible(false);
                p.setRange(electrodeMap.getBounds());

                double r = 1;

                ScatterPlot sp = new ScatterPlot();
                for (int electrode = 0; electrode < nElectrodes; electrode++) {
                    sp.add(electrodeMap.getXPosition(electrode),
                            electrodeMap.getYPosition(electrode), r);
                }

                ScatterPlotStyle style = new ScatterPlotStyle();
                style.setSymbolType(SymbolType.NONE);
                style.setErrorSymbolType(SymbolType.DISK);
                p.addData(sp, style);

                PlotUtil.showData("Neuron Locations", p);
            }
        });
        menu.add(item);

        item = new JMenuItem("Show Electrode Statistics");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                try {
                    int[] id = modelFile.getNeuronsForCleanning(cleaningLevel);
                    boolean[] neuronExtracted = new boolean[nElectrodes];
                    for (int i = 0; i < id.length; i++) {
                        neuronExtracted[neuronFile.getElectrode(id[i])] = true;
                    }

                    ElectrodeStatistics s = new ElectrodeStatistics(
                            electrodeMap, skipElectrode, skipElectrodeForever,
                            neuronExtracted);

                    JDialog f = new JDialog( (JFrame)null, "Electrode Status", true);
                    f.add(s);
                    if (nElectrodes == 65) {
                        f.setSize(450 + 8, 450 + 27 + 60);
                    } else if (nElectrodes == 513) {
                        f.setSize(900 + 8, 450 + 27 + 60);
                    }
                    f.setLocationRelativeTo(mainFrame);
                    f.setVisible(true);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        });
        menu.add(item);

        item = new JMenuItem("Show Spike Time Dependence");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                /*
                                 ScatterPlot sp = new ScatterPlot();
                                 int nNeurons = userInterface.getNeuronsCount();
                                 for (int neuron = 0; neuron < nNeurons; neuron++) {
                    int n = getNeuronSpikeTimes(neuron);
                    for (int sIndex = 0; sIndex < n; sIndex++) {
                        sp.add(exactNeuronSpikeTimes[sIndex], -neuron * 10);
                    }
                                 }
                                 ScatterPlotStyle style = new ScatterPlotStyle();
                                 style.setSymbolType(SymbolType.VERTICAL_LINE);
                                 style.setSymbolSize(50);
                                 PlotPanel p = new PlotPanel();
                                 p.addData(sp, style);
                                 p.autoscale();
                                 PlotUtilities.showData("Spike Time Dependence", p);
                 */

                double dt = DoubleDialog.showDoubleInputDialog(mainFrame, "Binning",
                        "Type the spike rate histogram binning in seconds", 1, 5e-5, 1e10);

                JFrame f = new JFrame("Spike Rate Time Dependence");
                f.getContentPane().setLayout(new GridLayout(0, 1));

                int nNeurons = userInterface.getNeuronsCount();
                for (int n = 0; n < nNeurons; n++) {
                    DoubleHistogram h = new DoubleHistogram("", 0, nSamples / 20000.0, dt);
                    double[] t = getNeuronSpikeTimes1(n);
                    for (int i = 0; i < t.length; i++) {
                        h.fill(t[i] / 20000, 1);
                    }
                    h.scale(1.0 / dt);

                    PlotPanel p = new PlotPanel();
                    p.setLabels("", PlotUtil.getColorName(n));
//					ArrayList l = new ArrayList();
//l.add(PlotUtilities.getColorName(i) + " - " +
//					PlotUtilities.getColorName(j) + " : " +
//					StringUtil.format(c * 100, 1) + "%");
//					p.setAdditionalLegend(l);
                    p.addData(h, new HistogramStyle());
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
                        ArrayList<String> l = new ArrayList<String>();
                        l.add(PlotUtil.getColorName(i) + " - " +
                                PlotUtil.getColorName(j) + " : " +
                                StringUtil.format(c * 100, 1) + "%");
                        p.setAdditionalLegend(l);
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

        item = new JMenuItem("PCs as a function of time");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ex) {
                JDialog f = new JDialog( (JFrame)null, "PCs as a function of time", true);
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

        item = new JMenuItem("Show Eigenvectors");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ex) {
                JDialog f = new JDialog( (JFrame)null, "Eigenvectors", false);
                f.getContentPane().setLayout(new GridLayout(0, 1));

                ScatterPlotStyle style = new ScatterPlotStyle();
                style.setConnectingPoints(true);
                style.setConnectionPeriod(nPoints + 0);
                for (int d = 0; d < nDimensions; d++) {
                    ScatterPlot sp = new ScatterPlot();
                    double[] eigenvector = pca.getEigenVector(d);
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

                if (nSamplesParameter.getValue() != nSamplesToAnalyze) {
                    nSamplesToAnalyze = nSamplesParameter.getValue();
                    nextElectrode(CalculateMode.RELOAD_DATA);
//					} else if (truncateDataParameter.getValue() != TRUNCATE_DATA) {
//	TRUNCATE_DATA = truncateDataParameter.getValue();
//nextElectrode(NEW_SPIKE_FINDING);
                } else {
                    userInterface.updateNeuronsInfo();
                    userInterface.updateView(true, true);
                }
            }
        });
        menu.add(item);

        return menu;
    }


    private JMenu makeThresholdMenu() {
        JMenu menu = new JMenu("SNF: Threshold");
        menu.setMnemonic('T');
        JMenuItem item;

        ActionListener listener1 = new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                floatingSigma = DoubleDialog.showDoubleInputDialog(
                        mainFrame, "Threshold",
                        "Input the Threshold Value In Standard Deviations:",
                        floatingSigma, minSpikeFindingSigma, 1e10);
                userInterface.lock();
                nextElectrode(CalculateMode.NEW_SPIKE_FINDING);
            }
        };
        item = new JMenuItem("In Standard Deviations...");
        item.addActionListener(listener1);
        menu.add(item);

        ActionListener listener2 = new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                double A = DoubleDialog.showDoubleInputDialog(
                        mainFrame, "Threshold",
                        "Input the Threshold Value In ADC Counts:",
                        sigmas[currentElectrode] * floatingSigma,
                        sigmas[currentElectrode] * minSpikeFindingSigma, 1e10);
                floatingSigma = A / sigmas[currentElectrode];
                userInterface.lock();
                nextElectrode(CalculateMode.NEW_SPIKE_FINDING);
            }
        };
        item = new JMenuItem("In ADC Counts...");
        item.addActionListener(listener2);
        menu.add(item);

        return menu;
    }


    private void calculateEI(int neuron) {
        double[] t = getNeuronSpikeTimes1(neuron);
        try { 
            averageEI[neuron] = WaveformCalculator1.calculateEI(
                    t, t.length, 100, nlPointsEI, nrPointsEI, originalRawDataFile);
        } catch (IOException e) {
            e.printStackTrace();
        }

        // remove the means form the EIs
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


    ActionListener viewListener = new ActionListener() {
        public void actionPerformed(ActionEvent e) {
            JComponent source = (JComponent) e.getSource();
            int neuron = Integer.parseInt(source.getName());
            showEI(neuron);
        }
    };

    /*
        WhiteNoiseMovie m1;
        ActionListener staListener = new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                JComponent source = (JComponent) e.getSource();
                int neuron = Integer.parseInt(source.getName());
//            try {
                    if (m1 == null) {
                        m1 = WhiteNoiseMovie.load(
                            "f:\\2003-08-16-0\\data016\\data016.movie");
                    }
                    int firstTTL = neuronFile.getTTLTimes()[0];
                    STACalculator sta = new STACalculator(40, firstTTL, 0.05,
                        new BufferedMovie(m1));
                    final int n = getNeuronSpikeTimes(neuron);
                    for (int i = 0; i < n; i++) {
                        sta.addSpike( (int) Math.round(exactNeuronSpikeTimes[i]));
                    }
                    sta.finish();
                    VisionUtilities.showSTA(
                        "" + PlotUtilities.getColorName(neuron), sta.getSTA(), null);
//            } catch (IOException ex) {
//                ex.printStackTrace();
//            }
            }
        };
     */

    class ElectrodeDialog
    extends JDialog {

        PlotPanel panel;
        DoubleHistogram h;
        JPanel northPanel, buttonPanel;
        final JButton acceptButton, skipButton, skipCompletellyButton;
        final JSpinner electrodeSpinBox;
        private JLabel statusLabel;
        ArrayList legend = new ArrayList();


        public ElectrodeDialog() {
            super( (JFrame)null, "Choose the next electrode", true);
            setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE);

            JMenuBar menuBar = new JMenuBar();
            JMenu menu = new JMenu("Actions");
            JMenuItem item;

            item = new JMenuItem("Do Cleaning");
            item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    int selection = JOptionPane.showOptionDialog(
                            mainFrame, "Do you really want to clean now?", "Cleanning ?",
                            JOptionPane.YES_NO_OPTION,
                            JOptionPane.QUESTION_MESSAGE, null, null, null);
                    if (selection == JOptionPane.YES_OPTION) {
                        electrodeSpinBox.setValue(new Integer( -1));
                        setVisible(false);
                        cleaning();
                    }
                }
            });
            menu.add(item);
            menu.add(new JSeparator());

            item = new JMenuItem("Exit SNF");
            item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    canExitSNF = true;
                    endCalculation();
                }
            });
            menu.add(item);
            menuBar.add(menu);
            setJMenuBar(menuBar);

            panel = new PlotPanel();
            panel.setLabels("Spike Amplitude (ADC)", "Spike Count");
//			panel.setBackground(Color.lightGray);
h = new DoubleHistogram("", 0, amplitudeHistogram[1].length, 1); ;

acceptButton = new JButton("Accept");
acceptButton.setMnemonic('A');
acceptButton.addActionListener(new ActionListener() {
    public void actionPerformed(ActionEvent e) {
        setVisible(false);
    }
});
skipButton = new JButton("Skip");
skipButton.setMnemonic('S');
skipButton.addActionListener(new ActionListener() {
    public void actionPerformed(ActionEvent e) {
        int electrode = ( (Integer) electrodeSpinBox.getValue()).intValue();
        skipElectrode[electrode] = true;

        int nextElectrode = getBestElectrode();
        if (nextElectrode == -1) {
            electrodeSpinBox.setValue(new Integer( -1));
            setVisible(false);
            doCleaning();
            return;
        } else {
            electrodeSpinBox.setValue(new Integer(nextElectrode));
        }
    }
});
skipCompletellyButton = new JButton("Skip Forever");
skipCompletellyButton.addActionListener(new ActionListener() {
    public void actionPerformed(ActionEvent e) {
        int electrode = ( (Integer) electrodeSpinBox.getValue()).intValue();
        skipElectrodeForever[electrode] = true;

        int nextElectrode = getBestElectrode();
        if (nextElectrode == -1) {
            electrodeSpinBox.setValue(new Integer( -1));
            setVisible(false);
            doCleaning();
            return;
        } else {
            electrodeSpinBox.setValue(new Integer(nextElectrode));
        }
    }
});
buttonPanel = new JPanel();
buttonPanel.add(acceptButton);
buttonPanel.add(skipButton);
buttonPanel.add(skipCompletellyButton);

electrodeSpinBox = new JSpinner(
        new SpinnerNumberModel(1, 1, nElectrodes - 1, 1));
electrodeSpinBox.addChangeListener(new ChangeListener() {
    public void stateChanged(ChangeEvent e) {
        int electrode = ( (Integer) electrodeSpinBox.getValue()).intValue();
        if (electrode == -1) {
            return;
        }

        h.clear();
        for (int i = 0; i < amplitudeHistogram[1].length; i++) {
            h.setBin(i, amplitudeHistogram[electrode][i]);
        }

        panel.removeAllData();
        panel.addData(h, new HistogramStyle());
        panel.autoscale();
        panel.setXRange(0, maxAmplitude.getValue());
        panel.setPreferredSize(new Dimension(400, 400));

        // estimate number od spike on this electrode
        int threshold =
            (int) Math.floor(sigmas[electrode] * minSpikeFindingSigma) -
            1;
        int nSpikes = 0;
        for (int i = threshold; i < amplitudeHistogram[electrode].length;
        i++) {
            nSpikes += amplitudeHistogram[electrode][i];
        }

        // create the legend showing the electrode statistics
        int nForever = MathUtil.countValues(true,
                skipElectrodeForever);
        int nSkipped = MathUtil.countValues(true, skipElectrode);
        int nDisconnected = electrodeMap.countDisconnectedElectrodes();
        legend.clear();
        legend.add("nSpikes ~ " + nSpikes / 1000 + "k");
        legend.add("Electrodes Left: " +
                (nElectrodes - 1 - nForever - nSkipped - nDisconnected));
        legend.add("Skipped: " + nSkipped);
        legend.add("Skipped Forever: " + nForever);
        legend.add("Disconnected: " + nDisconnected);
        panel.setAdditionalLegend(legend);

        if (electrodeMap.isDisconnected(electrode)) {
            acceptButton.setEnabled(false);
            statusLabel.setText(" is disconnected");
        } else if (skipElectrodeForever[electrode]) {
            acceptButton.setEnabled(true);
            statusLabel.setText(" is skipped forever");
        } else if (isIgnored[electrode]) {
            acceptButton.setEnabled(false);
            statusLabel.setText(" is ignored");
        } else if (skipElectrode[electrode]) {
            acceptButton.setEnabled(true);
            statusLabel.setText(" is skipped");
        } else {
            acceptButton.setEnabled(true);
            statusLabel.setText("");
        }
    }
});
northPanel = new JPanel();
northPanel.add(new JLabel("Electrode "));
northPanel.add(electrodeSpinBox);
statusLabel = new JLabel("");
northPanel.add(statusLabel);

this.add(northPanel, BorderLayout.NORTH);
this.add(panel, BorderLayout.CENTER);
this.add(buttonPanel, BorderLayout.SOUTH);
        }


        private void doCleaning() {
            String message =
                "\n Sorry, all the elecrodes are used up!" +
                "\n You cannot continue unless you remove the " +
                neuronsToRemoveList.size() + " waiting neurons.";
            System.out.println(message);
            JOptionPane.showMessageDialog(
                    mainFrame, message, "Removal Needed", JOptionPane.ERROR_MESSAGE);
            cleaning();
        }


        int getBestElectrode() {
            int N = 1000;
            int maxAmplitude = Integer.MIN_VALUE;
            int maxAmplitudeElectrode = -1;

            for (int electrode = 1; electrode < nElectrodes; electrode++) {
                if (electrodeMap.isDisconnected(electrode) ||
                        skipElectrode[electrode] ||
                        skipElectrodeForever[electrode] ||
                        isIgnored[electrode]) {
                    continue;
                }

                int sum = 0;
                int amplitude = amplitudeHistogram[electrode].length - 1;
                while (sum < N && amplitude >= 0) {
                    sum += amplitudeHistogram[electrode][amplitude];
                    amplitude--;
                }
                if (amplitude > maxAmplitude) {
                    maxAmplitude = amplitude;
                    maxAmplitudeElectrode = electrode;
                }
            }

            return maxAmplitudeElectrode;
        }


        public int geNextElectrode() {
            final int bestElectrode = getBestElectrode();
            if (bestElectrode == -1) {
                doCleaning();
                return -1;
            } else {
                electrodeSpinBox.setValue(new Integer(bestElectrode));
                setSize(400, 400);
                setLocationRelativeTo(mainFrame);
                setVisible(true);

                return (Integer) electrodeSpinBox.getValue();
            }
        }
    }


    public static int[][] loadAmplitudeHistograms(String fileName) throws
    IOException {

        DataInputStream dis = new DataInputStream(
                new BufferedInputStream(new FileInputStream(fileName), 1024));
        int nElectrodes = dis.readInt();
        int hSize = dis.readInt();
        int[][] amplitudeHistogram = new int[nElectrodes][hSize];

        for (int electrode = 0; electrode < nElectrodes; electrode++) {
            for (int i = 0; i < hSize; i++) {
                amplitudeHistogram[electrode][i] = dis.readInt();
            }
        }

        dis.close();
        return amplitudeHistogram;
    }


    public static class ElectrodeStatistics
    extends JPanel {

        ElectrodeMap map;
        int nElectrodes;
        boolean[] isSkipped, isSkippedForever, neuronExtracted;
        int currentElectrode = -1;
        FontMetrics metrics;


        public ElectrodeStatistics(ElectrodeMap map, boolean[] isSkipped,
                boolean[] isSkippedForever, boolean[] neuronExtracted) {

            this.map = map;
            this.nElectrodes = map.getNumberOfElectrodes();
            this.isSkipped = isSkipped;
            this.isSkippedForever = isSkippedForever;
            this.neuronExtracted = neuronExtracted;
        }


        public void drawColored(double x, double y, double r, Color color, String s,
                Graphics g) {
            g.setColor(color);
            g.fillOval( (int) Math.round(x - r), (int) Math.round(y - r),
                    (int) Math.round(2 * r), (int) Math.round(2 * r));
            g.setColor(Color.black);
            g.drawOval( (int) Math.round(x - r), (int) Math.round(y - r),
                    (int) Math.round(2 * r), (int) Math.round(2 * r));

            if (color == Color.black) {
                g.setColor(Color.white);
            }

            int w = metrics.stringWidth(s);
            g.drawString(s, (int) (x - w / 2.0), (int) (y + metrics.getAscent() / 2.0));

            g.setColor(Color.black);
        }


        final synchronized public void paintComponent(Graphics _g) {
            Graphics2D g = (Graphics2D) _g;
            g.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                    RenderingHints.VALUE_ANTIALIAS_ON);

            double pitch = map.getPitch();
            double[] bounds = map.getBounds();

            double x1 = 0;
            double x2 = getWidth();
            double y1 = 60;
            double y2 = getHeight();

            double x1p = bounds[0];
            double x2p = bounds[1];
            double y1p = bounds[2];
            double y2p = bounds[3];

            double ax = (x2 - x1) / (x2p - x1p);
            double bx = x1 - ax * x1p;
            double ay = (y2 - y1) / (y2p - y1p);
            double by = y1 - ay * y1p;

            double r = 0.95 * Math.min(ax * pitch, ay * pitch) / 2.0;

            g.clearRect(0, 0, getWidth(), getHeight());
            g.setFont(new Font("Dialog", Font.PLAIN, 11));
            metrics = g.getFontMetrics();

            for (int electrode = 0; electrode < nElectrodes; electrode++) {
                double x = ax * map.getXPosition(electrode) + bx;
                double y = y1 + y2 - (ay * map.getYPosition(electrode) + by);

                if (map.isDisconnected(electrode)) {
                    drawColored(x, y, r, Color.white, "" + electrode, g);
                } else if (neuronExtracted[electrode]) {
                    drawColored(x, y, r, Color.green, "" + electrode, g);
                } else if (isSkippedForever[electrode]) {
                    drawColored(x, y, r, Color.black, "" + electrode, g);
                } else if (isSkipped[electrode]) {
                    drawColored(x, y, r, Color.lightGray, "" + electrode, g);
                } else {
                    drawColored(x, y, r, Color.red, "" + electrode, g);
                }

                if (electrode == currentElectrode) {
                    g.setColor(Color.black);
                    g.drawOval( (int) Math.round(x - 1.1 * r),
                            (int) Math.round(y - 1.1 * r),
                            (int) Math.round(2 * 1.1 * r),
                            (int) Math.round(2 * 1.1 * r));

                }
            }

            g.drawLine( (int) Math.round(x1), (int) Math.round(y1),
                    (int) Math.round(x2), (int) Math.round(y1));
            g.setColor(Color.black);
            g.setFont(new Font("normal", Font.BOLD, 14));
            r = 10;
            int dx = 160;
            double x, y;

            x = 20;
            y = 1.5 * r;
            drawColored(x, y, r, Color.white, "", g);
            g.drawString("Disconnected", (int) (x + 2 * r), (int) (y + 0.5 * r));

            x = 20;
            y = 4 * r;
            drawColored(x, y, r, Color.lightGray, "", g);
            g.drawString("Skipped", (int) (x + 2 * r), (int) (y + 0.5 * r));

            x = 20 + dx;
            y = 1.5 * r;
            drawColored(x, y, r, Color.black, "", g);
            g.drawString("Skipped Forever", (int) (x + 2 * r), (int) (y + 0.5 * r));

            x = 20 + dx;
            y = 4 * r;
            drawColored(x, y, r, Color.green, "", g);
            g.drawString("Neuron Extracted", (int) (x + 2 * r), (int) (y + 0.5 * r));

            x = 20 + 2 * dx;
            y = 1.5 * r;
            drawColored(x, y, r, Color.red, "", g);
            g.drawString("Unchecked", (int) (x + 2 * r), (int) (y + 0.5 * r));

//			if (currentElectrode != -1) {
//			g.drawString("Electrode: " + currentElectrode, (int) (20 + 3 * dx),
//			(int) (y + 2 * r));
//			}

            g.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                    RenderingHints.VALUE_ANTIALIAS_OFF);
        }
    }


    public static double getCorrelation(double[] tList1, int n1, double[] tList2, int n2) throws
    IOException {

        DoubleHistogram ccH = CrossCorrelationCalculator.getCrossCorrelationHistogram(
                tList1, tList2, 1, 400);
        final int n = ccH.getBinCount();
        double maxV = Double.NEGATIVE_INFINITY;

        for (int k = 0; k < n - 2 * coincidenceTime; k++) {
            double v = 0;
            for (int m = 0; m < 2 * coincidenceTime; m++) {
                v += ccH.getBin(k + m);
            }
            if (v > maxV) {
                maxV = v;
            }
        }

        return maxV / Math.min(tList1.length, tList2.length);
    }


    public static double getMaxCorrelation(double[] tList1, int n1, NeuronFile neuronFile) throws
    IOException {

        double maxP = Double.NEGATIVE_INFINITY;
        int[] idList = neuronFile.getIDList();

        for (int i = 0; i < idList.length; i++) {
            double[] tList2 = neuronFile.getExactSpikeTimes(idList[i]);
            double p = getCorrelation(tList1, n1, tList2, tList2.length);
            if (p > maxP) {
                maxP = p;
            }
        }

        return maxP;
    }


}
