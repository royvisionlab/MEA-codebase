package edu.ucsc.neurobiology.vision.neuronviewer;

import java.io.*;
import java.lang.reflect.*;
import java.util.*;
import java.util.List;

import static java.awt.Color.*;

import java.awt.*;
import java.awt.event.*;
import java.awt.geom.*;

import javax.swing.*;
import javax.swing.FocusManager;
import javax.swing.event.*;
import javax.swing.tree.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.dataview.*;
import edu.ucsc.neurobiology.vision.electrodemap.*;
import edu.ucsc.neurobiology.vision.gui.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.io.chunk.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.neuronviewer.PlotMaker.KeybindingBuilder;
import edu.ucsc.neurobiology.vision.parameters.*;
import static edu.ucsc.neurobiology.vision.plot.SymbolType.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public final class NeuronViewer
extends JPanel implements TreeSelectionListener, ChangeListener, ActionListener,
ClassifierListener, Closable {


    private ScatterPlotStyle scatterPlotStyle = new ScatterPlotStyle(
            "Scatter Plot", CIRCLE, 6, black, false, black, 1);

    private static String base = "edu.ucsc.neurobiology.vision.neuronviewer.";
    private static final int TREE_WIDTH = 150;

    private CrossCorrelator correlator;
    private LinkedHashMap<String, String> parametersMap, plotParametersMap;
    private PlotMaker[] plotMakers;
    private List<KeybindingBuilder> keybindingBuilders;
    private ArrayList<CalculationAction> actions;
    private VisionEventQueue eventQueue;
    private ParametersTable plotsParamsTable, actionsParamsTable;
    private HashMap<String, Object> storedData = new HashMap<String, Object>();

    ElectrodeMap electrodeMap = null;
    Point[] arrayCorners;
    int newClassIndex = 0;
    Component mainFrame;
    JSplitPane leftSplitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, null, null);
    JSplitPane rightSplitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, null, null);
    String filePathRoot, datasetName, experimentName;
    Object rightSelection, leftSelection;
    JPanel leftHalfPanel, rightHalfPanel;
    Config configuration;


    // data files
    NeuronFile neuronFile;
    STACollection staCollection, stvCollection;
    ParametersFile paramsFile;
    PhysiologicalImagingFile imgFile;
    GlobalsFile globalsFile;


    public boolean showNeuronsInSubclasses = true;
    public InteractiveTree rightTree, leftTree;
    private boolean internal;



    public Object getStoredData(String name) {
        return storedData.get(name);
    }


    public void storeData(String name, Object data) {
        storedData.put(name, data);
    }


    public NeuronViewer(String fileNameOrPath, boolean internal, Config configuration) throws IOException {
        super(new BorderLayout());
        this.configuration = configuration;
        Object queue = Toolkit.getDefaultToolkit().getSystemEventQueue();
        if (queue == null || !(queue instanceof VisionEventQueue)) {
            // Vision is probably not started
            eventQueue = new VisionEventQueue();
            //FIXME no custom queue available in standalone mode
            // Toolkit.getDefaultToolkit().getSystemEventQueue().push(eventQueue);
        } else {
            // Vision seems to be started
            eventQueue = (VisionEventQueue) queue;
        }

        this.internal = internal;

        correlator = (CrossCorrelator) Vision.getInstance().getFrameOfClass(CrossCorrelator.class);

        parametersMap = configuration.getParameterList("Neuron Viewer 1");
        plotParametersMap = configuration.getParameterList("Neuron Viewer Plots");
        if (parametersMap == null) throw new Error("The group Neuron Viewer 1 does not exist in the config file");

        File ff = new File(fileNameOrPath);
        if (IOUtil.isValidFile(ff)) {
            datasetName = StringUtil.removeExtension(ff.getName());
            experimentName = ff.getParentFile().getParentFile().getName();
            filePathRoot = ff.getParentFile().getAbsolutePath() + File.separator +
            datasetName;
        } else if (IOUtil.isValidFolder(fileNameOrPath)) {
            datasetName = ff.getName();
            experimentName = ff.getParentFile().getName();
            filePathRoot = fileNameOrPath + File.separator + datasetName;
        } else {
            throw new IOException("Bad file: " + fileNameOrPath);
        }

        File neuronF = new File(filePathRoot + ".neurons");
        if (neuronF.exists() && neuronF.canRead()) neuronFile = new NeuronFile(filePathRoot + ".neurons");

        staCollection = STACollection.Factory.fromFilename(filePathRoot + VisionParams.STA_FILE_EXTENSION);
        
        File steF = new File(filePathRoot + ".stv");
        if (steF.exists() && steF.canRead()) stvCollection = new STAFile(filePathRoot + ".stv");

        paramsFile = new ParametersFile(filePathRoot + ".params");

        File gf = new File(filePathRoot + ".globals");
        if (gf.exists() && gf.isFile() && gf.canRead()) {
            globalsFile = new GlobalsFile(gf.getAbsolutePath(), GlobalsFile.READ);
            int createdBy = globalsFile.getCreatedBy();
            if (createdBy == GlobalsFile.UPDATER) {
                System.out.println("Globals file was created through an automatic file update.  Scaling should not be trusted.");
            }
            if (globalsFile.imageCalibrationParamsExists()) {
                GlobalsFile.ImageCalibrationParams calParams = globalsFile.getImageCalibrationParams();
                electrodeMap = ElectrodeMapFactory.getElectrodeMap(calParams.arrayID + (calParams.arrayPart << 16) 
                        + (calParams.arrayNParts << 24));
            }
        }

        File f = new File(filePathRoot + ".ei");
        if (f.exists() && f.isFile() && f.canRead()) {
            imgFile = new PhysiologicalImagingFile(f.getAbsolutePath());

            //prefer value from globals file, if it exists.
            if (electrodeMap == null) {
                electrodeMap = ElectrodeMapFactory.getElectrodeMap(imgFile.getArrayID());
            }
        }




        // try to read the corners file
        f = new File(new File(filePathRoot).getParentFile().getParentFile().
                getAbsolutePath() +
                File.separator + "stimuli" + File.separator + "corners.txt");
        File fImage = new File(filePathRoot + ".ei");
        if (f.exists() && f.canRead() && fImage.exists() && fImage.canRead()) {
            try {
                arrayCorners = readCorners(
                        f, (electrodeMap.getNumberOfElectrodes() == 513) ? 5 : 7);
            } catch (IOException ex) {
                System.err.println("Could not read corners: " + ex.getMessage());
            }
        }

        // read the parameters
        plotsParamsTable   = new ParametersTable(configuration.getGroupNode("Neuron Viewer Plots"));
        actionsParamsTable = new ParametersTable(configuration.getGroupNode("Neuron Viewer Actions"));

        makeTrees(paramsFile.getClassIDs());

        JComponent leftTreeComponent = new JScrollPane(leftTree);
        leftTreeComponent.setPreferredSize(new Dimension(TREE_WIDTH, 1000));
        this.add(leftTreeComponent, BorderLayout.WEST);

        JComponent rightTreeComponent = new JScrollPane(rightTree);
        rightTreeComponent.setPreferredSize(new Dimension(TREE_WIDTH, 1000));
        this.add(rightTreeComponent, BorderLayout.EAST);

        leftHalfPanel = new JPanel(new BorderLayout());
        leftHalfPanel.add(leftTreeComponent, BorderLayout.WEST);
        leftHalfPanel.add(leftSplitPane, BorderLayout.CENTER);
        leftSplitPane.setDividerLocation(Integer.parseInt(parametersMap.get(
        "leftDividerSize")));
        
        rightHalfPanel = new JPanel(new BorderLayout());
        rightHalfPanel.add(rightTreeComponent, BorderLayout.EAST);
        rightHalfPanel.add(rightSplitPane, BorderLayout.CENTER);
        rightSplitPane.setDividerLocation(Integer.parseInt(parametersMap.get(
        "rightDividerSize")));

        setTreeLayout(0);

        if (internal) {
            // plot params table
            ( (DefaultTreeCellRenderer) plotsParamsTable.getTree().getCellRenderer()).setLeafIcon(null);
            ( (DefaultTreeCellRenderer) plotsParamsTable.getTree().getCellRenderer()).setClosedIcon(null);
            ( (DefaultTreeCellRenderer) plotsParamsTable.getTree().getCellRenderer()).setOpenIcon(null);
            plotsParamsTable.expandAll();
            
            // action params table
            ( (DefaultTreeCellRenderer) actionsParamsTable.getTree().getCellRenderer()).setLeafIcon(null);
            ( (DefaultTreeCellRenderer) actionsParamsTable.getTree().getCellRenderer()).setClosedIcon(null);
            ( (DefaultTreeCellRenderer) actionsParamsTable.getTree().getCellRenderer()).setOpenIcon(null);
            actionsParamsTable.expandAll();
            actionsParamsTable.getColumn("Value").setWidth(30);
            
            // make the plots controller
            JPanel controller1 = new JPanel(new BorderLayout());
            JScrollPane scroll = new JScrollPane(plotsParamsTable);
            scroll.setPreferredSize(new Dimension(0, 400));
            scroll.setSize(new Dimension(0, 400));
            controller1.add(scroll, BorderLayout.CENTER);
            JButton button = new JButton("Apply");
            button.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    NeuronViewer.this.updateTrees();
                }
            });
            controller1.add(button, BorderLayout.NORTH);

            // make the actions controller
//			JPanel controller2 = new JPanel(new BorderLayout());
//			scroll = new JScrollPane(actionsParamsTable);
//			scroll.setPreferredSize(new Dimension(0, 400));
//			scroll.setSize(new Dimension(0, 400));
//			controller2.add(scroll, BorderLayout.CENTER);
//			button = new JButton("Apply");
//			button.addActionListener(new ActionListener() {
//			public void actionPerformed(ActionEvent e) {
//			NeuronViewer.this.updateActions();
//			}
//			});
//			controller2.add(button, BorderLayout.NORTH);

//			JTabbedPane pane = new JTabbedPane();
//			pane.add("Plots", controller1);
//			pane.add("Actions", controller2);

            Vision.getInstance().createFrame(
                    this, controller1, new JMenu[] {makeFileMenu(), makeEditMenu(),
                            makeViewMenu()},
                            "Neuron Viewer: " + filePathRoot + ".params");// + " " + experimentName + File.separator + datasetName);
            mainFrame = Vision.getInstance().getMainFrame();
            
            plotsParamsTable.percentFillColumns(new int[]{75, 25});
        } else {
            // make the UI
            final JFrame frame = new JFrame("Neuron Viewer: " + fileNameOrPath);
            
            frame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
            frame.addWindowListener(new WindowAdapter() {
                public void windowClosing(WindowEvent e) {
                    canClose();
                }});
            JMenuBar menuBar = new JMenuBar();
            menuBar.add(makeFileMenu());
            menuBar.add(makeViewMenu());
            frame.setJMenuBar(menuBar);
            frame.getContentPane().add(this);
            frame.setBounds(200, 0, 1680 - 200, 1050 - 30);
            frame.setVisible(true);
            mainFrame = frame;
        }
        
        updateTrees();

        leftTree.setSelectionRow(0);
        rightTree.setSelectionRow(0);		
    }


    ArrayList<KeyStroke> keyStrokes = new ArrayList<KeyStroke>();

    private void makeTrees(LinkedHashMap<Integer,? extends Object> classLabels) {
        DefaultTreeModel treeModel = InteractiveTree.makeModel(classLabels);

        leftTree = new InteractiveTree(treeModel);
        leftTree.addChangeListener(this);
        leftTree.addTreeSelectionListener(this);

        rightTree = new InteractiveTree(treeModel);
        rightTree.addChangeListener(this);
        rightTree.addTreeSelectionListener(this);
        
        // actionsParamsTable is derived from config.xml "Neuron Viewer Actions"
        actions = new ArrayList<CalculationAction>();
        for (int i = 0; i < actionsParamsTable.getParametersCount(); i++) {
            String actionClassName = actionsParamsTable.getParameter(i).getName();
            if (actionClassName.indexOf(".") > 0) continue;
            
            try {
                final CalculationAction action = 
                        (CalculationAction) Class.forName(base + actionClassName).newInstance();
                action.initialize(this);
                actions.add(action);
                
                KeyStroke keyStroke = action.getKeystroke();
                if (keyStroke != null) {
                    keyStrokes.add(keyStroke);
                    eventQueue.addKeyStroke(keyStroke, new ActionListener() {
                        public void actionPerformed(ActionEvent e) {
                            TreePath path = rightTree.getSelectionPath();
                            IntegerList list = new IntegerList();
                            InteractiveTree.getNeuronsInClass(
                                    (DefaultMutableTreeNode) path.getLastPathComponent(),
                                    list, showNeuronsInSubclasses);
                            action.doAction(list, path);
                        }
                    });
                }

                leftTree.addCalculationAction(action);
                rightTree.addCalculationAction(action);
            } catch (Exception ex) {
                System.err.println("Could not create " + actionClassName);
                ex.printStackTrace(Vision.getInstance().getFileStdOut());
            }
        }

        plotMakers = new PlotMaker[plotParametersMap.size()];
        int i = 0;
        for (String plotMakerName : plotParametersMap.keySet()) {
            if (plotMakerName.indexOf(".") >= 0) continue;

            try {
                char lastChar = plotMakerName.charAt(plotMakerName.length() - 1);
                if (Character.isDigit(lastChar)) {
                    plotMakers[i] = (PlotMaker) Class.forName(base +
                            plotMakerName.substring(0, plotMakerName.length() - 1)).
                            newInstance();
                } else {
                    plotMakers[i] = (PlotMaker) Class.forName(base + plotMakerName).
                    newInstance();
                }

                plotMakers[i].initialize(this);
//				String valueString = plotParametersMap.get(plotMakerName).trim();
//				plotMakers[i].setOptionString(valueString);

                if (Character.isDigit(lastChar)) {
                    plotMakers[i].name += lastChar;
                }
            } catch (Exception ex) {
                System.err.println("Could not create " + plotMakerName);
                ex.printStackTrace(Vision.getInstance().getFileStdOut());

                plotMakers[i] = null;
            }
            i++;
        }
        
        keybindingBuilders = new LinkedList<KeybindingBuilder>();
        for (PlotMaker p : plotMakers)
            if (p instanceof KeybindingBuilder) keybindingBuilders.add((KeybindingBuilder) p);			
    }


    public void setTreeLayout(int mode) {
        this.removeAll();

        if (mode == 0) {
            this.setLayout(new GridLayout(1, 2));
            this.add(leftHalfPanel);
            this.add(rightHalfPanel);
        } else if (mode == 1) {
            this.setLayout(new GridLayout(1, 1));
            this.add(rightHalfPanel);
        } else if (mode == 2) {
            this.setLayout(new GridLayout(1, 1));
            this.add(leftHalfPanel);
        }

        this.validate();
    }


    private JMenu makeFileMenu() {
        JMenu fileMenu = new JMenu("Neuron Viewer: File");

        JMenuItem item = new JMenuItem("Save File");
        item.addActionListener(this);
        item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S, KeyEvent.CTRL_DOWN_MASK));
        fileMenu.add(item);

        fileMenu.add(new JSeparator());

        item = new JMenuItem("Save Classification To File ...");
        item.addActionListener(this);
        fileMenu.add(item);

        item = new JMenuItem("Load Classification From File ...");
        item.addActionListener(this);
        fileMenu.add(item);

        return fileMenu;
    }


    private JMenu makeEditMenu() {
        JMenu fileMenu = new JMenu("Neuron Viewer: Edit");

        JMenuItem item = new JMenuItem("Copy Left Panel To Clipboard");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Clipboard.saveComponentToClipboard(leftSplitPane);
            }
        });
//		item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S, KeyEvent.CTRL_DOWN_MASK));
        fileMenu.add(item);

        item = new JMenuItem("Copy Right Panel To Clipboard");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Clipboard.saveComponentToClipboard(rightSplitPane);
            }
        });
//		item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S, KeyEvent.CTRL_DOWN_MASK));
        fileMenu.add(item);

        return fileMenu;
    }


    private JMenu makeViewMenu() {
        JMenu menu = new JMenu("Neuron Viewer: View");
        {
            ButtonGroup tileGroup = new ButtonGroup();
            JRadioButtonMenuItem item = new JRadioButtonMenuItem("Show Both Trees");
            item.setSelected(true);
            item.setMnemonic(KeyEvent.VK_T);
            item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    setTreeLayout(0);
                }
            });
            tileGroup.add(item);
            menu.add(item);

            item = new JRadioButtonMenuItem("Show Only Right Tree");
            item.setMnemonic(KeyEvent.VK_H);
            item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    setTreeLayout(1);
                }
            });
            tileGroup.add(item);
            menu.add(item);

            item = new JRadioButtonMenuItem("Show Only Left Tree");
            item.setMnemonic(KeyEvent.VK_V);
            item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    setTreeLayout(2);
                }
            });
            tileGroup.add(item);
            menu.add(item);
        }

        menu.add(new JSeparator());

        JMenuItem item;

        for (final String name : parametersMap.keySet()) {
            if (!name.startsWith("show")) {
                continue;
            }

            boolean s = Boolean.valueOf(parametersMap.get(name)).booleanValue();
            try {
                NeuronViewer.class.getField(name).set(NeuronViewer.this, s);
            } catch (Exception ex) {
                ex.printStackTrace();
            }

            item = new JCheckBoxMenuItem(name, s);
            item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    String v = parametersMap.get(name);
                    if (v.equals("true")) {
                        v = "false";
                    } else {
                        v = "true";
                    }
                    parametersMap.put(name, v);

                    try {
                        NeuronViewer.class.getField(name).set(
                                NeuronViewer.this, ( (JMenuItem) e.getSource()).isSelected());
                    } catch (Exception ex) {
                        ex.printStackTrace();
                    }

                    updateTrees();
                }
            });
            menu.add(item);
        }

        return menu;
    }


    public static boolean isValidFolder(TreePath classTreePath) {
        Object o = ( (DefaultMutableTreeNode) classTreePath.getLastPathComponent()).getUserObject();
        if (o instanceof String) {
            return true;
        } else {
            return false;
        }
    }


    public void updateActions() {
        ParameterModel model = (ParameterModel) actionsParamsTable.getTreeTableModel();
        ParameterModel.TreeTableNode root = (ParameterModel.TreeTableNode)
        model.getRoot();
        for (int i = 0; i < root.getChildCount(); i++) {
            ParameterModel.TreeTableNode plotNode = root.getChild(i);
            CalculationAction action = getActionOfType(plotNode.getName());
            if (action == null) {
                continue;
            }

            for (int j = 0; j < plotNode.getChildCount(); j++) {
                ParameterModel.TreeTableNode paramNode = plotNode.getChild(j);
                try {
                    Field f = action.getClass().getField(paramNode.getName());
                    f.set(action, paramNode.p.valueAsObject());
                } catch (Exception ex) {
                    ex.printStackTrace();
                }
            }
        }
    }


    public void updateTrees() {
        ParameterModel model = (ParameterModel) plotsParamsTable.getTreeTableModel();
        ParameterModel.TreeTableNode root = (ParameterModel.TreeTableNode)
        model.getRoot();
        for (int i = 0; i < root.getChildCount(); i++) {
            ParameterModel.TreeTableNode plotNode = root.getChild(i);
            PlotMaker plotMaker = getPlotMakerOfType(plotNode.getName());
            if (plotMaker == null) {
                continue;
            }

            try {
                Field f = plotMaker.getClass().getField("isSelected");
                f.set(plotMaker, plotNode.p.valueAsObject());
            } catch (Exception ex) {
                ex.printStackTrace();
            }

            for (int j = 0; j < plotNode.getChildCount(); j++) {
                ParameterModel.TreeTableNode paramNode = plotNode.getChild(j);
                try {
                    Field f = plotMaker.getClass().getField(paramNode.shortName);
                    if (f.getType().isEnum()) {
                        int v = (int) Double.parseDouble(paramNode.p.valueAsString());
                        f.set(plotMaker, f.getType().getEnumConstants()[v]);
                    } else {
                        f.set(plotMaker, paramNode.p.valueAsObject());
                    }
                } catch (Exception ex) {
                    ex.printStackTrace();
                }
            }
        }

        TreePath tp = rightTree.getSelectionPath();
        rightTree.setSelectionPath(null);
        rightTree.setSelectionPath(tp);

        tp = leftTree.getSelectionPath();
        leftTree.setSelectionPath(null);
        leftTree.setSelectionPath(tp);
    }


    public void actionPerformed(ActionEvent event) {
        String action = event.getActionCommand();

        if (action.equals("Save File")) {
            new Thread() {
                public void run() {
                    NeuronViewer.this.leftTree.setEnabled(false);
                    NeuronViewer.this.rightTree.setEnabled(false);
                    try {
                        paramsFile.flush();
                    } catch (IOException e) {
                        Vision.reportException(e);
                    }
                    NeuronViewer.this.leftTree.setEnabled(true);
                    NeuronViewer.this.rightTree.setEnabled(true);
                }
            }.start();
        } else if (action.equals("Save Classification To File ...")) {
        /*	JFileChooser fc = new JFileChooser(new File(filePathRoot).getParent());
            fc.setSelectedFile(new File(datasetName + ".classification.txt"));
            if (fc.showSaveDialog(mainFrame) == JFileChooser.APPROVE_OPTION) {
                try {
                    PrintWriter w = new PrintWriter(new FileWriter(fc.getSelectedFile()));
                    LinkedHashMap idMap = paramsFile.getClassIDs();
                    Object[] idList = idMap.keySet().toArray();
                //	Arrays.sort(idList);
                    for (int i = 0; i < idList.length; i++) {
                        Integer key = (Integer) idList[i];
                        w.println(key.intValue() + "  " + idMap.get(key));
                    }
                    w.close();
                } catch (IOException e) {
                    Vision.reportException(e);
                }		
            }*/
            
            JFileChooser fc = new JFileChooser(new File(filePathRoot).getParent());
            fc.setSelectedFile(new File(datasetName + ".classification.txt"));
            if (fc.showSaveDialog(mainFrame) == JFileChooser.APPROVE_OPTION) {
                try {
                    PrintWriter w = new PrintWriter(new FileWriter(fc.getSelectedFile()));
                    ParameterModel model = (ParameterModel) actionsParamsTable.getTreeTableModel();
                    DefaultMutableTreeNode root = (DefaultMutableTreeNode) rightTree.getModel().getRoot();
                    writeFolderNeurons(root, "", w);
                    w.close();
                }  catch (IOException e) {
                    System.exit(2);
                    Vision.reportException(e);
                }	
            }
            
        } else if (action.equals("Load Classification From File ...")) {
            JFileChooser fc = new JFileChooser(new File(filePathRoot).getParent());
            fc.setSelectedFile(new File(datasetName + ".classification.txt"));
            if (fc.showOpenDialog(mainFrame) == JFileChooser.APPROVE_OPTION) {
                loadClassification(fc.getSelectedFile().getAbsolutePath());
            }
        }
    }
    
    public void writeFolderNeurons(DefaultMutableTreeNode folder, String folderName, PrintWriter w) throws IOException {

           for (int i = 0; i < folder.getChildCount(); i++) {
                if (folder.getAllowsChildren()) {
                    writeFolderNeurons( (DefaultMutableTreeNode) folder.getChildAt(i), folderName + folder + "/", w);
                } 
            }
           
            if(!folder.getAllowsChildren()) {
                
                w.println(folder + "  " + folderName);
            }
    }


    /**
     * Gets called when one of the trees is modified.
     */
    public void valueChanged(TreeSelectionEvent event) {
        DefaultMutableTreeNode o = (DefaultMutableTreeNode) event.getPath().
        getLastPathComponent();
        Object selectedObject = o.getUserObject();
        int leftDivider = leftSplitPane.getDividerLocation();
        int rightDivider = rightSplitPane.getDividerLocation();
                
        // Make sure we don't get stale keybindings
        for (KeybindingBuilder k : keybindingBuilders) k.clearKeybindings();		
        
        // Build plots
        Component panelA = null, panelB = null;		
        if (selectedObject instanceof String) { // a folder was selected
            try {
                IntegerList list = InteractiveTree.getNeuronsInClass(o,
                        new IntegerList(), showNeuronsInSubclasses);

                if (plotMakers[0] != null) {
                    panelA = (Component) plotMakers[0].makePlot(
                            list, PlotMaker.CLASS_PLOT, event.getPath());
                } else {
                    panelA = new JLabel("No mosaic Available");
                }

                panelB = makeClassPanel(list, event.getPath());
            } catch (IOException e) {
                Vision.reportException(e);
            }
        } else { // a neuron was selected
            try {
                IntegerList list = new IntegerList();
                list.add( (Integer) selectedObject);

                if (plotMakers[1] != null) {
                    panelA = (Component) plotMakers[1].makePlot(
                            list, PlotMaker.NEURON_PLOT, event.getPath());
                } else {
                    panelA = new JLabel("No STA Available");
                }

                panelB = makeNeuronPanel(list, event.getPath());
            } catch (IOException ex) {
                ex.printStackTrace();
                return;
            }
        }

        // Collect any keybindings (AKA hotkeys) created for the plots
        List<PlotMaker.Keybinding> keybindings = new LinkedList<PlotMaker.Keybinding>();
        for (KeybindingBuilder k : keybindingBuilders)
            keybindings.addAll(k.getKeybindings());
        
        // Determine whether to operate on left or right panel depending on where event originated
        InteractiveTree sourceTree = (InteractiveTree) event.getSource();
        JSplitPane splitPane;
        JPanel halfPanel;
        if (sourceTree == leftTree) {
            splitPane = leftSplitPane;
            halfPanel = leftHalfPanel;
            leftSelection = selectedObject;
            splitPane.setDividerLocation(leftDivider);
        } else if (sourceTree == rightTree) {
            splitPane = rightSplitPane;
            halfPanel = rightHalfPanel;
            rightSelection = selectedObject;
            splitPane.setDividerLocation(rightDivider);
        } else {
            throw new Error("unknown sourceTree");
        }
        
        // Set up the appropriate panel
        splitPane.setLeftComponent(panelA);
        splitPane.setRightComponent(panelB);	
        for (PlotMaker.Keybinding k : keybindings) {
            halfPanel.getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(k.key, k.handle);
            halfPanel.getActionMap().put(k.handle, k.action);
        }
        halfPanel.validate();
        
        if (correlator != null && leftSelection instanceof Integer && rightSelection instanceof Integer)
            correlator.set( (Integer) leftSelection, (Integer) rightSelection);
    }


    /**
     * Called by the InteractiveTree whenever a change to the tree happens
     * @param e
     */
    public void stateChanged(ChangeEvent e) {
        HashMap<Integer, String> state = rightTree.getState();
        for (Integer key : state.keySet()) {
            paramsFile.setCell(key.intValue(), "classID", state.get(key));
        }
    }


    public void classificationHappened(int[] idList, TreePath classPath) {
        rightTree.setClassLabels(idList, classPath);
    }


    public void exit(boolean saveParams) {
        // remove the keystroke listeners
//		System.err.println(keyStrokes.size() + " Viewer keystrokes");

        for (KeyStroke s : keyStrokes) {
            eventQueue.removeKeyStroke(s);
        }

        // save NeuronViewer state
        try {
            // save the first group
            parametersMap.put("leftDividerSize", "" + leftSplitPane.getDividerLocation());
            parametersMap.put("rightDividerSize", "" + rightSplitPane.getDividerLocation());

            configuration.setParameterGroup("Neuron Viewer 1", parametersMap);
            configuration.setParameterGroup("Neuron Viewer Plots", plotsParamsTable);
        } catch (Exception ex1) {
            ex1.printStackTrace();
        }

        // save parameters if needed
        try {
            paramsFile.close(saveParams);
        } catch (IOException ex) {
            ex.printStackTrace();
        }

        // close all other files
        if (neuronFile != null) {
            try {
                neuronFile.close();
            } catch (IOException ex2) {
                System.err.println("Could not close neuron file");
            }
        }

        if (staCollection != null) {
            try {
                staCollection.close();
            } catch (IOException ex2) {
                System.err.println("Could not close sta file");
            }
        }

        if (stvCollection != null) {
            try {
                stvCollection.close();
            } catch (IOException ex2) {
                System.err.println("Could not close stv file");
            }
        }

        if (imgFile != null) {
            try {
                imgFile.close();
            } catch (IOException ex2) {
                System.err.println("Could not close ei file");
            }
        }

        // destroy the viewer
        if (internal) {
            FocusManager.getCurrentManager().clearGlobalFocusOwner();
        } else {
            ((JFrame) mainFrame).dispose();
        }
    }

    /**
     *  Closes the neuron viewer.  Should be used when the neuron viewer
     *  was opened outside of the vision gui.
     * 
     */

    public void externalClose() {
        try {
            if (paramsFile.needsFlush()) {
                int response = JOptionPane.showOptionDialog(
                        mainFrame,
                        "The Parameters File was Changed. Do you want to save it?",
                        "Warning",
                        JOptionPane.YES_NO_CANCEL_OPTION,
                        JOptionPane.WARNING_MESSAGE, null, null, null);
                if (response == JOptionPane.YES_OPTION) {
                    paramsFile.close(true);
                } else if (response == JOptionPane.NO_OPTION) {
                    paramsFile.close(false);
                }
            } else {
                paramsFile.close(false);
            }
            if (staCollection != null) {
                staCollection.close();
            }
            if (stvCollection != null) {
                stvCollection.close();
            }
            if (neuronFile != null) {
                neuronFile.close();
            }
            if (imgFile != null) {
                imgFile.close();
            }

        } catch (IOException ex) {
            ex.printStackTrace();
        }	
        ((JFrame) mainFrame).dispose();

    }


    /**
     * Returns true if the application can exit and false if it cannot exit because
     * there is unsaved data
     */
    public boolean canClose() {
        if (paramsFile.needsFlush()) {
            int response = JOptionPane.showOptionDialog(
                    mainFrame,
                    "The parameters file " + experimentName +
                    File.separator + datasetName + " was changed. Do you want to save it?",
                    "Warning",
                    JOptionPane.YES_NO_CANCEL_OPTION,
                    JOptionPane.WARNING_MESSAGE, null, null, null);
            if (response == JOptionPane.YES_OPTION) {
                exit(true);
                return true;
            } else if (response == JOptionPane.NO_OPTION) {
                exit(false);
                return true;
            }
            return false;
        } else {
            exit(false);
            return true;
        }
    }


    private JPanel makeNeuronPanel(IntegerList list, TreePath path) throws IOException {
        JPanel panel = new JPanel(new GridLayout(0, 1));

        // add the automatic plots
        for (int i = 2; i < plotMakers.length; i++) {
            if (plotMakers[i] != null && plotMakers[i].isSelected) {
                if (plotMakers[i].getPlotType() == PlotMaker.NEURON_PLOT ||
                        plotMakers[i].getPlotType() == PlotMaker.CLASS_AND_NEURON_PLOT) {
                    Object plots = null;
                    try {
                        plots = plotMakers[i].makePlot(list, PlotMaker.NEURON_PLOT, path);
                    } catch (Exception ex) {
                        System.err.println(
                                "Plot " + plotMakers[i].name + " failed to execute.\n" +
                        "See vision-output.txt for error message.");
                        ex.printStackTrace(Vision.getInstance().getFileStdOut());
                    }

                    if (plots instanceof Component) {
                        panel.add( (Component) plots);
                    } else if (plots instanceof Component[]) {
                        for (int j = 0; j < ( (Component[]) plots).length; j++) {
                            panel.add( ( (Component[]) plots)[j]);
                        }
                    }
                }
            }
        }

        return panel;
    }


    /**
     * Creates the panel containing all class information
     * (called when a folder is selected)
     */
    private JPanel makeClassPanel(IntegerList list, TreePath path) throws IOException {
        JPanel panel = new JPanel(new GridLayout(0, 1));

        // add the automatic plots
        for (int i = 2; i < plotMakers.length; i++) {
            if (plotMakers[i] != null && plotMakers[i].isSelected) {
                if (plotMakers[i].getPlotType() == PlotMaker.CLASS_PLOT ||
                        plotMakers[i].getPlotType() == PlotMaker.CLASS_AND_NEURON_PLOT) {

                    Object plots = null;
                    try {
                        plots = plotMakers[i].makePlot(list, PlotMaker.CLASS_PLOT, path);
                    } catch (Exception ex) {
                        System.err.println(
                                "Plot " + plotMakers[i].name + " failed to execute.\n" +
                        "See vision-output.txt for error message.");
                        ex.printStackTrace(Vision.getInstance().getFileStdOut());
                    }

                    if (plots instanceof Component) {
                        panel.add( (Component) plots);
                    } else if (plots instanceof Component[]) {
                        for (int j = 0;
                        j < ( (Component[]) plots).length; j++) {
                            panel.add( ( (Component[]) plots)[j]);
                        }
                    }
                }
            }
        }

        return panel;
    }


    public DefaultMutableTreeNode select(String selection) {
        boolean found = false;

        // select the class in the right tree
        String[] s = StringUtil.decomposeString(selection, "/");
        DefaultMutableTreeNode folder =
            (DefaultMutableTreeNode) rightTree.getModel().getRoot();
        for (int j = 0; j < s.length; j++) {
            found = false;
            for (int i = 0; i < folder.getChildCount(); i++) {
                DefaultMutableTreeNode c = (DefaultMutableTreeNode) folder.getChildAt(i);
                if (c.getUserObject().toString().equals(s[j])) {
                    folder = c;
                    found = true;
                    break;
                }
            }
        }
        rightTree.setSelectionPath(new TreePath(folder.getPath()));

        if (found == true) {
            return folder;
        } else {
            return null;
        }
    }


    public static Component getChild(Component c, final String name) {
        String n = c.getName();
//		System.err.println(n + ", " + c);
        if (n != null && n.equals(name)) {
            return c;
        } else if (c instanceof JSplitPane) {
            Component c1 = getChild( ( (JSplitPane) c).getLeftComponent(), name);
            if (c1 != null) {
                return c1;
            }
            return getChild( ( (JSplitPane) c).getRightComponent(), name);
        } else if (c instanceof Container) {
            for (int i = 0; i < ( (Container) c).getComponentCount(); i++) {
                Component c1 = getChild( ( (Container) c).getComponent(i), name);
                if (c1 != null) {
                    return c1;
                }
            }
        }
        return null;
    }


    public Component getPlot(String name) {
        return (Component) getChild(rightSplitPane, name);
    }


    public void loadClassification(String fileName) {
        LinkedHashMap<Integer, String> classLabels = new LinkedHashMap<Integer, String>();
        LineNumberReader r = null;
        String line = "";
        /*---- HAAAAAACK
         * This hack is to allow loading classifications that contain extra neuron IDs without crashing vision
         * As may happen during Spectra/mVision concatenated runs
         * Because empty neurons in 1 sub dataset are discarded, but may be kept in other sub-datasets.
         * @author Vincent Deo - Stanford University - 11/17/2015
         */
        int[] idList = paramsFile.getIDList();
        Arrays.sort(idList);
        //-----
        try {
            r = new LineNumberReader(new FileReader(fileName));
            
            while ( (line = r.readLine()) != null) {
                line = line.trim();
                int index = line.indexOf(' ');
                int id = Integer.parseInt(line.substring(0, index));
                String classLabel = line.substring(index, line.length()).trim();
                //---- HAAAAAACK
                if (Arrays.binarySearch(idList, id) >= 0) {
                //-----
                    classLabels.put(new Integer(id), classLabel);
                //
                }
                //
            }
        } catch (IOException e) {
            System.out.println(line);
            if (r != null) {
                try {
                    r.close();
                } catch (IOException e2) {}
            }
            Vision.reportException(e);
        }

        DefaultTreeModel treeModel = InteractiveTree.makeModel(classLabels);
        leftTree.setModel(treeModel);
        rightTree.setModel(treeModel);
    }


    public static Point[] readCorners(File f, int nPoints) throws IOException {
        Point[] p = new Point[nPoints];
        LineNumberReader fr = new LineNumberReader(new FileReader(f));
        for (int i = 0; i < nPoints - 1; i++) {
            String s = fr.readLine();
            int ind = s.indexOf(" ");
            if (ind < 0) {
                throw new IOException("No space in line " + fr.getLineNumber());
            }
            int x = -1, y = -1;
            try {
                x = Integer.parseInt(s.substring(0, ind).trim());
            } catch (NumberFormatException e) {
                throw new IOException("Not an integer: " + s.substring(0, ind));
            }
            try {
                y = Integer.parseInt(s.substring(ind, s.length()).trim());
            } catch (NumberFormatException e) {
                throw new IOException("Not an integer: " + s.substring(ind, s.length()));
            }
            p[i] = new Point(x, y);
        }
        p[nPoints - 1] = p[0];
        fr.close();
        return p;
    }


    public static PlotPanel makeHistogram(
            final NeuronViewer nv, String expression, double binWidth,
            final TreePath classTreePath) throws CannotEvaluateException {

        String classPath = InteractiveTree.pathToString(classTreePath);
        int[] id = nv.paramsFile.getNeuronsInClass(classPath);

        return makeHistogram(nv, expression, binWidth, Double.NaN, Double.NaN,
                classTreePath, id);
    }


    public static PlotPanel makeHistogram(
            final NeuronViewer nv, String expression, double binWidth, double x1,
            double x2, final TreePath classTreePath, int[] idList) throws
            CannotEvaluateException {

        final HashMap<Integer,Double> valueMap = nv.paramsFile.evaluate(expression, idList);
        return makeHistogram(nv, valueMap, binWidth, x1, x2, classTreePath, "");
    }


    public static PlotPanel makeHistogram(
            final NeuronViewer nv, final HashMap<Integer, Double> valueMap, double binWidth,
            double x1, double x2, final TreePath classTreePath, String valueName) throws
            CannotEvaluateException {

        if (binWidth == 0) {
            throw new CannotEvaluateException("binWidth cannot be zero");
        }

        if (x2 <= x1) {
            throw new CannotEvaluateException("x2 must be greater than x1");
        }


        if (valueMap == null) {
            throw new CannotEvaluateException("valueMap == null");
        }

        if (Double.isNaN(x1) || Double.isNaN(x2)) {
            double xMin = Double.POSITIVE_INFINITY;
            double xMax = Double.NEGATIVE_INFINITY;
            for (Integer id : valueMap.keySet()) {
                Double ox = valueMap.get(id);
                if (ox != null) {
                    double x = ox.doubleValue();
                    if (x > xMax) {
                        xMax = x;
                    }
                    if (x < xMin) {
                        xMin = x;
                    }
                }
            }

            if (Double.isInfinite(xMin) && Double.isInfinite(xMax)) {
                return null;
            }

            if (xMin != xMax) {
                x1 = xMin + ( -0.1) * (xMax - xMin);
                x2 = xMax + ( +0.1) * (xMax - xMin);
            } else {
                x1 = xMin + ( -0.1) * xMin;
                x2 = xMax + ( +0.1) * xMin;
            }
        }

        DoubleHistogram h = new DoubleHistogram("", x1, x2, binWidth);
        MeanVarianceCalculator mvc = new MeanVarianceCalculator();
        for (Integer id : valueMap.keySet()) {
            Double v = valueMap.get(id);
            if (v != null) {
                h.fill(v, 1);
                mvc.add(v);
            }
        }

        PlotPanel p = new PlotPanel("histogram");
        p.addToLegend(valueName + " = " +
                new Num(mvc.getMean(), mvc.getMeanVariance()).toString(3));
        if (nv != null) {
            p.addSelectionAction(new SelectionAction("New Class") {
                public void selectionPerformed(JComponent source, Selection selection) {
                    Rectangle2D b = selection.getSelection().getBounds2D();
                    double x1 = b.getX();
                    double x2 = x1 + b.getWidth();
                    IntegerList idList = new IntegerList();

                    for (Integer id : valueMap.keySet()) {
                        Double ox = valueMap.get(id);
                        if (ox != null) {
                            double x = ox.doubleValue();
                            if (x > x1 && x < x2) {
                                idList.add(id.intValue());
                            }
                        }
                    }

                    nv.classificationHappened(idList.toArray(), 
                            classTreePath.pathByAddingChild(nv.getNewClass()));
                }
            });
        }
        HistogramStyle hs = new HistogramStyle("style");
        hs.setOutlineThickness(0.5f);
        p.addData(h, hs);
        p.autoscale();

        return p;
    }

    public PlotPanel makeScatterPlot(String xExpression, String yExpression, String xLabel, String yLabel, TreePath classTreePath) throws
            CannotEvaluateException {
        return makeScatterPlot(xExpression, yExpression, xLabel, yLabel, classTreePath, 
                paramsFile.getNeuronsInClass(InteractiveTree.pathToString(classTreePath)));
    }

    public PlotPanel makeScatterPlot(
            String xExpression, String yExpression, TreePath classTreePath, int[] idList) throws
            CannotEvaluateException {

        return makeScatterPlot(
                xExpression, yExpression, xExpression, yExpression, classTreePath, idList);
    }


    public PlotPanel makeScatterPlot(
            String xExpression, String yExpression, String xLabel, String yLabel,
            TreePath classTreePath, int[] idList) throws CannotEvaluateException {

//		String cutExpression =
//			"classID==\"" + InteractiveTree.pathToString(classTreePath) + "\"";
//		HashMap<Integer, Double> xMap = paramsFile.evaluateEx(xExpression, cutExpression);
//		HashMap<Integer, Double> yMap = paramsFile.evaluateEx(yExpression, cutExpression);
        
        
        HashMap<Integer, Double> xMap = paramsFile.evaluate(xExpression, idList);
        HashMap<Integer, Double> yMap = paramsFile.evaluate(yExpression, idList);
        
        
        if (xMap == null || yMap == null) {
            return new PlotPanel();
        }

        return makeScatterPlot(xMap, yMap, classTreePath, xLabel, yLabel);
    }


    public PlotPanel makeScatterPlot(final HashMap<Integer, Double> xMap,
            final HashMap<Integer, Double> yMap,
            final TreePath classTreePath,
            String xLabel, String yLabel) {

        ScatterPlot sp = new ScatterPlot();
        for (Integer id : xMap.keySet()) {
            Double x0 = xMap.get(id);
            Double y0 = yMap.get(id);
            if (x0 != null && y0 != null) {
                double x = x0.doubleValue();
                double y = y0.doubleValue();
                if (!Double.isNaN(x) && !Double.isNaN(y)) {
                    sp.add(x, y);
                }
            }
        }
        PlotPanel p = new PlotPanel();
        p.setAntiallias(true);
        p.addSelectionAction(new SelectionAction("New Class") {
            public void selectionPerformed(JComponent source, Selection selection) {
                SelectionRegion r = selection.getSelection();
                IntegerList idList = new IntegerList();
                for (Integer id : xMap.keySet()) {
                    Double ox = xMap.get(id);
                    Double oy = yMap.get(id);
                    double x = (ox != null) ? ox.doubleValue() : Double.NaN;
                    double y = (oy != null) ? oy.doubleValue() : Double.NaN;
                    if (r.contains(x, y)) {
                        idList.add(id.intValue());
                    }
                }

                DefaultMutableTreeNode newClassNode = getNewClass();
                TreePath newClassPath = classTreePath.pathByAddingChild(newClassNode);
                classificationHappened(idList.toArray(), newClassPath);
            }
        });
        p.setLabels(xLabel, yLabel);
        p.addData(sp, scatterPlotStyle);
        p.autoscale();
        p.pad();

        return p;
    }

    DefaultMutableTreeNode getNewClass() { return getNewClass("nc"); }
    DefaultMutableTreeNode getNewClass(String prefix) {
        return new DefaultMutableTreeNode(prefix + newClassIndex++, true);
    }

    public PlotMaker getPlotMaker(String name) {
        for (int i = 0; i < plotMakers.length; i++) {
            if (plotMakers[i].getName().equals(name)) {
                return plotMakers[i];
            }
        }
        return null;
    }


    public PlotMaker getPlotMakerOfType(String type) {
        for (int i = 0; i < plotMakers.length; i++) {
            if (plotMakers[i] != null && plotMakers[i].getClass().getName().endsWith(type)) {
                return plotMakers[i];
            }
        }
        return null;
    }


    public CalculationAction getActionOfType(String type) {
        for (int i = 0; i < actions.size(); i++) {
            CalculationAction a = actions.get(i);
            if (a != null && a.getClass().getName().endsWith(type)) {
                return a;
            }
        }
        return null;
    }


    public ParametersFile getParametersFile() {
        return paramsFile;
    }
}