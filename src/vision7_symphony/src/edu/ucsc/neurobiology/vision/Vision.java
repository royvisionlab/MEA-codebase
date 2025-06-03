package edu.ucsc.neurobiology.vision;

import java.beans.*;
import java.io.*;
import java.net.*;
import java.util.*;
import java.util.jar.*;

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.plaf.*;

import org.jdesktop.swingx.*;
import edu.ucsc.neurobiology.vision.actions.*;
import edu.ucsc.neurobiology.vision.calculations.*;
import edu.ucsc.neurobiology.vision.gui.*;
import edu.ucsc.neurobiology.vision.gui.Desktop;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.tasks.*;
import edu.ucsc.neurobiology.vision.util.*;

import com.martiansoftware.jsap.*;


/**
 * The main VISION application class.<p></p>
 * This class contains a set of functionality which is vital for the whole application,
 * which includes:
 * <ul>
 *   <li> methods for creating, opening and manipulating VISION windows
 *   <li> methods for accessing the <tt>Config</tt> and <tt>DataManager</tt> objects.
 *   <li> methods for reporting exceptions.
 *   <li> methods for accessing the VISION status bar
 * </ul>
 * <p></p>
 * This is a SINGLETON class. Classes which want access to the methods of this
 * class must use the <tt>getInstance()</tt> method to get the real VISION instance.
 *
 * @author Charles A. Loomis, University of California, Santa Cruz<br>
 *         Dumitru Petrusca, University of California, Santa Cruz
 */
public class Vision
implements ActionListener, Desktop.SelectionListener {
    
    /**
     * This keeps track of the only instance of this class.
     */
    private static Vision singleton = null;
    private static Config config;
    private static StartupWindow startupWindow;

    private Map<JInternalFrame, JComponent> controllers = new HashMap<JInternalFrame, JComponent>();
    private Map<JInternalFrame, JMenu[]>    menus       = new HashMap<JInternalFrame, JMenu[]>();
    private ViewManager viewManager;
    private DataManager dataManager;
    private CalculationManager calculationsManager;
    private CalculationManagerGUI calculationsManagerGUI;
    private Desktop mainDesktop;
    private JMenuBar menuBar;
    private JMenu helpMenu;
    private JMenu fileMenu;
    private JInternalFrame currentFrame;
    private JFileChooser fileDialog;
    private JFrame mainFrame;
    private StatusBar statusBar;
    private static boolean guiEnabled = false;
    private HashMap<String, String> visionParametersGroup;
    private JXTaskPaneContainer controlPane;
    private JScrollPane terminalComponent;
    private JSplitPane splitPane;
    private int mainDividerSize = 300;
    private JXTaskPane userHelpPane, developerHelpPane;
    private PrintStream fileStandardOutput;


    /**
     * This constructor should not be called by user code.
     * Instead the factory method getInstance() should be used.
     * This will throw an IllegalStateException if the singleton class already exists.
     */
    private Vision() throws Exception {
        // Check to make sure that this hasn't already been called.
        if (singleton != null) {
            throw new IllegalStateException("Singleton already exists.");
        }

        // Really before doing anything else, mark this as being the singleton.
        singleton = this;
        calculationsManager = new CalculationManager();
    }

    //Must be run after getInstance, and before the config.xml file is accessed.
    public void setConfig(String configString) throws Exception {
        try {
            config = new Config(configString);
        } catch (Exception e) {
            System.out.println("Problem loading configuration file: " + configString);
            e.printStackTrace();
            throw new IllegalArgumentException("Could not read the config file.");
        }
    }


    public void updateTitle() {
        if (currentFrame == null) {
            mainFrame.setTitle("Vision");
        } else {
            mainFrame.setTitle("Vision - " + currentFrame.getTitle());
        }
    }


    public PrintStream getFileStdOut() {
        if (fileStandardOutput != null) {
            return fileStandardOutput;
        } else {
            return System.err;
        }
    }


    /**
     * Returns the <tt>JFrame</tt> which contains the whole application.
     */
    public JFrame getMainFrame() {
        return mainFrame;
    }


    /**
     * This method will display the given message in the left side of the application
     * status bar.
     */
    public void sendMessage(String message) {
        if (guiEnabled) {
            statusBar.setMessage(message);
        } else {
//			System.out.println("Calculation Status: " + message);
        }
    }


    public void startProgressBar() {
        if (!guiEnabled) {
            String s = "================= Progress Bar ===================";
            System.out.println(s);
        }
    }


    public void endProgressBar() {
        if (!guiEnabled) {
            System.out.println();
        }
    }


    /**
     * This method will change the value of the progress monitor in the right side of
     * the application status bar.
     *
     * @param progress the percentage progress value
     */
    int nPrintedDots = 0;
    public void setProgress(int progress) {
        if (guiEnabled) {
            statusBar.setProgress(progress);
        } else {
            int nDotsToPrint = progress / 2 - nPrintedDots;
            for (int i = 0; i < nDotsToPrint; i++) {
                System.out.print("*");
            }
            nPrintedDots += nDotsToPrint;
        }
    }


    /**
     * The Vision class should be created with this factory method. This
     * is a singleton class, so only one instance of this class will
     * be created. If this is called more than once, a reference to
     * the original instance will be returned.
     */
    public static Vision getInstance() {
        if (singleton == null) {
            try {
                singleton = new Vision();
            } catch (Exception e) {
                reportFatalException("Vision could not start because:", e);
            }
        }
        return singleton;
    }


    /**
     * Returns the <tt>Config</tt> object containing the information read from the
     * configuration file.
     */
    public static Config getConfig() {
        return config;
    }


    /**
     * Returns a String containing the VISION version.
     */
    public String getVersion() {
        return "Development Version for v1.0";
    }


    /**
     * Returns the application-wide data manager. The <tt>DataManager</tt> class
     * contains all the open data in VISION.
     */
    public DataManager getDataManager() {
        return dataManager;
    }


    public CalculationManager getCalculationManager() {
        return calculationsManager;
    }


    public CalculationManagerGUI getCalculationManagerGUI() {
        if (isConsoleBased()) {
            throw new IllegalStateException("Vision is in Console mode");
        }
        return calculationsManagerGUI;
    }


    public ViewManager getViewManager() {
        return viewManager;
    }


    private JMenu createViewMenu() {
        JMenu menu = new JMenu("View");
        menu.setMnemonic(KeyEvent.VK_V);

        /*
                item = new JCheckBoxMenuItem("View Controler", true);
                item.setMnemonic(KeyEvent.VK_C);
                item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_C,
                    ActionEvent.ALT_MASK + ActionEvent.SHIFT_MASK));
                item.addItemListener(new ItemListener() {
                    int n;
                    public void itemStateChanged(ItemEvent e) {
                        switch (e.getStateChange()) {
                            case ItemEvent.DESELECTED:
                                n = splitPane.getDividerLocation();
                                splitPane.setDividerLocation(0);
                                break;
                            case ItemEvent.SELECTED:
                                splitPane.setDividerLocation(n);
                                break;
                        }
                        mainFrame.validate();
                        mainDesktop.updateLayout();
                    }
                });
                menu.add(item);
         */
        {
            ButtonGroup tileGroup = new ButtonGroup();
            JRadioButtonMenuItem item = new JRadioButtonMenuItem("Show Output at Left");
            item.setSelected(true);
            item.setMnemonic(KeyEvent.VK_T);
            item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    setTerminal(1);
                }
            });
            tileGroup.add(item);
            menu.add(item);

            item = new JRadioButtonMenuItem("Show Output at Right");
            item.setMnemonic(KeyEvent.VK_H);
            item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    setTerminal(2);
                }
            });
            tileGroup.add(item);
            menu.add(item);

            item = new JRadioButtonMenuItem("Do Not Show Output");
            item.setMnemonic(KeyEvent.VK_V);
            item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    setTerminal(0);
                }
            });
            tileGroup.add(item);
            menu.add(item);
        }

        menu.add(new JSeparator());

        JCheckBoxMenuItem item = new JCheckBoxMenuItem("View Status Bar", true);
        item.setMnemonic(KeyEvent.VK_B);
        item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S,
                ActionEvent.CTRL_MASK + ActionEvent.SHIFT_MASK));
        item.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                switch (e.getStateChange()) {
                case ItemEvent.DESELECTED:
                    mainFrame.getContentPane().remove(statusBar);
                    break;
                case ItemEvent.SELECTED:
                    mainFrame.add(statusBar, BorderLayout.SOUTH);
                    break;
                }

                mainFrame.validate();
                mainDesktop.updateLayout();
            }
        });
        menu.add(item);

        return menu;
    }


    private JMenu createFileMenu() {
        JMenu menu = new JMenu("File");
        menu.setMnemonic(KeyEvent.VK_F);

        // Read from a file.
        JMenuItem item = new JMenuItem("Open Location", VisionImage.getIcon("open.gif"));
        item.setMnemonic(KeyEvent.VK_O);
        item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_O, ActionEvent.CTRL_MASK));
        item.addActionListener(this);
        menu.add(item);

        // refresh the folder
        item = new JMenuItem("Refresh Open Folder");
//		item.setMnemonic(KeyEvent.VK_O);
//		item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_O, ActionEvent.CTRL_MASK));
        item.addActionListener(this);
        menu.add(item);

        menu.add(new JSeparator());

        // Export options.
        item = new JMenuItem("Export View...");
        item.setMnemonic(KeyEvent.VK_V);
        item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_V, ActionEvent.CTRL_MASK));
        item.setActionCommand("ExportView");
        item.addActionListener(this);
        menu.add(item);

        item = new JMenuItem("Export Desktop...");
        item.setMnemonic(KeyEvent.VK_D);
        item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_D, ActionEvent.CTRL_MASK));
        item.setActionCommand("ExportDesktop");
        item.addActionListener(this);
        menu.add(item);

        item = new JMenuItem("Export ScreenShot...");
        item.setActionCommand("ExportScreenShot");
        item.addActionListener(this);
        menu.add(item);

        menu.add(new JSeparator());

        item = new JMenuItem("Exit");
        item.setMnemonic(KeyEvent.VK_X);
        item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_X, ActionEvent.CTRL_MASK));
        item.setActionCommand("Exit");
        item.addActionListener(this);
        menu.add(item);

        return menu;
    }


    private static javax.help.HelpSet openHelpSet(String helpset) {
        URL url = null;
        try {
            String name = System.getProperty("user.dir") + File.separator + helpset;
            if (new File(name).exists()) {
                url = new URL("file:///" + name);
            }
        } catch (Exception ex) {}

        if (url == null) {
            try {
                url = javax.help.HelpSet.findHelpSet(null, helpset);
            } catch (Exception ex) {}
        }

        if (url == null) {
            //never fully implemented.
//			System.err.println("Could not find " + helpset + " location");
            return null;
        }

        try {
            return new javax.help.HelpSet(null, url);
        } catch (Exception ex1) {
            System.err.println("Could not open " + helpset);
            return null;
        }
    }


    /** Make the menu. */
    private JMenuBar makeMenu() {
        // Make a menubar for the toolbar.
        JMenuBar bar = new JMenuBar();
        bar.setBorderPainted(false);

        // the file menu
        fileMenu = createFileMenu();
        bar.add(fileMenu);

        bar.add(createViewMenu());

        // add the calculations manager menu
        bar.add(calculationsManagerGUI.getMenu());

        bar.add(createToolsMenu());
        bar.add(createPluginsMenu());

        // Current Views menu.
        bar.add(mainDesktop.getCurrentViewsMenu());

        // the help menu
        helpMenu = createHelpMenu();
        bar.add(helpMenu);

        return bar;
    }


    private JMenu createPluginsMenu() {
        JMenu menu = new JMenu("Plugins");

        ArrayList<File> jars = IOUtil.findInClasspath(".jar");

        for (File file : jars) {
            try {
                JarFile f = new JarFile(file);
                Manifest m = f.getManifest();
                if (m == null) {
                    continue;
                }
                Attributes a = m.getMainAttributes();
                String mainClass = a.getValue("Main-Class");
                if (mainClass == null) {
                    continue;
                }

                URL url = new URL("file:///" + file.getAbsolutePath());
                URLClassLoader loader = new URLClassLoader(new URL[] {url},
                        ClassLoader.getSystemClassLoader());
                final Class c = loader.loadClass(mainClass);
                boolean isPlugin = false;
                for (Class i : c.getInterfaces()) {
                    if (i.equals(VisionPlugin.class)) {
                        isPlugin = true;
                    }
                }

                if (!isPlugin) {
                    continue;
                }

                JMenuItem item = new JMenuItem(
                        mainClass.substring(mainClass.lastIndexOf(".") + 1));
                item.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        try {
                            Object plugin = c.newInstance();
                        } catch (Exception ex) {
                            Vision.reportException("Could not open plugin", ex);
                            ex.printStackTrace();
                        }
                    }
                });
                menu.add(item);
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }

        return menu;
    }


    private JMenu createHelpMenu() {
        JMenu helpMenu = new JMenu("Help");
        helpMenu.setMnemonic(KeyEvent.VK_H);

        javax.help.HelpSet helpSet = openHelpSet("help/UserHelp/userhelp.hs");
        if (helpSet != null) {
            javax.help.JHelp helpViewer = new javax.help.JHelp(helpSet);
            JSplitPane split = (JSplitPane) helpViewer.getComponent(0);
            split.setOrientation(JSplitPane.VERTICAL_SPLIT);
            userHelpPane = new JXTaskPane();
            userHelpPane.setTitle("User Help");
            userHelpPane.add(helpViewer);

            JCheckBoxMenuItem item = new JCheckBoxMenuItem("User Help");
            item.addItemListener(new ItemListener() {
                public void itemStateChanged(ItemEvent e) {
                    if (e.getStateChange() == ItemEvent.SELECTED) {
                        controlPane.add(userHelpPane, 0);
                    } else {
                        controlPane.remove(userHelpPane);
                    }
                    controlPane.getParent().validate();
                }
            });
            helpMenu.add(item);
        }

        helpSet = openHelpSet("help/DeveloperHelp/developerhelp.hs");
        if (helpSet != null) {
            javax.help.JHelp helpViewer = new javax.help.JHelp(helpSet);
            JSplitPane split = (JSplitPane) helpViewer.getComponent(0);
            split.setOrientation(JSplitPane.VERTICAL_SPLIT);
            developerHelpPane = new JXTaskPane();
            developerHelpPane.setTitle("Developer Help");
            developerHelpPane.add(helpViewer);

            JCheckBoxMenuItem item = new JCheckBoxMenuItem("Developer Help");
            item.addItemListener(new ItemListener() {
                public void itemStateChanged(ItemEvent e) {
                    if (e.getStateChange() == ItemEvent.SELECTED) {
                        controlPane.add(developerHelpPane, 0);
                    } else {
                        controlPane.remove(developerHelpPane);
                    }
                    controlPane.getParent().validate();
                }
            });
            helpMenu.add(item);
        }

        helpMenu.add(new JSeparator());

        JMenuItem item = new JMenuItem("Print JVM memory status");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                for (int i = 0; i < 5; i++) {
                    Runtime.getRuntime().gc();
                    Runtime.getRuntime().runFinalization();
                }

                double Mb = 1.0 / 1024.0 / 1024.0;
                System.err.printf("Total memory: %.2f Mb\n",
                        Runtime.getRuntime().totalMemory() * Mb);
                System.err.printf("Free memory: %.2f Mb\n",
                        Runtime.getRuntime().freeMemory() * Mb);
                System.err.printf("Used memory: %.2f Mb\n",
                        (Runtime.getRuntime().totalMemory() -
                                Runtime.getRuntime().freeMemory()) * Mb);
            }
        });
        helpMenu.add(item);

        helpMenu.add(new JSeparator());

        item = new JMenuItem("About VISION");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                startupWindow.setProgress("Click to close...");
                startupWindow.show();
            }
        });
        helpMenu.add(item);

        return helpMenu;
    }


    private JMenu createToolsMenu() {
        JMenu toolsMenu = new JMenu("Tools");
        toolsMenu.setMnemonic(KeyEvent.VK_T);

        JMenuItem item;

        item = new JMenuItem("Start Command Receiver");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                try {
                    new CommandReceiver().start();
                } catch (IOException ex) {
                    reportException("Could not start the Command Reciever", ex);
                }
            }
        });
        toolsMenu.add(item);
        toolsMenu.add(new JSeparator());

        item = new JMenuItem(new NeuronIdentification());
        toolsMenu.add(item);

        item = new JMenuItem(new WhiteNoiseDatasetAnalysis());
        toolsMenu.add(item);
        
        item = new JMenuItem(new StaAnalysis());
        toolsMenu.add(item);

        item = new JMenuItem(new RawMovieDatasetAnalysis());
        toolsMenu.add(item);
        
        item = new JMenuItem(new MappingAnalysis());
        toolsMenu.add(item);
        
        item = new JMenuItem(new ConcatenatedMappingAnalysis());
        toolsMenu.add(item);

        toolsMenu.add(new JSeparator());

        item = new JMenuItem(new CreatePinkNoiseMovie());
        item.setMnemonic(KeyEvent.VK_P);
        toolsMenu.add(item);

        item = new JMenuItem(new CreateDriftingSinusoidMovie());
        toolsMenu.add(item);

        item = new JMenuItem(new OpenWaveformViewer());
        toolsMenu.add(item);



        item = new JMenuItem(new OpenNeuronWaveformViewer());
        toolsMenu.add(item);

        return toolsMenu;
    }


    /** This method is public as an implementation side effect. Do NOT call. */
    public void actionPerformed(ActionEvent event) {
        String command = event.getActionCommand();

        if (command.equals("ExportDesktop")) {
            GraphicsIO.saveImage(mainDesktop, mainFrame, "Export Desktop");
        } else if (command.equals("ExportView")) {
            JInternalFrame f = mainDesktop.getSelectedFrame();
            if (f != null) {
                GraphicsIO.saveImage(f, mainFrame, "Export View");
            }
        } else if (command.equals("ExportScreenShot")) {
            GraphicsIO.saveImage(mainFrame.getRootPane(), mainFrame, "Export Screenshot");
        } else if (command.equals("Open Location")) {
            int returnValue = fileDialog.showOpenDialog(
                    Vision.getInstance().getMainFrame());
            if (returnValue == JFileChooser.APPROVE_OPTION) {
                try {
                    String folder = fileDialog.getSelectedFile().getAbsolutePath();
                    dataManager.setFolder(folder);
                    visionParametersGroup.put("folder", folder);
                    config.setParameterGroup("Vision Parameters", visionParametersGroup);
                } catch (Exception e) {
                    reportException(e);
                }
            }
        } else if (command.equals("Refresh Open Folder")) {
            try {
                dataManager.setFolder( (String) visionParametersGroup.get("folder"));
            } catch (IOException ex) {
                reportException("Cannot refresh the current folder", ex);
            }
        } else if (command.equals("Exit")) {
            exitVision();
        }
    }


    public Component getFrameOfClass(Class c) {
        if (mainDesktop == null) {
            return null;
        }
        JInternalFrame[] f = mainDesktop.getAllFrames();
        for (int i = 0; i < f.length; i++) {
            if (f[i].getContentPane().getComponent(0).getClass().equals(c)) {
                return (Component) f[i].getContentPane().getComponent(0);
            }
        }
        return null;
    }


    /**
     * This method creates a new frame inside the VISION window. The new frame will
     * contain the specified <tt>view</tt> component. The <tt>Control</tt> tab in the
     * Vision UI will display the specified <tt>controller</tt> when this window gets
     * selected.
     *
     * @param view the content of the new frame
     * @param controller the controller of the frame
     * @param name the name of the new frame
     */
    public JInternalFrame createFrame(
            JComponent view, JComponent controller, JMenu[] menu, String name) {

        final JInternalFrame frame = new JInternalFrame(name, true, true, true, true);
        frame.getContentPane().add(view);
        frame.setDefaultCloseOperation(JInternalFrame.DO_NOTHING_ON_CLOSE);
        frame.addInternalFrameListener(new InternalFrameAdapter() {
            public void internalFrameClosing(InternalFrameEvent event) {
                removeFrame(frame);
            }
        });

        controllers.put(frame, controller);
        menus.put(frame, menu);

        mainDesktop.addFrame(frame, name);

        frame.setVisible(true);

        // Try to make it the selected one.
        try {
            frame.setSelected(true);
        } catch (PropertyVetoException e) {}

        mainDesktop.updateLayout();

        return frame;
    }


    /**
     * Removes a JInternalFrame from the Vision window. The method is thread safe.
     *
     * @param frame JInternalFrame
     */
    public void removeFrame(final JInternalFrame frame) {
        if (SwingUtilities.isEventDispatchThread()) {
            removeFrameImpl(frame);
        } else {
            try {
                SwingUtilities.invokeAndWait(new Runnable() {
                    public void run() {
                        removeFrameImpl(frame);
                    }
                });
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }
    }


    private void removeFrameImpl(final JInternalFrame frame) {
        if (frame.getContentPane().getComponentCount() == 1 &&
                frame.getContentPane().getComponent(0) instanceof Closable) {
            if (! ( (Closable) frame.getContentPane().getComponent(0)).canClose()) {
                return;
            }
        }

        JMenu[] menu = (JMenu[]) menus.get(currentFrame);
        if (menu != null) {
            for (int i = 0; i < menu.length; i++) {
                menuBar.remove(menu[i]);
            }
        }

        controllers.remove(frame);
        menus.remove(frame);
        mainDesktop.removeFrame(frame);

        frame.getContentPane().removeAll();
        frame.dispose();

        {
            MouseListener[] l = frame.getMouseListeners();
            for (int i = 0; i < l.length; i++) {
                frame.removeMouseListener(l[i]);
            }
        }
        {
            MouseMotionListener[] l = frame.getMouseMotionListeners();
            for (int i = 0; i < l.length; i++) {
                frame.removeMouseMotionListener(l[i]);
            }
        }
        {
            PropertyChangeListener[] l = frame.getPropertyChangeListeners();
            for (int i = 0; i < l.length; i++) {
                frame.removePropertyChangeListener(l[i]);
            }
        }
        {
            InternalFrameListener[] l = frame.getInternalFrameListeners();
            for (int i = 0; i < l.length; i++) {
                frame.removeInternalFrameListener(l[i]);
            }
        }

        frame.getContentPane().removeAll();

        // select another frame
        JInternalFrame[] f = mainDesktop.getAllFrames();
        if (f != null && f.length != 0) {
            selectionHappened(f[0]);
        } else {
            selectionHappened(null);
        }
    }


    public void setFrameTitle(JComponent controller, String newTitle) {
        for (JInternalFrame frame : controllers.keySet()) {
            if (controller == controllers.get(frame)) {
                frame.setTitle(newTitle);
                updateTitle();
                return;
            }
        }
        throw new IllegalArgumentException("Unknown controller.");
    }


    /**
     * Sets the currently selected <tt>JInternalFrame</tt>. The controller component
     * will update accordingly.
     */
    public void selectionHappened(JInternalFrame frame) {
        viewManager.setCurrentController(
                controllers.get(frame), frame == null ? "" : frame.getTitle());

        // remove the old menu (if present)
        if (currentFrame != null) {
            JMenu[] menu = (JMenu[]) menus.get(currentFrame);
            if (menu != null) {
                for (int i = 0; i < menu.length; i++) {
                    menuBar.remove(menu[i]);
                }
            }
        }

        // add the new menu (if present)
        JMenu[] menu = (JMenu[]) menus.get(frame);
        if (menu != null) {
            for (int i = 0; i < menu.length; i++) {
                menuBar.add(menu[i], menuBar.getComponentIndex(helpMenu) - 1);
            }
        }

        mainFrame.validate();
        mainFrame.repaint();

        currentFrame = frame;

        // set the title
        updateTitle();
    }


    int x, y, w, h;

    public void initGUI(StartupWindow startupWindow) {
        guiEnabled = true;

        startupWindow.setProgress("Set UI Defaults...");
        // set some UI Defaults
        UIDefaults defaults = UIManager.getDefaults();
        defaults.put("ScrollPane.background", defaults.get("Tree.textBackground"));
        defaults.put("Viewport.background", defaults.get("Tree.background"));
        defaults.put("Table.selectionBackground", new ColorUIResource(184, 207, 240));
        defaults.put("Table.selectionForeground", new ColorUIResource(0, 0, 0));
        defaults.put("ScrollBar.width", new Integer(15));

        startupWindow.setProgress("Redirect Standard Output...");
        // redirect the output to a folder-specific file
        FileOutputStream outputStream;
        try {
            // get the name of the output log file
            String fileName = VisionParams.VISION_OUTPUT_FILE;
            
            // write the output text file to the user's home directory
            File outputFile = 
                new File(System.getProperty("user.home") + File.separator + fileName);
            outputStream = new FileOutputStream(outputFile, true);
            fileStandardOutput = new PrintStream(outputStream, true);

            // create the terminal
            JTextArea outputTextArea = new JTextArea();
//			Color c = UIManager.getColor("Desktop.background");

            outputTextArea.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
            terminalComponent = new JScrollPane(outputTextArea);
            outputTextArea.setFont(new Font("Monospaced", Font.PLAIN, 11));
            System.setOut(new UIPrintStream(outputStream, outputTextArea));
            System.setErr(System.out);
        } catch (IOException e) {
            reportFatalException(
                    "The vision-output.txt file cannot be opened, created or written.", e);
        }

        startupWindow.setProgress("Read Vision Parameters...");
        visionParametersGroup = config.getParameterList("Vision Parameters");
        try {
            x = Integer.parseInt( (String) visionParametersGroup.get("X"));
            y = Integer.parseInt( (String) visionParametersGroup.get("Y"));
            w = Integer.parseInt( (String) visionParametersGroup.get("W"));
            h = Integer.parseInt( (String) visionParametersGroup.get("H"));
            PlotPanel.saveFolder = (String) visionParametersGroup.get("Plot Panel Folder");
    
        } catch (Exception e) {
            x = 20;
            y = 20;
            w = 920;
            h = 570;
        }
        
        //get maximum possible screen size, taking into account multiple monitors
        Rectangle totalBounds = new Rectangle();
        GraphicsEnvironment ge = GraphicsEnvironment.
        getLocalGraphicsEnvironment();
        GraphicsDevice[] gs =
            ge.getScreenDevices();
        for (int j = 0; j < gs.length; j++) { 
            GraphicsDevice gd = gs[j];
            GraphicsConfiguration[] gc =
                gd.getConfigurations();
            for (int i=0; i < gc.length; i++) {
                totalBounds =
                    totalBounds.union(gc[i].getBounds());
            }
        }
        
        //constrain to available area.
        if(x<0) x=0;
        if(y<0) y=0;
        if(x+w > totalBounds.width)
            w = Math.max(totalBounds.width - x, 1000);
        if(y+h > totalBounds.height)
            h = Math.max(totalBounds.height - y, 500);
        
        

        startupWindow.setProgress("Create the Data Manager...");
        dataManager = new DataManager();
        try {
            String oldFolder = (String) visionParametersGroup.get("folder");
            if (oldFolder != null && oldFolder.trim().length() != 0) {
                dataManager.setFolder(oldFolder);
            }
        } catch (IOException e) {
            reportException(e);
        }

        startupWindow.setProgress("Create the View Manager...");
        viewManager = new ViewManager();

        startupWindow.setProgress("Create the Calculation Manager GUI...");
        calculationsManagerGUI = new CalculationManagerGUI(calculationsManager);

        // Make the dialogs.
//		exportDialog = new ExportDialog();

        // create the file dialog
        startupWindow.setProgress("Create the File Dialog...");
        fileDialog = new JFileChooser();
        fileDialog.removeChoosableFileFilter(fileDialog.getAcceptAllFileFilter());
        fileDialog.setApproveButtonText("Select");
        fileDialog.setDialogTitle("Select");
        fileDialog.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        fileDialog.setApproveButtonText("Open");
        fileDialog.setDialogTitle("Open Location");

        // Make the main desktop.
        startupWindow.setProgress("Create the Main Desktop...");
        mainDesktop = new Desktop();
        mainDesktop.addSelectionListener(this);

        // make ther munu bar
        startupWindow.setProgress("Create Menus...");
        menuBar = makeMenu();

        startupWindow.setProgress("Create the Main Frame...");
        mainFrame = new JFrame();
        mainFrame.setIconImage(VisionImage.getIcon("icon.png").getImage());
        mainFrame.getContentPane().setLayout(new BorderLayout());
        mainFrame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
        mainFrame.setJMenuBar(menuBar);
        mainFrame.addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent e) {
                exitVision();
            }
        });

        controlPane = new JXTaskPaneContainer();
        controlPane.add(calculationsManagerGUI);
        controlPane.add(viewManager);
        controlPane.add(dataManager);

        splitPane = new JSplitPane(
                JSplitPane.HORIZONTAL_SPLIT, true, null, null);

        setTerminal(1);

        splitPane.setDividerLocation(mainDividerSize);
        splitPane.setOneTouchExpandable(true);
        splitPane.setContinuousLayout(true);
        splitPane.setDividerSize(6);

        // make the frame
        mainFrame.add(splitPane, BorderLayout.CENTER);
        statusBar = new StatusBar();
        mainFrame.add(statusBar, BorderLayout.SOUTH);
        mainFrame.setBounds(x, y, w, h);

        SwingUtilities.updateComponentTreeUI(mainFrame);
        mainFrame.setVisible(true);

        sendMessage("Welcome to VISION!");
    }

    private void setTerminal(int i) {
        int n = splitPane.getDividerLocation();

        if (i == 0) {
            splitPane.setLeftComponent(new JScrollPane(controlPane));
            splitPane.setRightComponent(mainDesktop);
        } else if (i == 1) {
            JSplitPane split = new JSplitPane(
                    JSplitPane.VERTICAL_SPLIT, true, new JScrollPane(controlPane),
                    terminalComponent);
            split.setOneTouchExpandable(true);
            split.setDividerLocation(h - 200);
            split.setResizeWeight(1);
            split.setDividerSize(6);

            splitPane.setLeftComponent(split);
            splitPane.setRightComponent(mainDesktop);
        } else if (i == 2) {
            JSplitPane split = new JSplitPane(
                    JSplitPane.VERTICAL_SPLIT, true, mainDesktop, terminalComponent);
            split.setOneTouchExpandable(true);
            split.setDividerLocation(h - 200);
            split.setResizeWeight(1);
            split.setDividerSize(6);

            splitPane.setLeftComponent(new JScrollPane(controlPane));
            splitPane.setRightComponent(split);
        }

        splitPane.setDividerLocation(n);
        mainFrame.validate();
    }


    public void exitVision() {
        visionParametersGroup.put("X", "" + mainFrame.getX());
        visionParametersGroup.put("Y", "" + mainFrame.getY());
        visionParametersGroup.put("W", "" + mainFrame.getWidth());
        visionParametersGroup.put("H", "" + mainFrame.getHeight());
        visionParametersGroup.put("Plot Panel Folder", PlotPanel.saveFolder);
        try {
            config.setParameterGroup("Vision Parameters", visionParametersGroup);
        } catch (Exception e) {
            reportException(e);
        }

        if (calculationsManager.isCurrentlyCalculating()) {
            int result = JOptionPane.showConfirmDialog(
                    mainFrame, "The " +
                    calculationsManager.getCurrentCalculationName() +
                    " calculation is currently going on ! " +
                    "Do you want to exit anyway?", "Warning",
                    JOptionPane.YES_NO_OPTION, JOptionPane.WARNING_MESSAGE);
            if (result == JOptionPane.YES_OPTION) {
                System.exit(0);
            }
        } else {
            JInternalFrame[] frames = mainDesktop.getAllFrames();
            for (int i = 0; i < frames.length; i++) {
                if (frames[i].getContentPane().getComponentCount() == 1 &&
                        frames[i].getContentPane().getComponent(0) instanceof Closable) {
                    if (! ( (Closable) frames[i].getContentPane().getComponent(0)).
                            canClose()) {
                        return;
                    }
                }
            }
            System.exit(0);
        }
    }


    public static boolean isGUIBased() {
        return guiEnabled;
    }


    public static boolean isConsoleBased() {
        return!guiEnabled;
    }


    public static final void reportIncorectParams(LinkedHashMap<String, String> params) {
        System.out.println(
                "Incorrect number of command line arguments, " + params.size() +
        " required:");
        int i = 1;
        for (String item : params.keySet()) {
            System.out.println(" " + i + ") " + item);
        }
    }


    public static void reportException(Exception e) {
        reportException("", e);
    }


    public static void reportException(String message, Exception e) {
        String name = "EXCEPTION: " + e.getClass().getName();
        String msg = message + "\n" + e.getClass().getName() + ":\n" + e.getMessage();

        if (Vision.isGUIBased()) {
            int result = JOptionPane.showOptionDialog(
                    Vision.getInstance().getMainFrame(),
                    msg, name,
                    JOptionPane.YES_OPTION,
                    JOptionPane.ERROR_MESSAGE,
                    null,
                    new Object[] {"Ok", "Print Stack Trace"}
                    , null);
            if (result == 1) {
                System.out.println("========== Begin Stack Trace ===========");
                e.printStackTrace();
                System.out.println("=========== End Stack Trace ============");
            }
        }

        if (Vision.isConsoleBased()) {
            System.out.println(name);
            System.out.println(message);

            System.out.println("========== Begin Stack Trace ===========");
            e.printStackTrace();
            System.out.println("=========== End Stack Trace ============");
        }
    }


    public static void reportFatalException(Exception e) {
        reportFatalException("", e);
    }


    public static void reportFatalException(String message, Exception e) {
        String name = "FATAL EXCEPTION: " + e.getClass().getName();
        String msg = message + "\n" + e.getClass().getName() + ":\n" + e.getMessage() +
        "\n\nVISION WILL EXIT";

        if (Vision.isGUIBased()) {
            int result = JOptionPane.showOptionDialog(
                    Vision.getInstance().getMainFrame(),
                    msg, name,
                    JOptionPane.YES_OPTION,
                    JOptionPane.ERROR_MESSAGE,
                    null,
                    new Object[] {"Ok", "Print Stack Trace"}
                    , null);
            if (result == 1) {
                System.out.println("========== Begin Stack Trace ===========");
                e.printStackTrace();
                System.out.println("=========== End Stack Trace ============");
            }
        }

        if (Vision.isConsoleBased()) {
            System.out.println(name);
            System.out.println(message);

            System.out.println("========== Begin Stack Trace ===========");
            e.printStackTrace();
            System.out.println("=========== End Stack Trace ============");
        }

        System.exit(1);
    }


    /**
     * The main() method which starts the application.
     */
    public static void main(String args[]) throws Exception {  	

        // Check the JDK version number.
        String version = System.getProperty("java.version");
        if ( (version == null) ||
                (version.startsWith("0.")) ||
                (version.startsWith("1.0")) ||
                (version.startsWith("1.1")) ||
                (version.startsWith("1.2")) ||
                (version.startsWith("1.3")) ||
                (version.startsWith("1.4")) ||
                (version.startsWith("1.5")) ){ //||
                //(version.startsWith("1.6")) ) {
            throw new IllegalStateException(
            "Vision requires Java JRE (or JDK) version 7.0 (1.7) or higher.");
        }

        VisionJSAP jsap = new VisionJSAP(Vision.class.getName(),
                new com.martiansoftware.jsap.Parameter[] {
                    new FlaggedOption("config", JSAP.STRING_PARSER, "config.xml", JSAP.NOT_REQUIRED, 'c', "config", "Configuration file.")
                }
        );

        JSAPResult parsedArgs = jsap.parse(args);

        startupWindow = new StartupWindow();
        startupWindow.show();

        startupWindow.setProgress("Set Look & Feel...");
        try {
            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
        } catch (Exception e) {
            e.printStackTrace();
        }

        // Linux sets LookAndFeel to GTK. Unfortunately, the window manager for GTK is not very compatible;
        // it fails to trigger callbacks in multiple ways. E.g. setSelection(true) doesn't seem to call 
        // activateFrame(). This messes up Vision's GUI interconnections. So revert to Metal for now.
        if (UIManager.getLookAndFeel().getClass().getName() == "com.sun.java.swing.plaf.gtk.GTKLookAndFeel") {
            UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
        }

        
        // for a period of time, we used quaqua to display the window settings
        // however, it was hard to maintain the UI to keep up with changes from Apple
        // consequently, we don't use it anymore. you can give it a try if you want...
        // instead, i recommend migrating to SWT for everything. but that would take
        // some time and messing around.
        //if (System.getProperty("os.name").matches("Mac OS X")) {
        //    UIManager.setLookAndFeel("ch.randelshofer.quaqua.QuaquaLookAndFeel");
        //}
        
        Toolkit.getDefaultToolkit().getSystemEventQueue().push(new VisionEventQueue());


        startupWindow.setProgress("Initialize...");

        // Now create the application.
        Vision vision = Vision.getInstance();

        // Start logging immediately to the vision file
        // This will get set again in initGui but for now it lets us
        // catch messages that happen between now and the vision.initGUI call
        FileOutputStream outputStream;
        try {
            // get the name of the output log file
            String fileName = VisionParams.VISION_OUTPUT_FILE;

            // write the output text file to the user's home directory
            File outputFile =
                    new File(System.getProperty("user.home") + File.separator + fileName);
            outputStream = new FileOutputStream(outputFile, true);
            System.setOut(new PrintStream(outputStream));
            System.setErr(System.out);
        } catch (IOException e) {
            reportFatalException(
                    "The vision-output.txt file cannot be opened, created or written.", e);
        }

        // Try reading the config file
        try {
            vision.setConfig(parsedArgs.getString("config"));
        } catch (IllegalArgumentException ex) {
            startupWindow.hide();
            System.exit(1);
            return;
        }

        // If everything is ok with the config, we can initialize the GUI
        vision.initGUI(startupWindow);
        vision.updateTitle();
        vision.controlPane.setBorder(null);

        startupWindow.setProgress("Welcome to Vision");

        try {
            Thread.sleep(1000);
        } catch (InterruptedException ex) {}
        startupWindow.hide();
    }

    
    public static String buildDateString() {
        try {
            Manifest mf = getManifest();
            if (mf == null) return null;    		
            Attributes atts = mf.getMainAttributes();
            return atts.getValue("Built-Date");
        } catch (IOException e) {
            return null;
        }
    }
    
    public static Manifest getManifest() throws IOException {
        URL res = Vision.class.getResource(Vision.class.getSimpleName() + ".class");
        URLConnection conn = res.openConnection();
        if (!(conn instanceof JarURLConnection)) return null;
        JarURLConnection jconn = (JarURLConnection) conn;
        return jconn.getManifest();
    }

}
