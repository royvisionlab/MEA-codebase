package edu.ucsc.neurobiology.vision.gui;

import java.io.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.event.*;
import javax.swing.*;
import javax.swing.tree.*;
import java.util.*;

import org.jdesktop.swingx.*;
import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.dataview.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.io.chunk.*;
import edu.ucsc.neurobiology.vision.neuronviewer.*;
import edu.ucsc.neurobiology.vision.stimulus.*;
import edu.ucsc.neurobiology.vision.util.*;



/**
 * This class manages a set of Objects which represent data sets.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 * @author Matthew Grivich, The Salk Institute
 */
public class DataManager
extends JXTaskPane {

    private DefaultMutableTreeNode rootNode;
    private DefaultTreeModel treeModel;
    private JTree dataTree;

    public void setFolder(String folder) throws IOException {
        File f = new File(folder);
        if (!IOUtil.isValidFolder(f)) {
            folder = System.getProperty("user.dir");
            f = new File(folder);
        }
        if (IOUtil.isValidFolder(f)) {
            rootNode = new DefaultMutableTreeNode(folder);
            treeModel = new DefaultTreeModel(rootNode, true);
            buildTree(new File(folder), treeModel, rootNode);
            dataTree.setModel(treeModel);
            dataTree.setToggleClickCount(1000); // do not expand at double-click
            dataTree.setShowsRootHandles(true);
        }
    }
    
    
    /**
     * Create the path tree components and add them to the GUI
     * 
     * @param directory
     * @param model
     * @param treeFolder
     * @throws IOException
     */
    public static void buildTree(File directory, DefaultTreeModel model,
            DefaultMutableTreeNode treeFolder) throws IOException {

        File[] children = directory.listFiles(new NoOSXCrap());
        Arrays.sort(children);
        for (int i = 0; i < children.length; i++) {
            String name = children[i].getName();
            if (children[i].isDirectory()) {
                DefaultMutableTreeNode newFolder = new DefaultMutableTreeNode(name, true);
                model.insertNodeInto(newFolder, treeFolder, treeFolder.getChildCount());
                
            } else if (children[i].isFile()) {
                model.insertNodeInto(
                        new DefaultMutableTreeNode(name, false), treeFolder,
                        treeFolder.getChildCount());
            }
        }
    }
    
    
    private static class NoOSXCrap implements FilenameFilter {
        public boolean accept(File dir, String name) {
            if (name.equals(".DS_Store")) return false;
            if (name.startsWith("._"))    return false;
            
            return true;
        }
    }


    private String getFileNameForPath(TreePath path) {
        String name = "";
        for (int i = 0; i < path.getPathCount(); i++) {
            name += path.getPathComponent(i).toString();
            if (i != path.getPathCount() - 1) {
                name += File.separator;
            }
        }
        return name;
    }


    public DataManager() {
        this.setTitle("Data Manager");
        this.setScrollOnExpand(true);

//		super(new BorderLayout());

        dataTree = new JTree();
        dataTree.setEditable(false);
//		dataTree.setCellRenderer(new R());
        dataTree.setRootVisible(true);
        dataTree.setShowsRootHandles(true);
        dataTree.setScrollsOnExpand(true);
        dataTree.setRowHeight(15);
        dataTree.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
        dataTree.getSelectionModel().
        setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
        ToolTipManager.sharedInstance().registerComponent(dataTree);

        TreePane pane = new TreePane(dataTree);
        dataTree.addTreeExpansionListener(pane);
        pane.setSize(0, 400);
        pane.setPreferredSize(pane.getSize());

        this.add(pane, BorderLayout.CENTER);

//		( (DefaultTreeCellRenderer) dataTree.getCellRenderer()).setLeafIcon(null);
//		( (DefaultTreeCellRenderer) dataTree.getCellRenderer()).setOpenIcon(null);
//		( (DefaultTreeCellRenderer) dataTree.getCellRenderer()).setClosedIcon(null);

//		JPanel buttonPanel = new JPanel();
//		buttonPanel.setBorder(null);
//		JButton b = new JButton(VisionImage.getIcon("RotateCW.gif"));
//		b.setBorder(null);
//		b.setPreferredSize(new Dimension(20, 20));
//		buttonPanel.add(b);

//		this.add(buttonPanel, BorderLayout.NORTH);

        dataTree.addKeyListener(new KeyAdapter() {
            public void keyPressed(KeyEvent e) {
                if (e.getKeyCode() == KeyEvent.VK_ENTER) {
                    TreePath path = dataTree.getSelectionPath();
                    if (path != null) {
                        open(path);
                    }
                }
            }
        });

        dataTree.addMouseListener(new MouseAdapter() {
            public void mouseClicked(MouseEvent event) {
                if (event.getClickCount() == 2) {
                    TreePath path = dataTree.getPathForLocation(event.getX(), event.getY());
                    if (path == null) {
                        return;
                    }
                    open(path);
                }
            }
        });

        dataTree.addMouseMotionListener(new MouseMotionAdapter() {
            public void mouseMoved(MouseEvent event) {
                TreePath path = dataTree.getPathForLocation(event.getX(), event.getY());
                if (path == null) {
                    return;
                }

                String name = path.getLastPathComponent().toString();
                if (name.endsWith(".sta")) {
//					dataTree.setToolTipText("Tool tip");
//					ToolTipManager.sharedInstance().setEnabled(true);
//					ToolTipManager.sharedInstance().setInitialDelay(0);
//					ToolTipManager.sharedInstance().setsetInitialDelay(0);
                } else {
                    dataTree.setToolTipText("");
                }
            }
        });
    }


    private void togglePath(TreePath path) {
        if (dataTree.isCollapsed(path)) {
            dataTree.expandPath(path);
        } else {
            dataTree.collapsePath(path);
        }
    }
    
    
    /**
     * Try to open the file, if a folder try to open something inside.  
     * If a folder without an obvious thing to open inside, toggle the tree path instead.
     * 
     * @param path
     */
    private void open(final TreePath path) {
        new Thread() {
            public void run() {
                String name = getFileNameForPath(path);
                Vision.getInstance().sendMessage(
                        "Loading " + new File(name).getName() + " ...");
                DefaultMutableTreeNode last = ( (DefaultMutableTreeNode) path.
                        getLastPathComponent());

                try {
                    if (last.getAllowsChildren()) { // Directory
                        String pFileName =
                            name + File.separator + new File(name).getName() + ".params";
                        if (IOUtil.isValidFile(pFileName)) { // Open params file if it's under directory
                            new NeuronViewer(pFileName, true, Vision.getConfig());
                        } else { // Otherwise open the directory in GUI
                            String rawFileName = 
                                name + File.separator + new File(name).getName() + "000.bin";
                            if (IOUtil.isValidFile(rawFileName)) {
                                new RawDataView(new File(name));
                            } else {
                                togglePath(path);
                            }
                        }
                    } else if (name.endsWith(".params")) {
                        new NeuronViewer(name, true, Vision.getConfig());
                    } else if (name.endsWith(".neurons") || name.endsWith(".neurons-raw")) {
                        NeuronFile f = new NeuronFile(name);
                        System.out.println(f);
                        System.err.println("Neurons: " + f.getNumberOfNeurons());
                        f.close();
                        new CrossCorrelator(name);
                    }  else if (name.endsWith(".prj")) {
                        ProjectionsFile f = new ProjectionsFile(name);
                        IOUtil.printPublicFields(f.getHeader());
                        f.close();
                    } else if (name.endsWith(".movie")) {
                        new MovieView(new WhiteNoiseMovie(name, StringUtil.removeExtension(name) + ".globals"));
                    } else if (name.endsWith(".rawMovie")) {
                        MovieView movieView = new MovieView(new RawMovie(name,
                                StringUtil.removeExtension(name) + ".globals"));
//						movieView.saveMovie();

//						MovieView movieView = new MovieView(new ReversingGratingMovie(320, 1, 1, 1,
//						0.48, 30, 0, 7.5));
//						movieView.saveFrames();

                    } else if (name.endsWith(".spikes") || name.endsWith(".spikes-noise")) {
                        SpikeFile f = new SpikeFile(name);
                        new SpikeView(f);
                        new SpikeListView(f);
                        System.err.println("Spike File: \n================");
                        IOUtil.printPublicFields(f.getHeader());
                    } else if (name.endsWith(".bin")) {
                        new RawDataView(new File(name));
                    } else if (name.endsWith(".model")) {
                        ClusteringModelFile f = new ClusteringModelFile(name);
                        System.err.println("Model File: \n================");
                        IOUtil.printPublicFields(f.getUserHeader());
                    } else if (name.endsWith(".cov") || name.endsWith(".wcov") || name.endsWith(".ncov")) {
                        CovarianceFile f = new CovarianceFile(name);
                        System.err.println("Covarinace File: \n================");
                        IOUtil.printPublicFields(f.getHeader());
                        new CovarianceView(f);
                    } else if (name.endsWith(".globals")){
                        GlobalsFile gFile = new GlobalsFile(name, GlobalsFile.READ);
                        System.out.println(gFile.toString());
                    } else{
                        System.out.println("Unknown file type.");
                    }
                } catch (IOException e) {
                    Vision.reportException(e);
                }

                Vision.getInstance().sendMessage("Done.");
            }
        }.start();
    }


    class TreePane extends JScrollPane implements TreeExpansionListener {

        public TreePane(Component component) {
            super(component);
        }

        public void treeExpanded(TreeExpansionEvent e) {
            TreePath path = e.getPath();
            File directory = new File(getFileNameForPath(path));
            DefaultMutableTreeNode node = (DefaultMutableTreeNode) path.getLastPathComponent();


            //if we have not already opened the branch
            if (!node.children().hasMoreElements()) {
                try {
                    buildTree(directory, treeModel, node);
                } catch (IOException ex) {
                    ex.printStackTrace();
                }
            }
        }

        public void treeCollapsed(TreeExpansionEvent e) {

        }
    }



}
