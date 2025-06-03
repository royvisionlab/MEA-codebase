package edu.ucsc.neurobiology.vision.neuronviewer;

import java.util.*;
import java.util.List;

import java.awt.*;
import java.awt.datatransfer.*;
import java.awt.dnd.*;
import java.awt.event.*;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.tree.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class InteractiveTree
    extends JTree implements DragSourceListener, DragGestureListener, DropTargetListener,
    MouseListener, KeyListener {
    
    private static final long serialVersionUID = 1L;

    public static final String ROOT_NAME = "All";
    
    private DefaultTreeModel model;
    private DragSource dragSource;
    private DropTarget dropTarget;
    private ChangeListener changeListener;
    private ArrayList<CalculationAction> calculationActions;
    
    public InteractiveTree(final DefaultTreeModel model) {
        super(model);
        super.setRootVisible(true);
        super.setShowsRootHandles(false);
        super.setScrollsOnExpand(true);
        super.setRowHeight(14);
        super.setEditable(true);
        super.getSelectionModel().
            setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
        ( (DefaultTreeCellRenderer)this.getCellRenderer()).setLeafIcon(null);
        ( (DefaultTreeCellRenderer)this.getCellRenderer()).setOpenIcon(null);
        ( (DefaultTreeCellRenderer)this.getCellRenderer()).setClosedIcon(null);
        this.getCellEditor().addCellEditorListener(new CellEditorListener() {
            public void editingStopped(ChangeEvent e) {
                if (changeListener != null) {
                    changeListener.stateChanged(new ChangeEvent(this));
                }
            }
            
            
            public void editingCanceled(ChangeEvent e) {
            }
        });
        
        this.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
        this.model = model;
        model.addTreeModelListener(new TreeModelListener() {
            public void treeNodesChanged(TreeModelEvent e) {
                if (changeListener != null) {
                    changeListener.stateChanged(new ChangeEvent(this));
                }

//                sortFolder(
//                    (DefaultMutableTreeNode) model.getRoot(),
//                    model, true
//                    );
            }


            public void treeNodesInserted(TreeModelEvent e) {}
            public void treeNodesRemoved(TreeModelEvent e)  {}

//            public void treeNodesStructureChanged(TreeModelEvent e) {
//                if (changeListener != null) {
//                    changeListener.stateChanged(new ChangeEvent(this));
//                }
//            }


            public void treeStructureChanged(TreeModelEvent e) {
                if (changeListener != null) {
                    changeListener.stateChanged(new ChangeEvent(this));
                }
            }
        });

        calculationActions = new ArrayList<CalculationAction>();
        super.addMouseListener(this);
        super.addKeyListener(this);

        // drag and drop
        dragSource = DragSource.getDefaultDragSource();
        dragSource.createDefaultDragGestureRecognizer(this, DnDConstants.ACTION_COPY_OR_MOVE, this);
        dropTarget = new DropTarget(this, DnDConstants.ACTION_COPY_OR_MOVE, this, true);
    }


    public synchronized void setModel(TreeModel model) {
        super.setModel(model);
        this.model = (DefaultTreeModel) model;

        if (changeListener != null) {
            changeListener.stateChanged(new ChangeEvent(this));
        }
    }

    
    /**
     * Looks for Node matching path, if it doesn't find it then it creates a new node and puts it at the end of 
     * the putative parent folder.
     * @param path
     * @param model
     * @return
     */
    private static DefaultMutableTreeNode insurePath(String path, DefaultTreeModel model) {
        DefaultMutableTreeNode node = (DefaultMutableTreeNode) model.getRoot();

        StringTokenizer st = new StringTokenizer(path, "/");
        st.nextToken(); // throw away the root

        while (st.hasMoreTokens()) {
            String nodeName = st.nextToken();

            boolean exists = false;
            for (int i = 0; i < node.getChildCount(); i++) {
                DefaultMutableTreeNode child =
                    (DefaultMutableTreeNode) node.getChildAt(i);
                if (child.getUserObject().equals(nodeName)) {
                    exists = true;
                    node = child;
                }
            }

            if (!exists) {
                DefaultMutableTreeNode newNode = new DefaultMutableTreeNode(nodeName);
                model.insertNodeInto(newNode, node, node.getChildCount());
                node = newNode;
            }
        }

        return node;
    }


    public static DefaultTreeModel makeModel(LinkedHashMap<Integer, ? extends Object> classLabels) {
        DefaultMutableTreeNode root = new DefaultMutableTreeNode(ROOT_NAME, true);
        DefaultTreeModel model = new DefaultTreeModel(root, true);

        for (Object key : classLabels.keySet()) {
            String classPath = (String) classLabels.get(key);
            DefaultMutableTreeNode parent = insurePath(classPath, model);
            model.insertNodeInto(new DefaultMutableTreeNode(key, false), parent, 0);
        }

        sortFolder(root, model, true);
        return model;
    }


    private static void sortFolder(DefaultMutableTreeNode folder, DefaultTreeModel model, boolean recursive) {
        
        ArrayList<DefaultMutableTreeNode> folderList = new ArrayList<DefaultMutableTreeNode>();
        ArrayList<Integer> leafList = new ArrayList<Integer>();
        ArrayList<DefaultMutableTreeNode> childList = new ArrayList<DefaultMutableTreeNode>();
        
        for (int i = 0; i < folder.getChildCount(); i++) {
            DefaultMutableTreeNode child = (DefaultMutableTreeNode) folder.getChildAt(i);
            Object userObject = child.getUserObject();
            if (userObject instanceof Integer) {
                leafList.add((Integer) userObject);
            } else if (userObject instanceof String) {
                folderList.add(child);
            }
            childList.add(child);
        }
        
        for (int i = 0; i < childList.size(); i++)
            model.removeNodeFromParent(childList.get(i));
        
//        Collections.sort(folderList);
        for (DefaultMutableTreeNode f : folderList) {
            model.insertNodeInto(f, folder, folder.getChildCount());
            if (recursive) sortFolder(f, model, true);
        }
        
        Collections.sort(leafList);
        for (int i = 0; i < leafList.size(); i++) {
            model.insertNodeInto(
                new DefaultMutableTreeNode(leafList.get(i), false), folder,
                folder.getChildCount());
        }
    }


    private void removeNeurons(int[] idList, DefaultMutableTreeNode node) {
        ArrayList<DefaultMutableTreeNode> nodesToRemove = new ArrayList<DefaultMutableTreeNode>();

        // collect the nodes that have to be removed
        int n = node.getChildCount();
        for (int i = 0; i < n; i++) {
            DefaultMutableTreeNode child = (DefaultMutableTreeNode) node.getChildAt(i);
            if (child.getAllowsChildren()) {
                removeNeurons(idList, child);
            } else {
                int id = (Integer) child.getUserObject();
                for (int j = 0; j < idList.length; j++) {
                    if (idList[j] == id) {
                        nodesToRemove.add(child);
                    }
                }
            }
        }

        // do the actual removal
        for (int i = 0; i < nodesToRemove.size(); i++) {
            model.removeNodeFromParent( (DefaultMutableTreeNode) nodesToRemove.get(i));
        }
    }


    public static String pathToString(TreePath folderPath) {
        String path = "";

        for (int i = 0; i < folderPath.getPathCount(); i++) {
            if (i != 0) path += "/";
            path += (String) ( (DefaultMutableTreeNode) folderPath.
                              getPathComponent(i)).getUserObject();
        }

        return path;
    }


    public static IntegerList getNeuronsInClass(
        DefaultMutableTreeNode folder, IntegerList neuronList, boolean includeSubclasses) {

        for (int i = 0; i < folder.getChildCount(); i++) {
            DefaultMutableTreeNode child = (DefaultMutableTreeNode) folder.getChildAt(i);
            if (child.getAllowsChildren()) {
                if (includeSubclasses) {
                    neuronList = getNeuronsInClass( (DefaultMutableTreeNode) folder.
                        getChildAt(i), neuronList, includeSubclasses);
                }
            } else {
                neuronList.add( (Integer) child.getUserObject());
            }
        }

        return neuronList;
    }


    public synchronized void setClassLabels(int[] idList, TreePath folderPath) {
        // The new folder has been created by the caller, but hasn't been placed into the tree model.
        // We must retrieve the new folder and parent folders from folderPath, and then add the new folder
        // to the model.  We add new folder to the top of the parent.
        TreePath parentPath = folderPath.getParentPath();
        DefaultMutableTreeNode newFolder    = (DefaultMutableTreeNode) folderPath.getLastPathComponent();
        DefaultMutableTreeNode parentFolder = (DefaultMutableTreeNode) parentPath.getLastPathComponent();
        model.insertNodeInto(newFolder, parentFolder, 0);
        
        // Take neurons out of old folder
        removeNeurons(idList, (DefaultMutableTreeNode) model.getRoot());

        // Put them in new folder
        for (int i = idList.length - 1; i >= 0; i--)
            model.insertNodeInto(new DefaultMutableTreeNode(new Integer(idList[i]), false), newFolder, 0);

        // sort the folder
        sortFolder(newFolder, model, false);

        if (changeListener != null) changeListener.stateChanged(new ChangeEvent(this));
    }


    public synchronized HashMap<Integer,String> getState() {
        DefaultMutableTreeNode root = (DefaultMutableTreeNode) model.getRoot();
        return getState(root, (String) root.getUserObject(), new HashMap<Integer,String>());
    }


    public synchronized HashMap<Integer,String> getState(
            DefaultMutableTreeNode folder, String pathName, HashMap<Integer,String> state) {

        final int nChilds = folder.getChildCount();

        for (int i = 0; i < nChilds; i++) {
            DefaultMutableTreeNode child =
                (DefaultMutableTreeNode) folder.getChildAt(i);

            if (child.getAllowsChildren()) {
                String folderName = (String) child.getUserObject();
                getState(child, pathName + "/" + folderName, state);
            } else {
                Integer id = (Integer) child.getUserObject();
                state.put(id, pathName);
            }
        }

        return state;
    }


    public void addChangeListener(ChangeListener changeListener) {
        this.changeListener = changeListener;
    }


    public void addCalculationAction(CalculationAction action) {
        calculationActions.add(action);
    }


    // DragSourceListener implementation
    public void dragEnter(DragSourceDragEvent dsde) {}
    public void dragOver(DragSourceDragEvent dsde) {}
    public void dropActionChanged(DragSourceDragEvent dsde) {}
    public void dragExit(DragSourceEvent dse) {}
    public void dragDropEnd(DragSourceDropEvent dsde) {}

    // DragGestureListener implementation
    public synchronized void dragGestureRecognized(DragGestureEvent dge) {
        try {
            Point p = dge.getDragOrigin();
            TreePath selectionPath = super.getPathForLocation(p.x, p.y);
            Transferable t = new ObjectTransferable(selectionPath);
            switch (dge.getDragAction()) {
                case DnDConstants.ACTION_COPY: {
                    dragSource.startDrag(dge, DragSource.DefaultCopyDrop, t, this);
                    break;
                }
                case DnDConstants.ACTION_MOVE:
                    try {
                        dragSource.startDrag(dge, DragSource.DefaultMoveDrop, t, this);
                    } catch (InvalidDnDOperationException e) {
                        System.out.println("dragging too fast");
                    }
                    break
                        ;
            }
        } catch (InvalidDnDOperationException e) {
            Vision.reportException(e);
        } catch (ClassNotFoundException e) {
            Vision.reportException(e);
        }
    }


    // DropTargetListener implementation
    public void dragEnter(DropTargetDragEvent dtde) {}

    public void dragOver(DropTargetDragEvent dtde) {
//        DropTargetContext context = dtde.getDropTargetContext();
        Point p = dtde.getLocation();
        TreePath path = super.getClosestPathForLocation(p.x, p.y);
        int row = super.getRowForPath(path);
//        DefaultMutableTreeNode node =
//            (DefaultMutableTreeNode) path.getLastPathComponent();
//        if (node.getAllowsChildren()) {
//            super.setSelectionPath(path);
//        }
        super.scrollRowToVisible(row - 1);
        super.scrollRowToVisible(row);
        super.scrollRowToVisible(row + 1);

//        if (node.getAllowsChildren()) {
//            context.setCursor(DragSource.DefaultCopyDrop);
//        } else {
//            context.setCursor(DragSource.DefaultCopyNoDrop);
//        }
    }


    public void dropActionChanged(DropTargetDragEvent dtde) {}
    public void dragExit(DropTargetEvent dte) {}


    public static DefaultMutableTreeNode[] getChildren(DefaultMutableTreeNode parent) {
        DefaultMutableTreeNode[] children = new DefaultMutableTreeNode[parent.getChildCount()];
        @SuppressWarnings("unchecked") Enumeration<DefaultMutableTreeNode> childEnum = parent.children();
        for (int i = 0; i < children.length; i++) children[i] = childEnum.nextElement();
        return children;
    }


    public synchronized void drop(DropTargetDropEvent dtde) {
        Point p = dtde.getLocation();
        TreePath targetPath = super.getClosestPathForLocation(p.x, p.y);
        if (targetPath == null) return;

        TreePath sourcePath = null;
        try {
            sourcePath = (TreePath) dtde.getTransferable().getTransferData(
                new DataFlavor(DataFlavor.javaJVMLocalObjectMimeType));
        } catch (Exception e) {
            Vision.reportException(e);
            return;
        }
        if (sourcePath == null) return;
        if (targetPath.equals(sourcePath)) return;

        DefaultMutableTreeNode lastTargetNode = (DefaultMutableTreeNode) targetPath.getLastPathComponent();
        if (!lastTargetNode.getAllowsChildren()) return;
        
        DefaultMutableTreeNode lastSourceNode = (DefaultMutableTreeNode) sourcePath.getLastPathComponent();
        model.removeNodeFromParent(lastSourceNode);

        if (dtde.getDropAction() == DnDConstants.ACTION_COPY) {
            // if a copy is going on
            if (lastSourceNode.getAllowsChildren()) {
                for (DefaultMutableTreeNode child : getChildren(lastSourceNode)) {
                    model.removeNodeFromParent(child);
                    if (child.getAllowsChildren()) model.insertNodeInto(child, lastTargetNode, 0);
                    else insertCellIntoSortedPosition(child, lastTargetNode);
                }
            } else {
                insertCellIntoSortedPosition(lastSourceNode, lastTargetNode);
            }
        } else if (dtde.getDropAction() == DnDConstants.ACTION_MOVE) {
            // if a move is going on
            if (lastSourceNode.getAllowsChildren()) {
                model.insertNodeInto(lastSourceNode, lastTargetNode, 0);
            } else {
                insertCellIntoSortedPosition(lastSourceNode, lastTargetNode);
            }
        }

        if (changeListener != null) changeListener.stateChanged(new ChangeEvent(this));
    }
    
    private void insertCellIntoSortedPosition(DefaultMutableTreeNode cell, DefaultMutableTreeNode folder) {
        Integer cellid = (Integer) cell.getUserObject();
        int pos = getCellSortedPosition(cellid, folder);
        model.insertNodeInto(cell, folder, pos);      
    }

    /**
     * Assumes parent is already sorted and puts into position relative to other cells, ignoring
     * folders.
     * @param cell
     * @param parent
     * @return
     */
    private static int getCellSortedPosition(Integer cell, TreeNode parent) {
        @SuppressWarnings("unchecked") Enumeration<DefaultMutableTreeNode> children = parent.children();
        int i = 0;
        while (children.hasMoreElements()) {
            DefaultMutableTreeNode child = children.nextElement();
            Object o = child.getUserObject();
            if (o instanceof Integer && (Integer) o >= cell) break;
            i++;
        }
        return i;
    }
    
    
    // MouseListener implementation
    public void mouseClicked(MouseEvent e) {}
    
    
    public synchronized void mousePressed(MouseEvent e) {
        if (!SwingUtilities.isRightMouseButton(e)) {
            return;
        }
        showPopupMenu(getPathForLocation(e.getPoint().x, e.getPoint().y));
    }

    
    public static void deleteClass(DefaultTreeModel model, DefaultMutableTreeNode folder, DefaultMutableTreeNode parent) {
        clearSubclasses(model, folder, folder);
        ArrayList<DefaultMutableTreeNode> toAdd = new ArrayList<DefaultMutableTreeNode>();
        for (DefaultMutableTreeNode child : getChildren(folder)) {
            toAdd.add(child);
            model.removeNodeFromParent(child);
        }
        for (DefaultMutableTreeNode child : toAdd)
            model.insertNodeInto(child, parent, 0);
        model.removeNodeFromParent(folder);
    }
    
    
    public static void clearSubclasses(DefaultTreeModel model, DefaultMutableTreeNode startFolder, DefaultMutableTreeNode folder) {
        ArrayList<DefaultMutableTreeNode> toAdd = new ArrayList<DefaultMutableTreeNode>();        
        for (DefaultMutableTreeNode child : getChildren(folder)) {
            if (child.getAllowsChildren()) 
                clearSubclasses(model, startFolder, child);
            else
                toAdd.add(child);
            model.removeNodeFromParent(child);
        }
        for (DefaultMutableTreeNode child : toAdd)
            model.insertNodeInto(child, startFolder, 0);
    }


    public synchronized void showPopupMenu(final TreePath tp) {
        JPopupMenu popup = new JPopupMenu();
        JMenuItem item;

        if (tp == null) {
            return;
        } else {
            // Menu items for folders
            
            setSelectionPath(tp);
            final DefaultMutableTreeNode lastPathNode = (DefaultMutableTreeNode) tp.getLastPathComponent();
            final IntegerList children = new IntegerList();
            int actionTypes = CalculationAction.UNSUPPORTED_TYPE;
            
            if (lastPathNode.getAllowsChildren()) {
                item = new JMenuItem("Rename");
                item.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        startEditingAtPath(tp);
                    }
                });
                popup.add(item);

                item = new JMenuItem("New Class");
                item.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        synchronized (InteractiveTree.this) {
                            DefaultMutableTreeNode newClass = new DefaultMutableTreeNode("New Class", true); 
                            model.insertNodeInto(newClass, lastPathNode, 0);
                            startEditingAtPath(new TreePath(newClass.getPath()));
                        }
                        if (changeListener != null) changeListener.stateChanged(new ChangeEvent(this));
                    }
                });
                popup.add(item);
                
                item = new JMenuItem("Delete class");
                item.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        synchronized (InteractiveTree.this) {
                            TreeNode p = lastPathNode.getParent();
                            if (p == null) return;
                            if (!(p instanceof DefaultMutableTreeNode)) return;
                            DefaultMutableTreeNode parentNode = (DefaultMutableTreeNode) p;
                            deleteClass(model, lastPathNode, parentNode);
                            sortFolder(parentNode, model, false);
                            setSelectionPath(tp.getParentPath());
                        }
                    }                	
                });
                popup.add(item);
                
                item = new JMenuItem("Clear subclasses");
                item.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        synchronized (InteractiveTree.this) {
                            clearSubclasses(model, lastPathNode, lastPathNode);
                            sortFolder(lastPathNode, model, false);
                        }
                        if (changeListener != null) {
                            changeListener.stateChanged(new ChangeEvent(this));
                        }
                    }
                });
                popup.add(item);

//                item = new JMenuItem("Group Rest");
//                item.addActionListener(new ActionListener() {
//                    public void actionPerformed(ActionEvent e) {
//                        synchronized (InteractiveTree.this) {
//                            clearClassification(model, lastPathNode, lastPathNode);
//                            sortFolder(lastPathNode, model, false);
//                        }
//                        if (changeListener != null) {
//                            changeListener.stateChanged(new ChangeEvent(this));
//                        }
//                    }
//                });
//                popup.add(item);

                popup.add(new JSeparator());

                for (int i = 0; i < lastPathNode.getChildCount(); i++) {
                    Object userObject = ( (DefaultMutableTreeNode) lastPathNode.getChildAt(i)).getUserObject();
                    if (userObject instanceof Integer) children.add( (Integer) userObject);
                }
                actionTypes = CalculationAction.CLASS_ACTION;
                
            } else {
                // Menu items for cells
                
                // FIXME: Pretty sure this implementation runs inefficiently at O(N2).  Would be pretty easy to
                // have O(N), but doesn't really matter.
                item = new JMenuItem("Collapse Above");
                item.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        DefaultMutableTreeNode parent = (DefaultMutableTreeNode) lastPathNode.getParent();
                        int thisIndex = parent.getIndex(lastPathNode);
                        List<DefaultMutableTreeNode> cellNodes = getCellNodes(parent, 0, thisIndex);
                        if (cellNodes.size() < 2) return; // Not worth it then
                        
                        DefaultMutableTreeNode collapsed = new DefaultMutableTreeNode("collapsed", true);
                        
                        // Now move the collected cellNodes into the collapsed node and lastly put the 
                        // collapsed node into parent.
                        synchronized (this) {
                            for (DefaultMutableTreeNode cellNode : cellNodes) {
                                model.removeNodeFromParent(cellNode);
                                model.insertNodeInto(cellNode, collapsed, collapsed.getChildCount());
                            }
                            model.insertNodeInto(collapsed, parent, 0);
                        }
                        if (changeListener != null) changeListener.stateChanged(new ChangeEvent(this));
                    }
                });
                popup.add(item);
                
                // See "Collapse Above" for details
                item = new JMenuItem("Collapse Below");
                item.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        DefaultMutableTreeNode parent = (DefaultMutableTreeNode) lastPathNode.getParent();
                        int thisIndex = parent.getIndex(lastPathNode);                   	
                        List<DefaultMutableTreeNode> cellNodes = getCellNodes(parent, thisIndex+1, parent.getChildCount());
                        if (cellNodes.size() < 2) return; // Not worth it then
                        
                        DefaultMutableTreeNode collapsed = new DefaultMutableTreeNode("collapsed", true);
                        
                      synchronized (this) {
                            for (DefaultMutableTreeNode cellNode : cellNodes) {
                                model.removeNodeFromParent(cellNode);
                                model.insertNodeInto(cellNode, collapsed, collapsed.getChildCount());
                            }
                          model.insertNodeInto(collapsed, parent, parent.getChildCount());
                      }

                      if (changeListener != null) changeListener.stateChanged(new ChangeEvent(this));
                    }
                });
                popup.add(item);
                
                popup.add(new JSeparator());

                children.add((Integer) lastPathNode.getUserObject());
                actionTypes = CalculationAction.NEURON_ACTION;
            }
            
            // Add additional calculationActions
            // Ultimately these come from config.xml "Neuron Viewer Actions".
            // They are dynamically loaded at runtime via class introspection in NeuronViewer.
            for (int i = 0; i < calculationActions.size(); i++) {
                final CalculationAction a =
                    (CalculationAction) calculationActions.get(i);
                if (a.getAction() == actionTypes ||	a.getAction() == CalculationAction.CLASS_AND_NEURON_ACTION) {
                    item = new JMenuItem(a.getName());
                    item.addActionListener(new ActionListener() {
                        public void actionPerformed(ActionEvent e) {
                            a.doAction(children, tp);
                        }
                    });
                    popup.add(item);
                }
            }
        }

        int row = getRowForPath(tp);
        Rectangle r = getRowBounds(row);
        popup.show(this, r.x + r.width / 2, r.y + r.height / 2);
    }
    
    private static List<DefaultMutableTreeNode> getCellNodes(TreeNode parent, int start, int end) {
        ArrayList<DefaultMutableTreeNode> cellNodes = new ArrayList<DefaultMutableTreeNode>();
        for (int i = start; i < end; i++) {
            DefaultMutableTreeNode child = (DefaultMutableTreeNode) parent.getChildAt(i);
            if (child.getUserObject() instanceof Integer) cellNodes.add(child);
        }
        return cellNodes;
    }
    
    public void keyReleased(KeyEvent e) {
        if (e.getKeyCode() == KeyEvent.VK_F2)
            startEditingAtPath(getSelectionPath());
        if (e.getKeyCode() == KeyEvent.VK_2 && (e.getModifiersEx() & InputEvent.SHIFT_DOWN_MASK) == InputEvent.SHIFT_DOWN_MASK)
            startEditingAtPath(getSelectionPath());

        if (e.getKeyCode() == KeyEvent.VK_F3)
            showPopupMenu(getSelectionPath());
    }
    public void mouseReleased(MouseEvent e) {}
    public void mouseEntered(MouseEvent e) {}
    public void mouseExited(MouseEvent e) {}
    public void keyTyped(KeyEvent e) {}
    public void keyPressed(KeyEvent e) {}

}
