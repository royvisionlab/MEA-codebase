package edu.ucsc.neurobiology.vision.parameters;

import java.util.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.table.*;
import javax.swing.tree.*;

import org.w3c.dom.*;

import edu.ucsc.neurobiology.vision.swing.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * An implementation of the AbstractTreeTableModel that allows the display of a list
 * of parameters as a tree-table. This model is used in many places in Vision where there
 * is the need to display a list of parameters.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ParameterModel
    extends AbstractTreeTableModel {

    private int counter = -1;

    // Names of the columns.
    static protected String[] cNames = {"Name", "Value"};

    // Types of the columns.
    static protected Class[] cTypes = {TreeTableModel.class, String.class};

    JTree tree;

    /**
     * Creates a FileSystemModel2 with the root being <code>rootPath</code>.
     * This does not load it, you should invoke
     * <code>reloadChildren</code> with the root to start loading.
     */
    public ParameterModel(Node node) {
        super(null);
        root = new TreeTableNode(null, node, "");
    }


    public ParameterModel() {
        super(null);
        root = new TreeTableNode(null, new BooleanParameter("root", "root", "", true), "");
        ( (TreeTableNode) root).noExpand = true;
    }


    public String getTooltipOf(Object node) {
        return "Tooltip " + node;
    }


    public int getParametersCount() {
        return counter;
    }


    public void setTree(JTreeTable treetable) {
        tree = treetable.getTree();
    }


    //
    // The TreeModel interface
    //

    /**
     * Returns the number of children of <code>node</code>.
     */
    public int getChildCount(Object node) {
//        System.err.println(node +" " + ( (TreeTableNode) node).getChildCount());
        return ( (TreeTableNode) node).getChildCount();
    }


    /**
     * Returns the child of <code>node</code> at index <code>i</code>.
     */
    public Object getChild(Object node, int i) {
        return ( (TreeTableNode) node).getChild(i);
    }


    /**
     * Returns true if the passed in object represents a leaf, false
     * otherwise.
     */
    public boolean isLeaf(Object node) {
        return ( (TreeTableNode) node).isLeaf();
    }


    //
    //  The TreeTableNode interface.
    //

    /**
     * Returns the number of columns.
     */
    public int getColumnCount() {
        return cNames.length;
    }


    /**
     * Returns the name for a particular column.
     */
    public String getColumnName(int column) {
        return cNames[column];
    }


    /**
     * Returns the class for the particular column.
     */
    public Class getColumnClass(int column) {
        return cTypes[column];
    }


    /**
     * Returns the value of the particular column.
     */
    public Object getValueAt(Object node, int column) {
//        TreeTableNode fn = (TreeTableNode) node;

//        switch (column) {
//            case 0:
//                return "";
//            case 1:
//                return "";
//        }

        return null;
    }


    public void addParameter(Parameter p) {
        TreeTableNode parent = (TreeTableNode) root;
        TreeTableNode currentNode = (TreeTableNode) root;
        String[] names = StringUtil.decomposeString(p.getName(), ".");

        for (int i = 0; i < names.length - 1; i++) {
            parent = currentNode;
            currentNode = currentNode.getChild(names[i]);

            if (currentNode == null) {
                currentNode = new TreeTableNode(
                    parent, new BooleanParameter(names[i], "", "", true), names[i]);
                parent.addChild(currentNode);
            }

//            System.err.println("currentNode = " + currentNode + " : " + currentNode.parent);
        }

        TreeTableNode child = new TreeTableNode(currentNode, p, names[names.length - 1]);
        currentNode.addChild(child);
//        System.err.println("currentNode = " + child + " : " + child.parent);

        TreeTableNode r = (TreeTableNode) root;
//        for (int i = 0; i < r.getChildCount(); i++) {
//            System.err.println("child " + r.getChild(i));
//        }

        fireTreeStructureChanged(this, r.getPath(), null, null);
    }


    public Parameter getParameter(String name) {
        TreeTableNode currentNode = (TreeTableNode) root;
//        System.err.println("currentNode = " + currentNode);
        String[] names = StringUtil.decomposeString(name, ".");
        for (int i = 0; i < names.length; i++) {
            currentNode = currentNode.getChild(names[i]);
//            System.err.println("currentNode = " + currentNode);
            if (currentNode == null) {
                return null;
            }
        }

        return currentNode.p;
    }


    public Parameter getParameter(final int index) {
        return getParameter(index, (TreeTableNode) root);
    }


    private Parameter getParameter(final int index, TreeTableNode node) {
//        System.err.println("checking " + node.name + ", " + node.getChildCount() + " children");

        for (int i = 0; i < node.getChildCount(); i++) {
            TreeTableNode child = node.getChild(i);
            if (child.index == index) {
                return child.p;
            }

            if (!child.isLeaf()) {
//                System.err.println("going into " + child.name);
                Parameter p = getParameter(index, child);
                if (p != null) {
                    return p;
                }
            }
        }

        return null;
    }
    
    
    public class TreeTableNode {
        private final int index;
        private TreeTableNode parent;
        private ArrayList<TreeTableNode> children = new ArrayList<TreeTableNode> ();
        private boolean noExpand = false;
        public Parameter p;
        public String shortName;

        public void addChild(TreeTableNode c) {
            children.add(c);
            c.parent = this;
        }


        protected TreeTableNode(TreeTableNode parent, Parameter p, String shortName) {
            counter++;
            index = counter - 1; // the root is index -1
            this.parent = parent;
            this.p = p;
            this.shortName = shortName;
        }


        protected TreeTableNode(TreeTableNode parent, Node node, String longName) {
            counter++;
            index = counter - 1; // the root is index -1
            this.parent = parent;

            NamedNodeMap attributes = node.getAttributes();
            String name = attributes.getNamedItem("name").getNodeValue().trim();
            shortName = name;
            String value, screenName;
            if (attributes.getNamedItem("value") != null) {
                value = attributes.getNamedItem("value").getNodeValue().trim();
            } else {
                value = "true";
                noExpand = true;
            }
            Node screenNameNode = attributes.getNamedItem("screenName");
            if (screenNameNode != null) {
                screenName = screenNameNode.getNodeValue();
            } else {
                screenName = name;
            }

//            System.err.println(screenName);

            String toolTip = screenName;

            if (!node.getNodeName().equals("ParametersGroup")) {
                String type = node.getNodeName().trim();
                Node tooltipNode = attributes.getNamedItem("toolTip");
                if (tooltipNode != null) {
                    toolTip = tooltipNode.getNodeValue();
                }

                if (type.equals("FileParameter")) {
                    String ext = attributes.getNamedItem("extension").getNodeValue().
                                 trim();
                    p = new FileParameter(longName + name, screenName, toolTip, value,
                                          ext);
                } else if (type.equals("FileListParameter")) {
                    String ext = attributes.getNamedItem("extension").getNodeValue().
                                 trim();
                    p = new FileListParameter(longName + name, screenName, toolTip,
                                              value, ext);
                } else if (type.equals("IntegerParameter")) {
                    p = new IntegerParameter(
                        longName + name, screenName, toolTip, Integer.parseInt(value),
                        Integer.MIN_VALUE, Integer.MAX_VALUE);
                } else if (type.equals("EnumeratorParameter")) {
                    String values = attributes.getNamedItem("values").getNodeValue();
                    EnumeratorParameter ep = new EnumeratorParameter(longName + name,
                        screenName, toolTip);
                    StringTokenizer tokenizer = new StringTokenizer(values, ":", false);
                    while (tokenizer.hasMoreTokens()) {
                        String text = tokenizer.nextToken().trim();
                        int number = Integer.parseInt(tokenizer.nextToken().trim());
                        ep.addChoice(number, text);
                    }
                    ep.setValue(Double.parseDouble(value));
                    p = ep;
                } else if (type.equals("BooleanParameter")) {
                    p = new BooleanParameter(longName + name, screenName, toolTip,
                                             Boolean.valueOf(value).booleanValue());
                } else if (type.equals("DoubleParameter")) {
                    p = new DoubleParameter(longName + name, screenName, toolTip,
                                            Double.parseDouble(value));
                } else if (type.equals("StringParameter")) {
                    p = new StringParameter(longName + name, screenName, toolTip, value);
                } else if (type.equals("ColorParameter")) {
                    p = new ColorParameter(longName + name, screenName, toolTip,
                                           Color.black);
                } else {
                    throw new Error("Unknown parameter type: " + type);
                }
            } else {
                p = new BooleanParameter(longName + name, screenName, toolTip,
                                         Boolean.valueOf(value).booleanValue());
                ( (JCheckBox) p.getValueComponent()).addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        fireTreeStructureChanged(this, getPath(), null, null);
                        if ( ( (BooleanParameter) p).isSelected()) {
                            tree.expandPath(new TreePath(getPath()));
                        }
                    }
                });
            }

            // load children
            NodeList paramNodeList = node.getChildNodes();
            for (int i = 0; i < paramNodeList.getLength(); i++) {
                Node paramNode = paramNodeList.item(i);
                if (paramNode.getNodeType() == Node.ELEMENT_NODE) {
                    TreeTableNode n = new TreeTableNode(this, paramNode,
                        index == -1 ? "" : longName + name + ".");
                    children.add(n);
                }
            }

//            System.err.println(p.getName());
        }


        public String getName() {
            return p.getName();
        }


        public TreeTableNode getChild(String name) {
            for (TreeTableNode node : children) {
                if (node.shortName.equals(name)) {
                    return node;
                }
            }
            return null;
        }


        /**
         * Returns the the string to be used to display this leaf in the JTree.
         */
        public String toString() {
            return p.getScreenName();
        }
        

        /**
         * Returns the parent of the receiver.
         */
        public TreeTableNode getParent() {
            return parent;
        }


        /**
         * Returns true if the receiver represents a leaf, that is it is
         * isn't a directory.
         */
        public boolean isLeaf() {
            return children.size() == 0;
        }


        public int getChildCount() {
            if (p instanceof BooleanParameter && ! ( (BooleanParameter) p).isSelected()) {
//                System.err.println("sp[");
                return 0;
            } else {
//                System.err.println("st");
                return children.size();
            }
        }


        public TreeTableNode getChild(int i) {
            return children.get(i);
        }


        /**
         * Gets the path from the root to the receiver.
         */
        public TreeTableNode[] getPath() {
            return getPathToRoot(this, 0);
        }


        protected TreeTableNode[] getPathToRoot(TreeTableNode aNode, int depth) {
            TreeTableNode[] retNodes;

            if (aNode == null) {
                if (depth == 0) {
                    return null;
                } else {
                    retNodes = new TreeTableNode[depth];
                }
            } else {
                depth++;
                retNodes = getPathToRoot(aNode.getParent(), depth);
                retNodes[retNodes.length - depth] = aNode;
            }
            return retNodes;
        }
    }


    public void setValueAt(Object aValue, Object node, int column) {}


    public TableCellRenderer getRenderer() {
        return new CR();
    }


    public TableCellEditor getEditor() {
        return new CR();
    }


    private class CR
        extends AbstractCellEditor implements TableCellRenderer, TableCellEditor {

        TreeTableNode editingNode;

        public CR() {}


        public Component getTableCellRendererComponent(JTable table,
            Object value, boolean isSelected, boolean hasFocus, int row, int column) {

            TreePath path = tree.getPathForRow(row);
            TreeTableNode node = (TreeTableNode) path.getLastPathComponent();
            if (node.noExpand) {
                return new JLabel();
            } else {
                JComponent c = node.p.getValueComponent();
                //This block causes conflicts with extensions of JSpinner.  Not sure if it is needed. -- mgrivch 2007-08-22
//                if (isSelected) {
//                    if (c instanceof JSpinner) {
//                        ( (JComponent) c.getComponent(2)).getComponent(0).setBackground(
//                            UIManager.getColor("Table.selectionBackground"));
//                        ( (JComponent) c.getComponent(2)).getComponent(0).setForeground(
//                            UIManager.getColor("Table.selectionForeground"));
//                    } else {
//                        c.setBackground(UIManager.getColor("Table.selectionBackground"));
//                        c.setForeground(UIManager.getColor("Table.selectionForeground"));
//                    }
//                } else {
//                    if (c instanceof JSpinner) {
//                        ( (JComponent) c.getComponent(2)).getComponent(0).setBackground(
//                            Color.white);
//                        ( (JComponent) c.getComponent(2)).getComponent(0).setForeground(
//                            Color.black);
//                    }
//
//                    c.setBackground(Color.white);
//                    c.setForeground(Color.black);
//                }

                return c;
            }
        }


        public Component getTableCellEditorComponent(
            JTable table, Object value, boolean isSelected, int row, int column) {

            TreePath path = tree.getPathForRow(row);
            editingNode = (TreeTableNode) path.getLastPathComponent();
            if (editingNode.noExpand) {
                return new JLabel();
            } else {
                JComponent c = editingNode.p.getValueComponent();

                if (c instanceof JSpinner) {
                    ( (JSpinner) c).getEditor().setForeground(UIManager.getColor(
                        "Table.selectionForeground"));
                    ( (JSpinner) c).getEditor().setBackground(UIManager.getColor(
                        "Table.selectionForeground"));
                } else {
                    c.setForeground(UIManager.getColor("Table.selectionForeground"));
                    c.setBackground(UIManager.getColor("Table.selectionBackground"));
                }

                return c;
            }
        }


        public Object getCellEditorValue() {
            return editingNode.p.valueAsString();
        }
    }

}
