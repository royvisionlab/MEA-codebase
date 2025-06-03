package edu.ucsc.neurobiology.vision.parameters;

import java.awt.Dimension;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.ArrayList;
import java.util.Enumeration;

import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.UIDefaults;
import javax.swing.UIManager;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.table.DefaultTableColumnModel;
import javax.swing.table.TableColumn;

import org.w3c.dom.Node;

import edu.ucsc.neurobiology.vision.swing.JTreeTable;


/**
 * This class represents a two column table, each row containing one parameter.
 * The left column contains labels with parameters screen names, the right one
 * contains the component used to display and change the value of the parameter.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ParametersTable
    extends JTreeTable implements PropertyChangeListener {

    private ArrayList<ChangeListener> changeListeners;

    public ParametersTable() {
        super(new ParameterModel());
        ( (ParameterModel) treeTableModel).setTree(this);
        setRowHeight(18);
        getColumn("Value").setCellRenderer( ((ParameterModel) treeTableModel).getRenderer());
        getColumn("Value").setCellEditor(   ((ParameterModel) treeTableModel).getEditor());
        getTree().setEditable(true);

        changeListeners = new ArrayList<ChangeListener>();
    }


    /**
     * Creates a ParametersTable which is initially empty.
     */
    public ParametersTable(Node groupNode) {
        super(new ParameterModel(groupNode));
        ( (ParameterModel) treeTableModel).setTree(this);

        getColumn("Value").setCellRenderer( ( (ParameterModel) treeTableModel).
                                           getRenderer());
        getColumn("Value").setCellEditor( ( (ParameterModel) treeTableModel).getEditor());
        setRowHeight(18);
        getTree().setEditable(true);

        changeListeners = new ArrayList<ChangeListener>();
    }


    public void addStateChangeListener(ChangeListener listener) {
        changeListeners.add(listener);
    }


    /**
     * Adds a new <code>Parameter</code> to this table. The parameter will appear
     * immediately in the UI.
     *
     * @param p the parameter to be added
     */
    public void addParameter(Parameter p) {
        ( (ParameterModel) treeTableModel).addParameter(p);
        expandAll();
    }


    /**
     * Adds a new <code>Parameter</code> to this table. The <code>PropertyChangeListener
     * </code> listener will be added to this parameter before the paremeter is added
     * to the table. The parameter will appear immediately in the UI.
     *
     * @param p the parameter to be added
     * @param listener the listener to be added to this parameter
     */
    public void addParameter(Parameter p, PropertyChangeListener listener) {
        p.addPropertyChangeListener(listener);
        p.addPropertyChangeListener(this);
        addParameter(p);
    }


    /**
     * Finds and returns the parameter which has the specified internal name.
     * Returns <b>null</b> if the parameter does not exist.
     */
    public Parameter getParameter(String name) {
        return ( (ParameterModel) treeTableModel).getParameter(name);
    }


    public int getIntParameter(String name) {
        return ( (IntegerParameter) getParameter(name)).getValue();
    }


    public String getStringParameter(String name) {
        return ( (StringParameter) getParameter(name)).getValue();
    }


    public double getDoubleParameter(String name) {
        return ( (DoubleParameter) getParameter(name)).getValue();
    }

    public double getEnumeratorParameter(String name) {
        return ((EnumeratorParameter) getParameter(name)).getValue();
    }

    public boolean getBooleanParameter(String name) {
        return ( (BooleanParameter) getParameter(name)).getValue();
    }


    public String getFileParameter(String name) {
        return ( (FileParameter) getParameter(name)).getValue();
    }


    public Parameter getParameter(final int index) {
        return ( (ParameterModel) treeTableModel).getParameter(index);
    }


    /**
     * This method is public as an implementations side effect, do NOT call.
     */
    public void propertyChange(PropertyChangeEvent e) {
        this.repaint();

        ChangeEvent event = new ChangeEvent(this);
        for (int i = 0; i < changeListeners.size(); i++) {
            ChangeListener l = (ChangeListener) changeListeners.get(i);
            l.stateChanged(event);
        }
    }


    public boolean isValid() {
        return true;
    }


    public int getParametersCount() {
        return ( (ParameterModel) treeTableModel).getParametersCount();
    }


    public void expandAll() {
        for (int i = 0; i < this.getRowCount(); i++) this.getTree().expandRow(i);
    }

    
    public Dimension getPreferredScrollableViewportSize() {
        return getPreferredSize();
    }


    public void percentFillColumns(int[] percents) {
        DefaultTableColumnModel colModel = (DefaultTableColumnModel) getColumnModel();
        int totalWidth = getParent().getWidth();
        for (int i = 0; i < percents.length; i++) {
            TableColumn col = colModel.getColumn(i);
            col.setPreferredWidth(totalWidth * percents[i] / 100);
        }
        doLayout();
    }
  
    
    public static void main(String[] args) throws Exception {
        UIDefaults defaults = UIManager.getDefaults();
        defaults.put("ScrollPane.background", defaults.get("Tree.textBackground"));
        defaults.put("Viewport.background", defaults.get("Tree.background"));

        Enumeration<Object> e = defaults.keys();
        while (e.hasMoreElements()) {
            Object o = e.nextElement();
            System.err.println(o + " : " + defaults.get(o));
        }

//        Config config = new Config("config.xml");
//        ParametersTable t = new ParametersTable(config.getGroupNode("Spike Finding"));

        ParametersTable t = new ParametersTable();
        t.addParameter(new IntegerParameter("int1", "int1", "", 12, 0, 100));
        t.addParameter(new IntegerParameter("int2", "int2", "", 12, 0, 100));
        t.addParameter(new IntegerParameter("int3", "int3", "", 12, 0, 100));

        JFrame f = new JFrame();
        JScrollPane sp = new JScrollPane(t);
//        sp.setBackground(Color.red);
//        sp.setForeground(Color.red);
//        sp.getViewport().setBackground(Color.red);
        f.add(sp);
        f.setBounds(100, 100, 400, 400);
        f.setVisible(true);

        System.err.println(t.getPreferredSize());
        System.err.println(sp.getPreferredSize());
    }

}
