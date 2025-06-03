package edu.ucsc.neurobiology.vision.gui;

import java.beans.*;
import java.util.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;


/**
 *
 * Dumitru Petrusca, University of California, Santa Cruz
 */
public class Desktop
    extends JPanel implements ActionListener {

    /**
     * An interface to declare a selection contract.
     *
     * @author Dumitru Petrusca, University of California, Santa Cruz
     */
    public static interface SelectionListener {
        public void selectionHappened(JInternalFrame frame);
    }


    JDesktopPane mainDesktop;
    private GridDesktopManager manager;
    WindowBar buttonPanel = new WindowBar();
    HashMap<Integer, JInternalFrame> id2frame = new HashMap<Integer, JInternalFrame> ();
    HashMap<JInternalFrame, Integer> frame2id = new HashMap<JInternalFrame, Integer> ();
    ArrayList<SelectionListener> selectionListeners = new ArrayList<SelectionListener> ();


    public Desktop() {
        super(new BorderLayout());

        mainDesktop = new JDesktopPane();
        manager = new GridDesktopManager(this);
        mainDesktop.setDesktopManager(manager);
        mainDesktop.addContainerListener(manager);
        mainDesktop.addComponentListener(manager);

        add(mainDesktop, BorderLayout.CENTER);
        add(buttonPanel, BorderLayout.SOUTH);

        buttonPanel.addActionListener(this);
    }


    public void actionPerformed(ActionEvent e) {
        int id = e.getID();
        JInternalFrame frame = id2frame.get(id);
        if (e.getActionCommand().equals("button")) {
            try {
                frame.setSelected(true);
            } catch (PropertyVetoException ex) {
                ex.printStackTrace();
            }
        } else if (e.getActionCommand().equals("check")) {
            try {
//                frame.setVisible( ( (JCheckBox) e.getSource()).isSelected());
                frame.setIcon(! ( (JCheckBox) e.getSource()).isSelected());
            } catch (PropertyVetoException ex) {
                ex.printStackTrace();
            }
        }
        manager.updateLayout();
    }


    public void addSelectionListener(SelectionListener a) {
        if (a != null) {
            selectionListeners.add(a);
        } else {
            throw new NullPointerException("null SelectionListener provided");
        }
    }


    public void setCurrentFrame(JInternalFrame frame) {
        if (frame != null) {
            Integer id = frame2id.get(frame);
            buttonPanel.setSelection(id.intValue());
        }

        for (SelectionListener list : selectionListeners) {
            list.selectionHappened(frame);
        }
    }


    public void updateLayout() {
        manager.updateLayout();
    }


    public JMenu getCurrentViewsMenu() {
        return manager.getCurrentViewsMenu();
    }


    public JInternalFrame getSelectedFrame() {
        return mainDesktop.getSelectedFrame();
    }


    public JInternalFrame[] getAllFrames() {
        return mainDesktop.getAllFrames();
    }


    public void addFrame(final JInternalFrame frame, String title) {
        int id = buttonPanel.addButton(title);
        id2frame.put(id, frame);
        frame2id.put(frame, id);

        frame.addInternalFrameListener(new InternalFrameAdapter() {
            public void internalFrameIconified(InternalFrameEvent e) {
                int id = frame2id.get(frame);
                buttonPanel.setChecked(id, false);
            }


            public void internalFrameDeiconified(InternalFrameEvent e) {
                int id = frame2id.get(frame);
                buttonPanel.setChecked(id, true);
            }
        });

        frame.setTitle(title);
        mainDesktop.add(frame);
        frame.setVisible(true);
    }


    public void removeFrame(final JInternalFrame frame) {
        Object id = frame2id.remove(frame);
        id2frame.remove(id);

        // FIXME
        if (id != null) {
            buttonPanel.removeButton( (Integer) id);
        }

        mainDesktop.remove(frame);
        mainDesktop.setSelectedFrame(null);
    }
}
