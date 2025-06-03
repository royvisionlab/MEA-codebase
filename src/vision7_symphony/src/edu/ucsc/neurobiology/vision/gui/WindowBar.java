package edu.ucsc.neurobiology.vision.gui;

import java.util.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.*;


/**
 *
 * Dumitru Petrusca, University of California, Santa Cruz
 */
public class WindowBar
    extends JPanel implements ActionListener {

    int oldSelectedID = -1;
    ArrayList<ActionListener> actionListeners = new ArrayList<ActionListener> ();


    public WindowBar() {
        super(new GridLayout(1, 0));

        setPreferredSize(new Dimension(0, 22));
    }


    public void addActionListener(ActionListener a) {
        if (a != null) {
            actionListeners.add(a);
        } else {
            throw new NullPointerException("null ActionListener provided");
        }
    }


    public void actionPerformed(ActionEvent e) {
        int selectedTabID = e.getID();

        if (e.getActionCommand().equals("button")) {
            if (oldSelectedID != -1 && selectedTabID == oldSelectedID) {
                return;
            } else {
                for (int i = 0; i < this.getComponentCount(); i++) {
                    Tab tab = (Tab)this.getComponent(i);
                    if (tab.tabID != selectedTabID) {
                        tab.setSelected(false);
                    }
                }

                oldSelectedID = selectedTabID;

                for (ActionListener list : actionListeners) {
                    list.actionPerformed(new ActionEvent(e.getSource(), selectedTabID,
                        e.getActionCommand()));
                }
            }
        } else if (e.getActionCommand().equals("check")) {
            for (ActionListener list : actionListeners) {
                list.actionPerformed(new ActionEvent(e.getSource(), selectedTabID,
                    e.getActionCommand()));
            }
        }
    }


    public int addButton(String text) {
        Tab tab = new Tab(text);
        tab.addActionListener(this);
        this.add(tab);
        return tab.tabID;
    }


    public void removeButton(int id) {
        for (int i = 0; i < this.getComponentCount(); i++) {
            Tab tab = (Tab)this.getComponent(i);
            if (tab.tabID == id) {
                this.remove(tab);
            }
        }
    }


    public void setSelection(int id) {
        for (int i = 0; i < this.getComponentCount(); i++) {
            Tab tab = (Tab)this.getComponent(i);
            tab.setSelected(tab.tabID == id);
        }
    }


    public void setChecked(int id, boolean checked) {
        for (int i = 0; i < this.getComponentCount(); i++) {
            Tab tab = (Tab)this.getComponent(i);
            if (tab.tabID == id) {
                tab.check.setSelected(checked);
            }
        }
    }

}


class Tab
    extends JPanel {

    private static int staticTabID = 0;
    JCheckBox check;
    JLabel button;
    Border selectedBorder = BorderFactory.createLoweredBevelBorder();
    Border deselectedBorder = BorderFactory.createRaisedBevelBorder();
    final int tabID;
    ArrayList<ActionListener> actionListeners = new ArrayList<ActionListener> ();


    public Tab(String text) {
        super(new BorderLayout());
        this.setBorder(deselectedBorder);

        tabID = staticTabID;
        staticTabID++;

        check = new JCheckBox("", true);
        check.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                for (ActionListener list : actionListeners) {
                    list.actionPerformed(new ActionEvent(check, tabID, "check"));
                }
            }
        });
        button = new JLabel(text);
        button.addMouseListener(new MouseAdapter() {
            public void mousePressed(MouseEvent e) {
                setSelected(true);
                for (ActionListener list : actionListeners) {
                    list.actionPerformed(new ActionEvent(button, tabID, "button"));
                }
            }
        });

        add(check, BorderLayout.WEST);
        add(button, BorderLayout.CENTER);
    }


    public void addActionListener(ActionListener a) {
        if (a != null) {
            actionListeners.add(a);
        } else {
            throw new NullPointerException("null ActionListener provided");
        }
    }


    public void setSelected(boolean selected) {
        if (selected) {
            this.setBorder(selectedBorder);
            button.setFont(new Font("Tahoma", Font.BOLD, 11));
        } else {
            this.setBorder(deselectedBorder);
            button.setFont(new Font("Tahoma", Font.PLAIN, 11));
        }
    }
}
