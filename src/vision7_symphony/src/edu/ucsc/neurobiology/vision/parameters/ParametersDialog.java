package edu.ucsc.neurobiology.vision.parameters;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

import edu.ucsc.neurobiology.vision.Vision;


/**
 * This class represents a two column table, each row containing one parameter.
 * The left column contains labels with parameters screen names, the right one
 * contains the component used to display and change the value of the parameter.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ParametersDialog
    extends JDialog implements ChangeListener, ActionListener {

    public enum DialogSelection {
        CANCEL, OK, SAVE;
    }


    private DialogSelection result = DialogSelection.CANCEL;
    private JButton cancelButton, okButton, saveButton;

    public void actionPerformed(ActionEvent e) {
        if (e.getSource().equals(okButton)) {
            result = DialogSelection.OK;
        } else if (e.getSource().equals(saveButton)) {
            result = DialogSelection.SAVE;
        } else if (e.getSource().equals(cancelButton)) {
            result = DialogSelection.CANCEL;
        } else {
            throw new Error("Unknown source of event");
        }

        setVisible(false);
    }


    /**
     * Creates a ParametersDialog which is initially empty.
     */
    public ParametersDialog(Component owner, String title, ParametersTable table) {
        super( (JFrame) SwingUtilities.getAncestorOfClass(JFrame.class, owner), title, true);
//        this.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);

        JScrollPane pane = new JScrollPane(table,
                                           JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
                                           JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

        pane.setBorder(null);
        getContentPane().add(pane, BorderLayout.CENTER);

        okButton = new JButton("Ok/Run");
        okButton.setMnemonic(KeyEvent.VK_O);
        okButton.addActionListener(this);
        saveButton = new JButton("Save");
        saveButton.setMnemonic(KeyEvent.VK_C);
        saveButton.addActionListener(this);
        cancelButton = new JButton("Cancel");
        cancelButton.setMnemonic(KeyEvent.VK_C);
        cancelButton.addActionListener(this);

        JPanel p = new JPanel();
        p.add(okButton);
        p.add(saveButton);
        p.add(cancelButton);
        this.add(p, BorderLayout.SOUTH);
        
        pack(); // To get native header size accounted for
        setPreferredSize(new Dimension(calcPreferredWidth(), getHeight()));
        pack(); // Repack to get JDialog and ContentPane to nice width
        table.percentFillColumns(new int[]{40, 60}); // Now adjust table columns to percents of ContentPane
        setLocationRelativeTo(owner);
        setAlwaysOnTop(true);

//        table.addStateChangeListener(this);
    }

    
    private int calcPreferredWidth() {
        int mainFrameWidth  = Vision.getInstance().getMainFrame().getSize().width;
        return mainFrameWidth * 60 / 100;
    }

    public boolean isOkSelected() {
        return (result == DialogSelection.OK);
    }


    public boolean isCancelSelected() {
        return (result == DialogSelection.CANCEL);
    }


    public boolean isSaveSelected() {
        return (result == DialogSelection.SAVE);
    }


//    public DialogSelection getResult() {
//        return result;
//    }


    /**
     * This method is public as an implementations side effect, do NOT call.
     */
    public void stateChanged(ChangeEvent e) {
//        if (table.isValid()) {
//            okButton.setEnabled(true);
//        } else {
//            okButton.setEnabled(false);
//        }
    }


    public void dispose() {
        this.getContentPane().removeAll();

        okButton.removeActionListener(this);
        okButton = null;
        saveButton.removeActionListener(this);
        saveButton = null;
        cancelButton.removeActionListener(this);
        cancelButton = null;

        super.dispose();
        try {
            super.finalize();
        } catch (Throwable ex) {
        }
    }
}
