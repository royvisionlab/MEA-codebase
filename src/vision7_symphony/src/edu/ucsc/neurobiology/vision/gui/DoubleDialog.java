package edu.ucsc.neurobiology.vision.gui;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;


/**
 *
 * Dumitru Petrusca, University of California, Santa Cruz
 */
@SuppressWarnings("serial")
public class DoubleDialog
    extends JDialog {

    GridBagLayout gridBagLayout1 = new GridBagLayout();
    JButton okButton = new JButton();
    JLabel messageLabel = new JLabel();
    JTextField inputTextField = new JTextField();
    String message;
    double defaultValue, minValue, maxValue;


    private DoubleDialog(Component parent, String title, String message,
                         double defaultValue, double minValue, double maxValue) {

        super( (JFrame)null, title, true);
        this.message = message;
        this.defaultValue = defaultValue;
        this.minValue = minValue;
        this.maxValue = maxValue;
        try {
        jbInit();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    public static double showDoubleInputDialog(
        Component parent, String title, String message, double defaultValue,
        double minValue, double maxValue) {

        DoubleDialog d = new DoubleDialog(parent, title, message, defaultValue, minValue,
                                          maxValue);
        d.setBounds(400, 150, 500, 200);
        d.setVisible(true);
        return Double.parseDouble(d.inputTextField.getText());
    }


    private void jbInit() throws Exception {
        this.getContentPane().setLayout(gridBagLayout1);
        okButton.setMnemonic('O');
        okButton.setText("Ok");
        okButton.addActionListener(new DoubleDialog_okButton_actionAdapter(this));
        messageLabel.setText(message);
        inputTextField.setToolTipText("Enter The Value Here");
        inputTextField.setSelectionStart(11);
        inputTextField.setText("" + defaultValue);
        this.setDefaultCloseOperation(javax.swing.WindowConstants.DO_NOTHING_ON_CLOSE);
        this.setForeground(Color.black);
        this.setModal(true);


        this.add(okButton, new GridBagConstraints(0, 2, 1, 1, 0.5, 0.5
                                                  , GridBagConstraints.SOUTH,
                                                  GridBagConstraints.NONE,
                                                  new Insets(5, 5, 5, 5),
                                                  0, 0));
        this.add(messageLabel, new GridBagConstraints(0, 0, 1, 1, 1.0, 0.5
            , GridBagConstraints.SOUTH, GridBagConstraints.HORIZONTAL,
            new Insets(0, 10, 0, 10), 0, 0));
        this.add(inputTextField, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            , GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL,
            new Insets(10, 10, 10, 10), 0, 0));
        inputTextField.getDocument().addDocumentListener(new DocumentListener() {
            public void insertUpdate(DocumentEvent e) {
                updateState();
            }


            public void removeUpdate(DocumentEvent e) {
                updateState();
            }


            public void changedUpdate(DocumentEvent e) {
                updateState();
            }


            private void updateState() {
                try {
                    double v = Double.parseDouble(inputTextField.getText());
                    if (v < minValue || v > maxValue) {
                        okButton.setEnabled(false);
                    } else {
                        okButton.setEnabled(true);
                    }
                } catch (NumberFormatException ex) {
                    okButton.setEnabled(false);
                }
            }
        });
    }


    void okButton_actionPerformed(ActionEvent e) {
        this.setVisible(false);
        this.dispose();
    }
    
    public static void main(String args[]) {

        DoubleDialog.showDoubleInputDialog(null, "Threshold", "Input the Threshold", 4, 1, 1e10);
    }

}


class DoubleDialog_okButton_actionAdapter
    implements java.awt.event.ActionListener {
    DoubleDialog adaptee;

    DoubleDialog_okButton_actionAdapter(DoubleDialog adaptee) {
        this.adaptee = adaptee;
    }


    public void actionPerformed(ActionEvent e) {
        adaptee.okButton_actionPerformed(e);
    }
    

}
