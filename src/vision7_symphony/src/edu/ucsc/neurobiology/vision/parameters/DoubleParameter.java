package edu.ucsc.neurobiology.vision.parameters;

import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;


/**
 * This class represents a parameter that represents a primitive double value.
 * It can be shown on the screen using the ParametersTable class.
 * The value is displayed in a simple JTextField. When the value gets changed
 * a PropertyChangeEvent is send to all regitered listeners.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class DoubleParameter
    extends AbstractParameter {

    private JTextField textField;
    private double value;
    private double oldValue;


    /**
     * Creates an instance of IntegerParameter which will display by default
     * the given double value.
     *
     * @param name The internal name of the parameter
     * @param screenName The name which will be displayed on the screen
     *                   (can contain UNICODE characters)
     * @param toolTip The ToolTip associated with this parameter.
     * @param value The initial value of this parameter
     */
    public DoubleParameter(String name, String screenName, String toolTip, double value) {
        super(name, screenName, toolTip);
        this.value = value;
        this.oldValue = value;
        initialize();
    }


    /**
     * Changes the internal double value of the parameteer and updaes the UI
     * representation. A PropertyChangeEvent is send to all regitered listeners.
     *
     * @param newValue The new value of the parameter
     */
    public void setValue(double newValue) {
        textField.setText("" + newValue);
        textField.repaint();
    }


    /**
     * Retrives the internal value of the parameter by parsing the text typed into
     * the text field by the user. Returns Double.NaN if the text does not represent
     * a legal double value.
     *
     * @return The double value or NaN
     */
    public double getValue() {
        return value;
    }


    protected void initialize() {
        super.initialize();

        textField = new JTextField();
        textField.setToolTipText(toolTip);
        textField.getDocument().addDocumentListener(new DocumentListener() {
            public void insertUpdate(DocumentEvent e) {
                applyInput();
            }


            public void removeUpdate(DocumentEvent e) {
                applyInput();
            }


            public void changedUpdate(DocumentEvent e) {
                applyInput();
            }
        });

        textField.addFocusListener(new FocusAdapter() {
            public void focusLost(FocusEvent e) {
                applyInput();
            }
        });

        setValue(value);
    }


    private void applyInput() {
        if (isValid()) {
            textField.setForeground(normalForeground);
            textField.setBackground(normalBackground);

            value = Double.parseDouble(textField.getText());
            if (value != oldValue) {
                oldValue = value;
                fireParameterChange();
            }
        } else {
            textField.setForeground(errorForeground);
            textField.setBackground(errorBackground);
        }
    }


    /**
     * Returns true if the specified number is acceptable for this parameter
     * and false otherwise.
     */
    public boolean isValid() {
        try {
            Double.parseDouble(textField.getText());
            return true;
        } catch (NumberFormatException nfe) {
            return false;
        }
    }


    public JComponent getValueComponent() {
        return textField;
    }


    public String valueAsString() {
        return Double.toString(value);
    }


    public Object valueAsObject() {
        return value;
    }

}
