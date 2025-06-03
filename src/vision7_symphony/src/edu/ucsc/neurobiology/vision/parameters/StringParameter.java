package edu.ucsc.neurobiology.vision.parameters;

import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;


/**
 * This class represents a parameter that represents a String.
 * It can be shown on the screen using the ParametersTable class.
 * The value is displayed in a simple JTextField. When the string gets changed
 * a PropertyChangeEvent is send to all regitered listeners.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class StringParameter
    extends AbstractParameter {
    private JTextField textField;


    /**
     * Creates an instance of StringParameter which will display by default
     * the given string.
     *
     * @param name The internal name of the parameter
     * @param screenName The name which will be displayed on the screen
     *                   (can contain UNICODE characters)
     * @param toolTip The ToolTip associated with this parameter.
     * @param value The initial string displayed by this parameter
     */
    public StringParameter(String name, String screenName, String toolTip, String value) {
        super(name, screenName, toolTip);
        initialize(value);
    }


    /**
     * Retrives the text string of this parameter.
     * a legal double value.
     *
     * @return The double value or NaN
     */
    public String getValue() {
        return textField.getText();
    }


    protected void initialize(String initialValue) {
        super.initialize();

        textField = new JTextField(initialValue);
        textField.setToolTipText(toolTip);
        textField.getDocument().addDocumentListener(new DocumentListener() {
            public void insertUpdate(DocumentEvent e) {
                fireParameterChange();
            }


            public void removeUpdate(DocumentEvent e) {
                fireParameterChange();
            }


            public void changedUpdate(DocumentEvent e) {
                fireParameterChange();
            }
        });
        textField.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                fireParameterChange();
            }
        });

        textField.addFocusListener(new FocusAdapter() {
            public void focusLost(FocusEvent e) {
                fireParameterChange();
            }
        });
    }


    public JComponent getValueComponent() {
        return textField;
    }


    public String valueAsString() {
        return textField.getText();
    }


    public Object valueAsObject() {
        return textField.getText();
    }

}
