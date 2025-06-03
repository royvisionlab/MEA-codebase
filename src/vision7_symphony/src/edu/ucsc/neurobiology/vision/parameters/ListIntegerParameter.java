package edu.ucsc.neurobiology.vision.parameters;

import java.awt.event.*;
import javax.swing.*;


/**
 * This class implements the <code>Parameter</code> interface to represent a
 * list of Java integer values in a selectable and editable ComboBox. It also
 * provides getter and setter funtions for the value of the parameter. Every time
 * the value changes a <code>PropetryChangeEvent</code> is dispatched to all listeners.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ListIntegerParameter
    extends AbstractParameter {
    private JComboBox comboBox;
    private int value;


    /**
     * Creates an instance of ListIntegerParameter which will display by default
     * the given value. The needed possible choises have to be added later.
     *
     * @param name the internal name of the parameter
     * @param screenName the screen name associated with the parameter
     * @param toolTip the explanatory tooltop associated with the parameter
     * @param value the default value represented by the parameter
     */
    public ListIntegerParameter(String name, String screenName, String toolTip, int value) {
        super(name, screenName, toolTip);
        this.value = value;
        initialize();
    }


    public void initialize() {
        super.initialize();

        comboBox = new JComboBox();
        comboBox.setEditable(true);
        comboBox.setToolTipText(toolTip);
        /*
                comboBox.addgetDocument().addDocumentListener(new DocumentListener() {
                    public void insertUpdate(DocumentEvent e) { valueChanged(); }
                    public void removeUpdate(DocumentEvent e) { valueChanged(); }
                    public void changedUpdate(DocumentEvent e) { valueChanged(); }
                });
         */
        comboBox.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                if (e.getStateChange() == ItemEvent.SELECTED) {
                    valueChanged();
                    fireParameterChange();
                }
            }
        });
        /*
                 comboBox.addFocusListener(new FocusAdapter() {
            public void focusLost(FocusEvent e) {
                valueChanged();
            }
                 });
         */
    }


    /**
     * Adds the given integer to the possible selections of this parameter.
     * The UI updates respectively.
     */
    public void addInteger(int value) {
        comboBox.addItem("" + value);
    }


    /**
     * Returns the integer value currently displaying in the parameter.
     */
    public int getValue() {
        return value;
    }


    /**
     * Sets the integer value currently selected in the parameter.
     */
    public void setValue(int value) {
        comboBox.setSelectedItem("" + value);
    }


    private void valueChanged() {
        try {
            value = Integer.parseInt( (String) comboBox.getSelectedItem());
            fireParameterChange();
        } catch (NumberFormatException e) {
        }
    }


    public JComponent getValueComponent() {
        return comboBox;
    }


    public String valueAsString() {
        return Integer.toString(value);
    }


    public Object valueAsObject() {
        return value;
    }
}
