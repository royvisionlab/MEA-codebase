package edu.ucsc.neurobiology.vision.parameters;

import java.util.*;

import java.awt.event.*;
import javax.swing.*;


/**
 * This class implements the <code>Parameter</code> interface to represent a
 * list of String to double mappings available in a ComboBox. It also provides
 * getter and setter functions for the value of the parameter. Every time the value
 * changes a <code>PropetryChangeEvent</code> is dispatched to all listeners.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class EnumeratorParameter
    extends AbstractParameter {

    private JComboBox comboBox;
    private ArrayList<String> texts;
    private ArrayList<Double> values;


    /**
     * Creates an instance of EnumeratorParameter which will be empty by default.
     * The needed String-Integer associations have to be added later using the
     * <code>addChoice()</code> method.
     *
     * @param name the internal name of the parameter
     * @param screenName the screen name associated with the parameter
     * @param toolTip the explanatory tool tip associated with the parameter
     */
    public EnumeratorParameter(String name, String screenName, String toolTip) {
        super(name, screenName, toolTip);

        texts = new ArrayList<String>();
        values = new ArrayList<Double>();
        initialize();
    }


    /**
     * This method adds a String-Integer association to this parameter. The string
     * will be displayed in the UI, the double is the corresponding value of the
     * parameter.
     *
     * @param value an double number to represent the value of the parameter
     * @param text a String to be displayed in the UI for this number
     */
    public void addChoice(double value, String text) {
        texts.add(text);
        values.add(new Double(value));

        comboBox.addItem(text);
    }


    public String getTextFor(double value) {
        int index = values.indexOf(new Double(value));
        return (String) texts.get(index);
    }


    private double getValueFor(String text) {
        int index = texts.indexOf(text);
        return ( (Double) values.get(index)).doubleValue();
    }


    /**
     * Returns the double number associated with the String currently selected.
     */
    public double getValue() {
        return getValueFor( (String) comboBox.getSelectedItem());
    }


    /**
     * Sets the String currently selected to be the one associated with the given
     * double number.
     */
    public void setValue(double value) {
        comboBox.setSelectedItem(getTextFor(value));
    }


    public void initialize() {
        super.initialize();

        comboBox = new JComboBox();
        comboBox.setToolTipText(toolTip);
        comboBox.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                if (e.getStateChange() == ItemEvent.SELECTED) {
                    fireParameterChange();
                }
            }
        });
    }


    public JComponent getValueComponent() {
        return comboBox;
    }


    public String valueAsString() {
        return Double.toString(getValue());
    }


    public Object valueAsObject() {
        return getValue();
    }

}
