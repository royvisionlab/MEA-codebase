package edu.ucsc.neurobiology.vision.parameters;

import javax.swing.*;
import javax.swing.event.*;


/**
 * This class implements the <code>Parameter</code> interface to represent a
 * Java integer in a JSpinner. It also provides getter and setter funtions for the
 * value of the parameter. Every time the value changes a
 * <code>PropetryChangeEvent</code> is dispatched to all listeners.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class IntegerParameter
    extends AbstractParameter {
    private JSpinner spin;


    /**
     * Creates an instance of IntegerParameter which will display by default
     * the given value.
     *
     * @param name the internal name of the parameter
     * @param screenName the screen name associated with the parameter
     * @param toolTip the explanatory tooltop associated with the parameter
     * @param value the default value of the parameter
     * @param min int the minimum allowed value
     * @param max int the maximum allowed value
     */
    public IntegerParameter(String name, String screenName, String toolTip, int value,
                            int min, int max) {
        super(name, screenName, toolTip);
        initialize(min, max, value);
    }


    /**
     * Returns the value currently stored by this parameter.
     */
    public int getValue() {
        return (Integer) spin.getValue();
    }


    public boolean isValid() {
        return true;
    }


    public void initialize(int min, int max, int value) {
        super.initialize();

        spin = new JSpinner(new SpinnerNumberModel(value, min, max, 1));
        ( (javax.swing.JSpinner.DefaultEditor) spin.getEditor()).getTextField().
            setHorizontalAlignment(JFormattedTextField.LEFT);

        spin.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                fireParameterChange();
            }
        });
    }


    public JComponent getValueComponent() {
        return spin;
    }


    public String valueAsString() {
        return spin.getValue().toString();
    }


    public Object valueAsObject() {
        return spin.getValue();
    }

}
