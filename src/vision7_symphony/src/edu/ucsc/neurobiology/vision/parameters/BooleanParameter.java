package edu.ucsc.neurobiology.vision.parameters;

import java.awt.event.*;
import javax.swing.*;


/**
 * This class implements the <code>Parameter</code> interface to represent a
 * boolean value (<b>true</b> of <b>false</b>) as a CheckBox. It also provides
 * getter and setter funtions for the value of the parameter. Every time the value
 * changes a <code>PropetryChangeEvent</code> is dispatched to all listeners.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class BooleanParameter
    extends AbstractParameter {

    private JCheckBox status;
    private boolean defaultStatus;


    /**
     * Creates an instance of BooleanParameter which will display by default
     * the given status (selected or unselected).
     *
     * @param name the internal name of the parameter
     * @param screenName the screen name associated with the parameter
     * @param toolTip the explanatory tooltop associated with the parameter
     * @param defaultStatus the default status represented by the parameter
     */
    public BooleanParameter(String name, String screenName, String toolTip,
                            boolean defaultStatus) {

        super(name, screenName, toolTip);
        this.defaultStatus = defaultStatus;

        initialize();
    }


    /**
     * Returns <b>true</b> if the parameter is selected and <b>false</b> otherwise.
     */
    public boolean isSelected() {
        return status.isSelected();
    }


    /**
     * Sets the selection state of the parameter.
     *
     * @param selected boolean the selection state. If <b>true</b> the parameter
     * will be selected.
     */
    public void setSelected(boolean selected) {
        status.setSelected(selected);
    }


    public void initialize() {
        super.initialize();

        status = new JCheckBox("", defaultStatus);
        status.setOpaque(true);
        status.setToolTipText(toolTip);
        LookAndFeel.installColorsAndFont(status, "CheckBox.background",
                                         "CheckBox.foreground", "CheckBox.font");
        status.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                fireParameterChange();
            }
        });
    }


    public JComponent getValueComponent() {
        return status;
    }


    public String valueAsString() {
        return new Boolean(status.isSelected()).toString();
    }


    public boolean getValue() {
        return new Boolean(status.isSelected()).booleanValue();
    }


    public Object valueAsObject() {
        return status.isSelected();
    }
}
