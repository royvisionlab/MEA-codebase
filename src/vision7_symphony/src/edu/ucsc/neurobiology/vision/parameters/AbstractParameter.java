package edu.ucsc.neurobiology.vision.parameters;

import java.beans.*;
import java.util.*;

import java.awt.*;
import javax.swing.*;


/**
 * The class is used as a default implementation of the Parameter interface making
 * it easier to create a concrete parameter. All usual parameters may extend this
 * class except fot the unusul ones which may implement Parameter directly.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public abstract class AbstractParameter
    implements Parameter {

    public static final Color normalForeground = Color.black;
    public static final Color normalBackground = Color.white;
    public static final Color errorForeground = Color.white;
    public static final Color errorBackground = Color.red;


    /**
     * This name is displayed in the user interface.
     */
    protected String name;
    protected String screenName;

    /**
     * The explanatory tooltip text associated with the parameter.
     */
    protected String toolTip;


    private JLabel label;
    private boolean isInitialized;
    private ArrayList<PropertyChangeListener> listeners = new ArrayList<PropertyChangeListener>();


    /**
     * Creates an instance of AbstractParameter with a given internalName,
     * screenName and toolTip.
     */
    protected AbstractParameter(String name, String screenName, String toolTip) {
        if ( (name == null) || (name.trim().length() == 0)) {
            throw new IllegalArgumentException("wrong name");
        } else {
            this.name = name.trim();
        }

        if ( (screenName == null) || (screenName.trim().length() == 0)) {
            this.screenName = new String(this.name);
        } else {
            this.screenName = screenName;
        }

        if ( (toolTip == null) || (toolTip.trim().length() == 0)) {
            this.toolTip = new String(this.screenName);
        } else {
            this.toolTip = toolTip;
        }
    }


    public void setEnabled(boolean enabled) {
        this.getNameComponent().setEnabled(enabled);
        this.getValueComponent().setEnabled(enabled);
    }


    public void addPropertyChangeListener(PropertyChangeListener listener) {
        listeners.add(listener);
    }


    /**
     * This method is called to notify all the listeners about the change in
     * the value of the parameter
     */
    protected void fireParameterChange() {
        PropertyChangeEvent event = new PropertyChangeEvent(this, name, this, this);
        for (int i = 0; i < listeners.size(); i++) {
            PropertyChangeListener listener = (PropertyChangeListener) listeners.get(i);
            listener.propertyChange(event);
        }
    }


    /**
     * Called by the User Interface in order to initialize the graphics
     * of the parameter (if it will be displayed in the UI).
     */
    protected void initialize() {
        label = new JLabel(screenName);
        label.setOpaque(true);
        label.setToolTipText(toolTip);

//        LookAndFeel.installColorsAndFont(label,
//            "CheckBox.background","CheckBox.foreground","CheckBox.font");

        isInitialized = true;
    }


    /**
     * Returns <b>true</b> if the UI part of this parameter is created,
     * <b>false</b> otherwise.
     */
    public boolean isInitialized() {
        return this.isInitialized;
    }


    public String getName() {
        return name;
    }


    public String getScreenName() {
        return screenName;
    }


    public String getToolTip() {
        return toolTip;
    }


    public JComponent getNameComponent() {
        return label;
    }


    public boolean isValid() {
        return true;
    }


    public String toString() {
        return name;
    }

}
