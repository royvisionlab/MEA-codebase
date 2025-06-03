package edu.ucsc.neurobiology.vision.parameters;

import java.beans.*;

import javax.swing.*;


/**
 * This interface defines the minimum set of functionality required for any parameter.
 * This include the getter methods for various names of the parameter and for the two
 * graphical components used to display the Name and the Value of the parameter. Any
 * should accept PropertyChangeListeners to be added to it and be notified
 * of changes of the parameter. However the getter functions for the value of the
 * parameter are left to the concrete implementations.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public interface Parameter {

    /**
     * Returns the screen name of the parameter as a String. This name is displayed
     * in the UI (can contain Unicode)
     */
    public String getName();


    public String getScreenName();


    /**
     * Returns the explanatory tooltip text associated with the parameter.
     */
    public String getToolTip();


    /**
     * Returns the graphical component used to display the Name of this parameter
     * in the UI. This method is for internal use and should usually not be called
     * by the user.
     */
    public JComponent getNameComponent();


    /**
     * Returns the graphical component used to display the Value of this parameter
     * in the UI. This method is for internal use and should usually not be called
     * by the user.
     */
    public JComponent getValueComponent();


    /**
     * Adds a specified PropertyChangeListener to this parameter in order to be
     * notified when any change of the parameter happens.
     */
    public void addPropertyChangeListener(PropertyChangeListener listener);


    /**
     * Sets whether the parameter is enabled (true) or disabled (false).
     */
    public void setEnabled(boolean enabled);


    /**
     * Sets whether the parameter is enabled (true) or disabled (false).
     */
    public boolean isValid();


    /**
     * Returns a String representation of the value of the parameter.
     */
    public String valueAsString();


    public Object valueAsObject();
}
