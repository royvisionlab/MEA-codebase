package edu.ucsc.neurobiology.vision.parameters;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;


/**
 * This class implements the <code>Parameter</code> interface to represent a
 * set of colors in a ComboBox. The parameter willl contain by default all the
 * named colors in the <code>java.awt.Color</code> class.It also provides getter
 * and setter funtions for the current color stored in the parameter. Every time
 * the value changes a <code>PropetryChangeEvent</code> is dispatched to all listeners.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ColorParameter
    extends AbstractParameter {
    private JComboBox colorComboBox;
    private Color defaultColor;


    /**
     * Creates an instance of ColorParameter which will display as selected
     * the given color. The parameter willl contain by default all the named colors
     * in the <code>java.awt.Color</code> class.
     *
     * @param name the internal name of the parameter
     * @param screeName the screen name associated with the parameter
     * @param toolTip the explanatory tooltop associated with the parameter
     * @param color the default color selected in the parameter
     */
    public ColorParameter(String name, String screeName, String toolTip,
                          Color color) {
        super(name, screeName, toolTip);

        this.defaultColor = color;
        initialize();
    }


    public void initialize() {
        super.initialize();

        colorComboBox = new JComboBox(ColorComboBoxRenderer.colorMap.keySet().toArray());
        colorComboBox.setRenderer(new ColorComboBoxRenderer());
        colorComboBox.setSelectedItem(defaultColor);
        colorComboBox.setToolTipText(toolTip);

        colorComboBox.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                if (e.getStateChange() == ItemEvent.SELECTED) {
                    fireParameterChange();
                }
            }
        });
    }


    /**
     * Returns the current color selected in this parameter as a
     * <code>java.awt.Color</code> object.
     */
    public Color getColor() {
        return (Color) colorComboBox.getSelectedItem();
    }


    /**
     * Sets the current color selected in this parameter to be the given
     * <code>java.awt.Color</code> object.
     */
    public void setColor(Color color) {
        colorComboBox.setSelectedItem(color);
    }


    public JComponent getValueComponent() {
        return colorComboBox;
    }


    public String valueAsString() {
        Color c = getColor();
        return c.getRed() + ":" + c.getGreen() + ":" + c.getBlue();
    }


    public Object valueAsObject() {
        return getColor();
    }
}
