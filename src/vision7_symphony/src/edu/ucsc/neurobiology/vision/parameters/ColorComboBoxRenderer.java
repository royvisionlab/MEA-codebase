package edu.ucsc.neurobiology.vision.parameters;

import java.util.*;

import java.awt.*;
import javax.swing.*;
import javax.swing.border.*;


/**
 * A swing renderer to draw the colors in the ColorParameters drop-down list.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
@SuppressWarnings("serial")
public class ColorComboBoxRenderer
    extends JLabel implements ListCellRenderer {

    public static Map<Color, String> colorMap;
    static {
        colorMap = new HashMap<Color, String>();
        colorMap.put(Color.black, "Black");
        colorMap.put(Color.blue, "Blue");
        colorMap.put(Color.cyan, "Cyan");
        colorMap.put(Color.darkGray, "DarkGray");
        colorMap.put(Color.gray, "Gray");
        colorMap.put(Color.green, "Green");
        colorMap.put(Color.lightGray, "LightGray");
        colorMap.put(Color.magenta, "Magenta");
        colorMap.put(Color.orange, "Orange");
        colorMap.put(Color.pink, "Pink");
        colorMap.put(Color.red, "Red");
        colorMap.put(Color.white, "White");
        colorMap.put(Color.yellow, "Yellow");
    }


    private ColorIcon icon = new ColorIcon();
    private Border selectedBorder;
    private Border deselectedBorder =
        BorderFactory.createEmptyBorder(2, 2, 2, 2);


    public Component getListCellRendererComponent(
        JList list, Object value, int index, boolean isSelected, boolean cellHasFocus) {

        Color color = (Color) value;
        icon.setColor(color);
        setIcon(icon);

        if (isSelected) {
            selectedBorder = BorderFactory.createLineBorder(
                list.getSelectionBackground(), 2);
            setBorder(selectedBorder);

            ToolTipManager.sharedInstance().setInitialDelay(1);
            ToolTipManager.sharedInstance().setReshowDelay(1);

            String message = (String) colorMap.get(color);
            if (message != null) {
                message += ": ";
            }
            list.setToolTipText(message +
                                color.getRed() + ", " + color.getGreen() + ", " +
                                color.getBlue());
        } else {
            setBorder(deselectedBorder);
        }

        return this;
    }

}


class ColorIcon
    implements Icon {

    private int w, h;
    private Color color;

    public ColorIcon() {
        this(Color.gray, 50, 12);
    }


    public ColorIcon(Dimension d) {
        this.w = d.width;
        this.h = d.height;
    }


    public ColorIcon(Color color, int w, int h) {
        this.color = color;
        this.w = w;
        this.h = h;
    }


    public void paintIcon(Component c, Graphics g, int x, int y) {
        g.setColor(Color.black);
        g.drawRect(x, y, w - 1, h - 1);
        g.setColor(color);
        g.fillRect(x + 1, y + 1, w - 2, h - 2);
    }


    public Color getColor() {
        return color;
    }


    public void setColor(Color color) {
        this.color = color;
    }


    public int getIconWidth() {
        return w;
    }


    public int getIconHeight() {
        return h;
    }

}
