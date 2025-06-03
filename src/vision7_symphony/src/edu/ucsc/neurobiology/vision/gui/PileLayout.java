package edu.ucsc.neurobiology.vision.gui;

import java.awt.*;


/**
 *
 * Dumitru Petrusca, University of California, Santa Cruz
 */
public class PileLayout
    implements LayoutManager {


    public PileLayout() {
    }


    public void addLayoutComponent(String parm1, Component parm2) {
    }


    public void removeLayoutComponent(Component parm1) {
    }


    public Dimension minimumLayoutSize(Container parm1) {
        return new Dimension(0, 0);
    }


    public Dimension preferredLayoutSize(Container parm1) {
        return new Dimension(0, 0);
    }


    public void layoutContainer(Container container) {
        int y = 0;
        for (int i = 0; i < container.getComponentCount() - 1; i++) {
            Component c = container.getComponent(i);
            Dimension d = c.getPreferredSize();
            c.setBounds(0, y, container.getWidth(), d.height);
            y += d.height;
        }

        Component c = container.getComponent(container.getComponentCount() - 1);
        c.setBounds(0, y, container.getWidth(), container.getHeight() - y);
    }

}
