package edu.ucsc.neurobiology.vision.gui;

import java.util.*;

import java.awt.*;


/**
 *
 * Dumitru Petrusca, University of California, Santa Cruz
 */
public class PlaceLayout
    implements LayoutManager2 {

    private HashMap<Component, PlaceC> constr = new HashMap();


    public void addLayoutComponent(Component comp, Object constraints) {
        if (! (constraints instanceof PlaceC)) {
            throw new Error("wrong constraint");
        }

        constr.put(comp, (PlaceC) constraints);
    }


    public void layoutContainer(Container parent) {
        int w = parent.getWidth();
        int h = parent.getHeight();

        for (int i = 0; i < parent.getComponentCount(); i++) {
            Component comp = parent.getComponent(i);
            comp.setSize(comp.getPreferredSize());
            PlaceC c = constr.get(comp);
            if (c.wx != -1) {
                comp.setSize( (int) Math.round( (c.wx * w)), comp.getHeight());
            }
            if (c.wy != -1) {
                comp.setSize(comp.getWidth(), (int) Math.round( (c.wy * h)));
            }

            comp.setLocation(
                (int) Math.round(c.x * w - c.ax * comp.getWidth()),
                (int) Math.round(c.y * h - c.ay * comp.getHeight()));
        }
    }


    public Dimension minimumLayoutSize(Container parent) {
        return new Dimension(0, 0);
    }


    public Dimension preferredLayoutSize(Container parent) {
        return parent.getSize();
    }


    public float getLayoutAlignmentX(Container target) {
        return 0;
    }


    public float getLayoutAlignmentY(Container target) {
        return 0;
    }


    public void invalidateLayout(Container target) {
    }


    public Dimension maximumLayoutSize(Container target) {
        return new Dimension(0, 0);
    }


    public void addLayoutComponent(String name, Component comp) {
        throw new Error("not implemented");
    }


    public void removeLayoutComponent(Component comp) {
        constr.remove(comp);
    }
}
