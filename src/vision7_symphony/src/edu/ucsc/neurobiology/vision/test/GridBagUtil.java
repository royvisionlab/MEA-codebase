package edu.ucsc.neurobiology.vision.test;

import java.awt.*;
import javax.swing.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */

public class GridBagUtil {

    public static void addVSpaceEater(Container cont, int gridx, int gridy) {
        JLabel l = new JLabel();
        GridBagConstraints c = new GridBagConstraints();
        c.gridx = gridx;
        c.gridy = gridy;
        c.weighty = 1;
        ( (GridBagLayout) cont.getLayout()).setConstraints(l, c);
        cont.add(l);
    }


    public static void add(Container cont, Component comp, int gridx, int gridy) {
        GridBagConstraints c = new GridBagConstraints();
        c.gridx = gridx;
        c.gridy = gridy;
        ( (GridBagLayout) cont.getLayout()).setConstraints(comp, c);
        cont.add(comp);
    }


    public static void add(Container cont, Component comp,
                           int gridx, int gridy, int gridwidth, int gridheight) {

        GridBagConstraints c = new GridBagConstraints();
        c.gridx = gridx;
        c.gridy = gridy;
        c.gridwidth = gridwidth;
        c.gridheight = gridheight;
        ( (GridBagLayout) cont.getLayout()).setConstraints(comp, c);
        cont.add(comp);
    }


    public static void add(Container cont, Component comp,
                           int gridx, int gridy, int gridwidth, int gridheight,
                           double weightx, double weighty) {

        GridBagConstraints c = new GridBagConstraints();
        c.gridx = gridx;
        c.gridy = gridy;
        c.gridwidth = gridwidth;
        c.gridheight = gridheight;
        c.weightx = weightx;
        c.weighty = weighty;
        ( (GridBagLayout) cont.getLayout()).setConstraints(comp, c);
        cont.add(comp);
    }


    public static void add(Container cont, Component comp,
                           int gridx, int gridy, int gridwidth, int gridheight,
                           double weightx, double weighty, int anchor, int fill) {

        GridBagConstraints c = new GridBagConstraints();
        c.gridx = gridx;
        c.gridy = gridy;
        c.gridwidth = gridwidth;
        c.gridheight = gridheight;
        c.weightx = weightx;
        c.weighty = weighty;
        c.anchor = anchor;
        c.fill = fill;
        ( (GridBagLayout) cont.getLayout()).setConstraints(comp, c);
        cont.add(comp);
    }


    public static void addLine(Container cont, Component comp2, int gridy) {
        GridBagLayout layout = (GridBagLayout) cont.getLayout();
        GridBagConstraints c = new GridBagConstraints();

        c.gridy = gridy;
        c.gridx = 0;
        c.anchor = GridBagConstraints.WEST;
        c.weightx = 0;
        c.gridwidth = 2;
        layout.setConstraints(comp2, c);
        cont.add(comp2);
    }


    public static void addLine(Container cont, String text, Component comp2, int gridy) {
        GridBagLayout layout = (GridBagLayout) cont.getLayout();
        GridBagConstraints c = new GridBagConstraints();
        c.gridy = gridy;

        c.gridx = 0;
        c.anchor = GridBagConstraints.WEST;
        c.weightx = 1;
        JLabel label = new JLabel(text);
        layout.setConstraints(label, c);
        cont.add(label);

        c.gridx = 1;
        c.anchor = GridBagConstraints.EAST;
        c.weightx = 0;
        layout.setConstraints(comp2, c);
        cont.add(comp2);
    }

}
