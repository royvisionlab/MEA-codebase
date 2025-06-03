package edu.ucsc.neurobiology.vision.test;

import java.io.*;
import java.net.*;
import java.util.*;
import javax.help.*;

import javax.swing.*;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class SwingTest
    extends JPanel {


    public static void testUI() throws Exception {
        UIDefaults defaults = UIManager.getDefaults();
        Enumeration<Object> e = defaults.keys();
        while (e.hasMoreElements()) {
            Object o = e.nextElement();
            Object v = defaults.get(o);
            if (v != null) {
                System.err.println(o + " : " + v.getClass() + " : " + v);
            }
        }
    }


    private static URL findHelpSet(String helpset) {
        URL url = null;
        try {
            String dir = System.getProperty("user.dir");
            url = new URL("file:///" + dir + File.separator + helpset);
        } catch (Exception ex) {}

        if (url == null) {
            try {
                url = javax.help.HelpSet.findHelpSet(null, helpset);
            } catch (Exception ex) {}
        }
        return url;
    }


    public static void testHelp() {
        URL url = findHelpSet("UserHelp/userhelp.hs");
        try {
            javax.help.HelpSet helpSet = new javax.help.HelpSet(null, url);
            javax.help.JHelp helpViewer = new javax.help.JHelp(helpSet);

            JFrame f = new JFrame("Vision Help");
            ( (JSplitPane) helpViewer.getComponent(0)).setOrientation(
                JSplitPane.VERTICAL_SPLIT);

            f.add(helpViewer);

            f.setBounds(100, 100, 500, 500);
            f.setVisible(true);
        } catch (HelpSetException ex1) {
            ex1.printStackTrace();
        }
    }


    static void testFileDialog() throws Exception {
        String saveFolder = "c:\\";
    }


    public static void main(String[] args) throws Exception {
        try {
            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
        } catch (Exception e) {
            e.printStackTrace();
        }

        testFileDialog();
    }
}
