package edu.ucsc.neurobiology.vision.gui;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import javax.swing.*;

import edu.ucsc.neurobiology.vision.Vision;
import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * Dumitru Petrusca, University of California, Santa Cruz
 * owner: Matthew Grivich, The Salk Institute
 */
public class StartupWindow {
    JWindow w;
    JLabel progress;


    public StartupWindow() {
        GraphicsEnvironment ge =
            GraphicsEnvironment.getLocalGraphicsEnvironment();
        GraphicsDevice gd = ge.getScreenDevices()[0];
        DisplayMode mode = gd.getDisplayMode();

        int screenWidth = mode.getWidth();
        int screenHeight = mode.getHeight();
        int width = 450;
        int height = 300;

        w = new JWindow();
        w.setAlwaysOnTop(true);
        w.setBounds(screenWidth / 2 - width / 2, screenHeight / 2 - height / 2, width,
                    height);
        w.getContentPane().setLayout(new BorderLayout());

        JPanel pCenter = new JPanel(new GridLayout(0, 1));
        pCenter.setBackground(Color.WHITE);
        ImageIcon splashImage = VisionImage.getIcon("splash.gif");
        
        String buildDateString = Vision.buildDateString();
        if (buildDateString == null) {
            buildDateString = "";
        } else {
            buildDateString = "Build: " + buildDateString + "<br>";
        }
        
        JLabel title = new JLabel(
                "<html><br><br><br><br><br><br><br><br><br><br><br><br>" + 
                "Version: " + VisionParams.versionString() + "<br>"
                + buildDateString
                + "Created by Dumitru Petrusca, Matthew Grivich, and Tim Machado          " + "<br>"
                + "Copyright Regents of the University of California, Santa Cruz</html>");
        title.setIcon(splashImage);
        pCenter.add(title);
        title.setHorizontalTextPosition(JLabel.CENTER);
        title.setVerticalTextPosition(JLabel.CENTER);

        w.add(pCenter, BorderLayout.CENTER);

        progress = new JLabel("Starting up...", JLabel.CENTER);
        progress.setForeground(Color.WHITE);
        JPanel pSouth = new JPanel(new GridLayout(0, 1));
        
        pSouth.setBackground(Color.BLACK);
        pSouth.add(progress);
        w.add(pSouth, BorderLayout.SOUTH);
        
        w.addMouseListener(new MouseListener() {
            public void mouseClicked(MouseEvent e) {
                hide();
            }
            public void mouseEntered(MouseEvent e)  {};
            public void mouseExited(MouseEvent e)   {};
            public void mouseReleased(MouseEvent e) {};
            public void mousePressed(MouseEvent e)  {};
        });
    }



    public void show() {
        w.setVisible(true);
    }


    public void setProgress(String text) {
        progress.setText(text);
    }


    public void hide() {
        w.setVisible(false);
    }
}
