package edu.ucsc.neurobiology.vision.gui;

import java.awt.*;
import javax.swing.*;

import org.jdesktop.swingx.*;


/**
 *
 * @author Charles A. Loomis, University of California, Santa Cruz
 * Dumitru Petrusca, University of California, Santa Cruz
 */
public class ViewManager
    extends JXTaskPane {

    private JComponent activeController;


    public ViewManager() {
        this.setTitle("Controller");
        this.setScrollOnExpand(true);
    }


    public void setCurrentController(Object controller, String name) {
        if (activeController != null) {
            this.remove(activeController);
            this.setTitle("Controller");
        }

        if (controller == null) {
            activeController = null;
            this.setTitle("Controller");
            repaint();
        } else if (controller instanceof JComponent) {
            activeController = (JComponent) controller;
            this.add(activeController, BorderLayout.CENTER);
            this.setTitle("Controller: " + name);
            repaint();
        } else {
            throw new IllegalArgumentException(
                "Wrong controller: " + controller.getClass());
        }
    }

}
