package edu.ucsc.neurobiology.vision.gui;

import javax.swing.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class VisionImage {

    // only static methods
    private VisionImage() {
    }


    public static ImageIcon getIcon(String name) {
        return new ImageIcon(VisionImage.class.getResource(name));

    }

}
