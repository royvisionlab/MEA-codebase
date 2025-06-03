package edu.ucsc.neurobiology.vision.gui;

import java.io.*;

import java.awt.*;
import java.awt.datatransfer.*;
import java.awt.image.*;


/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2002</p>
 *
 * <p>Company: </p>
 *
 * @author not attributable
 * @version 1.0
 */
public class Clipboard {

    private Clipboard() {
    }


    public static void saveComponentToClipboard(Component c) {
        BufferedImage bi = new BufferedImage(
            c.getWidth(), c.getHeight(), BufferedImage.TYPE_INT_RGB);
        Graphics2D g2 = bi.createGraphics();
        c.paint(g2);
        setClipboard(bi);
    }


    /**
     * This method writes a image to the system clipboard. Otherwise it returns null.
     */
    public static void setClipboard(Image image) {
        ImageSelection imgSel = new ImageSelection(image);
        Toolkit.getDefaultToolkit().getSystemClipboard().setContents(imgSel, null);
    }


    // This class is used to hold an image while on the clipboard.
    public static class ImageSelection
        implements Transferable {
        private Image image;

        public ImageSelection(Image image) {
            this.image = image;
        }


        // Returns supported flavors
        public DataFlavor[] getTransferDataFlavors() {
            return new DataFlavor[] {DataFlavor.imageFlavor};
        }


        // Returns true if flavor is supported
        public boolean isDataFlavorSupported(DataFlavor flavor) {
            return DataFlavor.imageFlavor.equals(flavor);
        }


        // Returns image
        public Object getTransferData(DataFlavor flavor) throws
            UnsupportedFlavorException, IOException {
            if (!DataFlavor.imageFlavor.equals(flavor)) {
                throw new UnsupportedFlavorException(flavor);
            }
            return image;
        }
    }
}
