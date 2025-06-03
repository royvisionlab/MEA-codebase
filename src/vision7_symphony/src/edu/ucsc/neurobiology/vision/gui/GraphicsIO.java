package edu.ucsc.neurobiology.vision.gui;

import java.io.*;
import javax.imageio.*;

import java.awt.*;
import java.awt.image.*;
import javax.swing.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.parameters.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.plot.eps.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class GraphicsIO {
    static ImageSaveDialog fc = new ImageSaveDialog(null, Vision.getConfig());

    private GraphicsIO() {
    }


    public static void saveImage(JComponent c, Component parent, String title) {
        fc.setName(title);
        if (fc.showDialog(parent, title) == JFileChooser.APPROVE_OPTION) {
            File f = fc.getSelectedFile();
            SimpleFileFilter filter = (SimpleFileFilter) fc.getFileFilter();

            if (filter.getExtension().equals("eps")) {
                JScrollPane pane = (JScrollPane) fc.getAccessory();
                ParametersTable params = (ParametersTable) pane.getViewport().getView();
                try {
                    double width = params.getDoubleParameter("Width (inches)");
                    double height = params.getDoubleParameter("Height (inches)");
                    GraphicsIO.saveAsEPS(
                        c, f.getAbsolutePath(), width, height, false);
                } catch (IOException ex1) {
                    ex1.printStackTrace();
                }
            } else if (filter.getExtension().equals("png")) {
                try {
                    GraphicsIO.saveAsPNG(c, f, -1, -1);
                } catch (IOException ex1) {
                    ex1.printStackTrace();
                }
            }
        }
    }


    public static void saveAsPNG(Component c, File file, int w, int h) throws IOException {
        if (w == -1 && h == -1) {
            w = c.getWidth();
            h = c.getHeight();
        }

        Dimension d = c.getSize();
        c.setSize(w, h);

        BufferedImage bi = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);
        Graphics2D g2 = bi.createGraphics();
        c.paint(g2);
        ImageIO.write(bi, "PNG", file);

        c.setSize(d);
    }


    /*
     public void saveAsEPS(File file, double width, double height, boolean boxSize) throws
                IOException {

                Rectangle oldBounds = getBounds();
                this.replotAllData();

                // 72 is default dpi.
                if (boxSize) {
                    Insets i = this.getInsets();
                    setSize( (int) Math.round(width * 72.0 + i.left + i.right),
                            (int) Math.round(height * 72.0 + i.top + i.bottom));
                } else {
     setSize( (int) Math.round(width * 72.0), (int) Math.round(height * 72.0));
                }

                EpsGraphics2D g = new EpsGraphics2D(
                    "image", getXAxis().getLabelFont(), 0, 0, getWidth(), getHeight());
                g.setAccurateTextMode(true);
                g.setColorDepth(g.CYMK);

                paintBorder(g);
                spiPaintComponent(g, null);

                FileWriter writer = new FileWriter(file);
                writer.write(g.toString());
                writer.flush();
                writer.close();

                setBounds(oldBounds);
            }
     */

    public static void saveAsEPS(
        final JComponent c, final String name, final double width, final double height,
        boolean insetsSize) throws
        IOException {

        boolean isDoubleBuffered = false;
        ;
        if (c instanceof PlotPanel) {
            isDoubleBuffered = ( (PlotPanel) c).isDoubleBuffered();
            ( (PlotPanel) c).setDoubleBuffered(false);
        }

        Rectangle oldBounds = c.getBounds();

        if (insetsSize) {
            Insets i = c.getInsets();
            c.setSize( (int) Math.round(width * 72.0 + i.left + i.right),
                      (int) Math.round(height * 72.0 + i.top + i.bottom));
        } else {
            c.setSize( (int) Math.round(width * 72.0), (int) Math.round(height * 72.0));
        }
        c.validate();

        RepaintManager.currentManager(c).setDoubleBufferingEnabled(false);
        EpsGraphics2D g = new EpsGraphics2D(
            "", Font.decode("Arial ITALIC 8"), 0, 0, c.getWidth(), c.getHeight());
        g.setAccurateTextMode(true);
        c.paint(g);
        RepaintManager.currentManager(c).setDoubleBufferingEnabled(true);

        FileWriter writer = new FileWriter(
            StringUtil.removeExtension(name) + ".eps");
        writer.write(g.toString());
        writer.flush();
        writer.close();

        c.setBounds(oldBounds);
        c.validate();

        if (c instanceof PlotPanel) {
            ( (PlotPanel) c).setDoubleBuffered(isDoubleBuffered);
        }
    }


    public static void saveComponentToPNG(
        Component c, String fileName) throws IOException {
        saveComponentToFile(c, fileName, "PNG", -1, -1);
    }


    public static void saveComponentToFile(
        Component c, String fileName, String type, int w, int h) throws IOException {

        Dimension d = c.getSize();
        if (w != -1) {
            c.setSize(w, h);
        }

        BufferedImage bi = new BufferedImage(
            c.getWidth(), c.getHeight(), BufferedImage.TYPE_INT_RGB);
        Graphics2D g2 = bi.createGraphics();
        c.paint(g2);
        ImageIO.write(bi, type, new File(fileName));

        c.setSize(d);
    }


}
