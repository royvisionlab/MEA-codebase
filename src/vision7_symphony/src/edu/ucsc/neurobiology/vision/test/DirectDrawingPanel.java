package edu.ucsc.neurobiology.vision.test;

import java.util.*;

import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import javax.swing.*;
import javax.swing.Timer;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
@SuppressWarnings("serial")
public class DirectDrawingPanel
    extends JPanel {
    private BufferedImage bim; //this is the image
    private final int[] bank; //this is the buffer in which to draw (row-major)


    public DirectDrawingPanel() {
        //create an image with the specified size
        bim = new BufferedImage(400, 200, BufferedImage.TYPE_INT_RGB);
        Raster raster = bim.getWritableTile(0, 0);
        DataBufferInt dataBufferByte = (DataBufferInt) raster.getDataBuffer();
        bank = dataBufferByte.getData(0);
    }


    public void paint(Graphics g) {
        ( (Graphics2D) g).drawImage(bim, null, 0, 0);
    }


    public int[] getScreenBuffer() {
        return bank;
    }


    public static void main(String arg[]) {
        DirectDrawingPanel directPanel = new DirectDrawingPanel();

        JFrame f = new JFrame();
        f.add(directPanel);
        f.setBounds(100, 100, 450, 250);
        f.setVisible(true);

        // start animating the image
        Timer timer = new Timer(10, new Animator(directPanel));
        timer.start();
    }
}


class Animator
    implements ActionListener {
    private char color = 0;
    private char colorIncrement = +1;
    private DirectDrawingPanel directPanel;
    private int[] screenBuffer;

    public Animator(DirectDrawingPanel directPanel) {
        this.directPanel = directPanel;
        this.screenBuffer = directPanel.getScreenBuffer();
    }


    public void actionPerformed(ActionEvent e) {
        color += colorIncrement;
        if ( (color >= 255) || (color <= 0)) {
            colorIncrement *= -1;
        }
        Arrays.fill(screenBuffer, (byte) color);
        directPanel.repaint();
    }
}
