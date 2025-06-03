package edu.ucsc.neurobiology.vision.gui;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;


/**
 *
 * Dumitru Petrusca, University of California, Santa Cruz
 */
public class MyGlassPane
    extends JComponent implements MouseListener, MouseMotionListener {

    Toolkit toolkit;
//    JMenuBar menuBar;
    Container contentPane;
    Font font = new Font("sans", Font.ITALIC, 24);


    public MyGlassPane(Container contentPane) {
//        this.menuBar = frame.getJMenuBar();
        this.contentPane = contentPane;

        toolkit = Toolkit.getDefaultToolkit();
        addMouseListener(this);
        addMouseMotionListener(this);
    }


    protected void paintComponent(Graphics g) {
        if (this.isVisible()) {
            g.setFont(font);
            FontMetrics fm = g.getFontMetrics(font);
            String s = "The user interface is temporarily disabled.";

            int w = fm.stringWidth(s);
            int h = fm.getHeight();
            int x = getWidth() / 2 - w / 2;
            int y = getHeight() / 2 - h / 2;

            float c = 0.8f;
            g.setColor(new Color(c, c, c, 1f));
            g.fillRect(x - 10, y, w + 20, h);
//            g.fillRect(0, 0, getWidth(), getHeight());
            g.setColor(Color.black);
            g.drawString("The user interface is temporarily disabled.", x, y + h - 7);
        }
    }


    public void mouseMoved(MouseEvent e) {
        redispatchMouseEvent(e, false);
    }


    public void mouseDragged(MouseEvent e) {
        redispatchMouseEvent(e, false);
    }


    public void mouseClicked(MouseEvent e) {
        redispatchMouseEvent(e, false);
    }


    public void mouseEntered(MouseEvent e) {
        redispatchMouseEvent(e, false);
    }


    public void mouseExited(MouseEvent e) {
        redispatchMouseEvent(e, false);
    }


    public void mousePressed(MouseEvent e) {
        redispatchMouseEvent(e, false);
    }


    public void mouseReleased(MouseEvent e) {
        redispatchMouseEvent(e, true);
    }


    //A more finished version of this method would
    //handle mouse-dragged events specially.
    private void redispatchMouseEvent(MouseEvent e, boolean repaint) {
        Point glassPanePoint = e.getPoint();
        Container container = contentPane;
        Point containerPoint = SwingUtilities.convertPoint(
            this, glassPanePoint, contentPane);

        /*
                if (containerPoint.y < 0) { //we're not in the content pane
                    if (containerPoint.y + menuBar.getHeight() >= 0) {
                        //The mouse event is over the menu bar.
                        //Could handle specially.
                    } else {
                        //The mouse event is over non-system window
                        //decorations, such as the ones provided by
                        //the Java look and feel.
                        //Could handle specially.
                    }

                } else {
                    //The mouse event is probably over the content pane.
                    //Find out exactly which component it's over.
//            Component component =
//                SwingUtilities.getDeepestComponentAt(
//                container,
//                containerPoint.x,
//                containerPoint.y);

                    if (component != null && component.equals(liveButton)) {
                        //Forward events over the check box.
                        Point componentPoint = SwingUtilities.convertPoint(
                            glassPane,
                            glassPanePoint,
                            component);
                        component.dispatchEvent(new MouseEvent(component,
                            e.getID(),
                            e.getWhen(),
                            e.getModifiers(),
                            componentPoint.x,
                            componentPoint.y,
                            e.getClickCount(),
                            e.isPopupTrigger()));
                    }
                }
         */

        //Update the glass pane if requested.
        if (repaint) {
            repaint();
        }
    }
}
