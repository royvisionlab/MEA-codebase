package edu.ucsc.neurobiology.vision.gui;

import java.beans.*;
import java.util.*;
import java.util.List;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;


/**
 *
 * @author Charles A. Loomis, University of California, Santa Cruz
 * Dumitru Petrusca, University of California, Santa Cruz
 */
public class GridDesktopManager
    extends DefaultDesktopManager implements ContainerListener, ActionListener,
    ComponentListener {

    enum LayoutPolicy {
        HORIZONTAL_TILING, VERTICAL_TILING, GRID_TILING, STACK, FREE
    }


    private Desktop parent;

    // This keeps a mapping of an internal frame to a menu item and vice
    // versa.  Consequently, two entries exist for each mapping.
    private Map<Component,Component> viewMenuHash;

    // This is the menu which keeps track of all of the current frames.
    private JMenu viewMenu;

    // This is the button group which keeps track of the selected frame.
    private ButtonGroup buttonGroup;

    private ActionListener menuListener;

    private LayoutPolicy layoutPolicy = LayoutPolicy.GRID_TILING;


    public static String formatString(String s) {
        StringBuffer sb = new StringBuffer(s.trim());

        for (int i = 0; i < sb.length(); i++) {
            if (i == 0) {
                sb.setCharAt(i, Character.toUpperCase(sb.charAt(i)));
            } else if (sb.charAt(i) == '_') {
                sb.setCharAt(i, ' ');

                sb.setCharAt(i + 1, Character.toUpperCase(sb.charAt(i + 1)));
                i++;
            } else {
                sb.setCharAt(i, Character.toLowerCase(sb.charAt(i)));
            }
        }

        return sb.toString();
    }


    public GridDesktopManager(Desktop parent) {
        this.parent = parent;
        menuListener = new MenuListener();

        // Make those objects necessary to keep track of the current views.
        viewMenu = new JMenu("Window");
        viewMenu.setMnemonic(KeyEvent.VK_W);
//        viewMenu.setEnabled(false);

        ButtonGroup tileGroup = new ButtonGroup();
        for (LayoutPolicy v : LayoutPolicy.values()) {
            String name = formatString(v.toString());
            JMenuItem item = new JRadioButtonMenuItem(name);
            item.setMnemonic(name.charAt(0));
            item.setAccelerator(KeyStroke.getKeyStroke(
                name.charAt(0), KeyEvent.ALT_MASK + KeyEvent.CTRL_MASK));
            item.addActionListener(this);
            tileGroup.add(item);
            viewMenu.add(item);

            if (v == layoutPolicy) {
                item.setSelected(true);
            }
        }

        viewMenu.add(new JSeparator());
        buttonGroup = new ButtonGroup();
        viewMenuHash = new HashMap();
    }


    /**
     * Called when the layout policy changes
     */
    public void actionPerformed(ActionEvent e) {
        String command = e.getActionCommand();

        for (LayoutPolicy v : LayoutPolicy.values()) {
            if (command.equals(formatString(v.toString()))) {
                layoutPolicy = v;
            }
        }

        updateLayout();
    }


    private static JInternalFrame getClosestTo(List<JInternalFrame> frames, int x, int y) {
        if (frames == null || frames.size() == 0) {
            throw new IllegalArgumentException("The list of frames is null or empty");
        }

        double minDist2 = Double.POSITIVE_INFINITY;
        int minDistIndex = -1;

        for (int i = 0; i < frames.size(); i++) {
            Point p = frames.get(i).getLocation();
            double dist2 = (p.x - x) * (p.x - x) + (p.y - y) * (p.y - y);
            if (dist2 < minDist2) {
                minDist2 = dist2;
                minDistIndex = i;
            }
        }

        return frames.get(minDistIndex);
    }


    public void updateLayout() {
        JInternalFrame[] a = parent.mainDesktop.getAllFrames();
        LinkedList<JInternalFrame> frames = new LinkedList<JInternalFrame>();
        for (JInternalFrame f : a) {
            if (!f.isIcon() && f.isVisible()) {
                frames.add(f);
            }
        }
        int n = frames.size();
        if (n == 0) {
            return;
        }

        int w = parent.mainDesktop.getWidth();
        int h = parent.mainDesktop.getHeight();

        switch (layoutPolicy) {
            case HORIZONTAL_TILING:
                int dw = w / n;
                for (int i = 0; frames.size() != 0; i++) {
                    JInternalFrame f = getClosestTo(frames, i * dw, 0);
                    frames.remove(f);
                    f.reshape(i * dw, 0, dw, h);
                }
                break;

            case VERTICAL_TILING:
                int dh = h / n;
                for (int i = 0; frames.size() != 0; i++) {
                    JInternalFrame f = getClosestTo(frames, 0, i * dh);
                    frames.remove(f);
                    f.reshape(0, i * dh, w, dh);
                }
                break;

            case GRID_TILING:
                int nx = (int) Math.ceil(Math.sqrt(n));
                int dx = w / nx;
                int ny = (int) Math.ceil( ( (double) n) / nx);
                int dy = h / ny;

                for (int i = 0; frames.size() != 0; i++) {
                    int y = i / nx;
                    int x = i % nx;
                    JInternalFrame f = getClosestTo(frames, x * dx, y * dy);
                    frames.remove(f);
                    f.reshape(x * dx, y * dy, dx, dy);
                }
                break;

            case STACK:
                for (int i = 0; i < frames.size(); i++) {
                    JInternalFrame f = frames.get(i);
                    f.reshape(0, 0, w, h);
                }
                break;
        }

    }


    public void activateFrame(JInternalFrame frame) {
        parent.setCurrentFrame(frame);

        // select the corresponding menu item
        if (viewMenuHash.containsKey(frame)) {
            JRadioButtonMenuItem item = (JRadioButtonMenuItem) viewMenuHash.get(frame);
            item.setSelected(true);
        }
        if (frame != null) {
            super.activateFrame(frame);
        }
    }


    public void closeFrame(JInternalFrame frame) {
        parent.setCurrentFrame(null);

        // See if this item is in frame menu.
        if (viewMenuHash.containsKey(frame)) {
            JRadioButtonMenuItem item = (JRadioButtonMenuItem) viewMenuHash.get(frame);
            viewMenu.remove(item);
            buttonGroup.remove(item);

            viewMenuHash.remove(frame);
            viewMenuHash.remove(item);

            // Disable the menu if nothing is left.
//            if (viewMenuHash.size() == 0) {
//                viewMenu.setEnabled(false);
//            }
        }

        super.closeFrame(frame);
        updateLayout();
    }


    /**
     * This desktop manager is also a container listener so that it can detect
     * when new internal frames are added and update the current view menu.
     * Note that the JLayeredPane moves a particular container to the front by
     * first removing it and then adding it again.  As a result, the container
     * additions must be handled here, and the removals should be handled in
     * the closeFrame() method. */
    public void componentAdded(ContainerEvent event) {
        // Get the component which has been added to the desktop.
        Component c = event.getChild();

        // Add this to the frame menu.
        if (c instanceof JInternalFrame) {
            JInternalFrame jif = (JInternalFrame) c;
            
            // Only do something if this component doesn't
            // already exist in the hash table.
            if (!viewMenuHash.containsKey(jif)) {
                String description = jif.getTitle();

                JRadioButtonMenuItem item = new JRadioButtonMenuItem();
                item.setText(description);
                item.addActionListener(menuListener);

                viewMenu.add(item);
                buttonGroup.add(item);
                viewMenuHash.put(jif, item);
                viewMenuHash.put(item, c);
//                viewMenu.setEnabled(true);

                if (layoutPolicy == LayoutPolicy.FREE) {
                    ( (JInternalFrame) c).reshape(0, 0, parent.getWidth(),
                                                  parent.getHeight());
                }
            }
        }
    }


    /**
     * This method does nothing.  When internal frames are removed from the
     * desktop, the menu is updated in the closeFrame() method. */
    public void componentRemoved(ContainerEvent event) {
        if (event.getChild() instanceof JInternalFrame) {
            closeFrame( (JInternalFrame) event.getChild());
        }
    }


    /**
     * Return a reference to the menu which keeps a list of the current
     * views. */
    public JMenu getCurrentViewsMenu() {
        return viewMenu;
    }


    /**
     * This inner class handles the actionPerformed() method calls from the menu
     * items associated with each internal frame.
     */
    private class MenuListener
        implements ActionListener {
        public void actionPerformed(ActionEvent event) {
            JRadioButtonMenuItem source = (JRadioButtonMenuItem) event.getSource();
            JInternalFrame frame = (JInternalFrame) viewMenuHash.get(source);
            try {
                if (frame != null) {
                    frame.setSelected(true);
                }
            } catch (PropertyVetoException pve) {}
        }
    }


    public void componentHidden(ComponentEvent e) {
        updateLayout();
    }


    public void componentMoved(ComponentEvent e) {}


    public void componentShown(ComponentEvent e) {
        updateLayout();
    }


    public void componentResized(ComponentEvent e) {
        updateLayout();
    }

}
