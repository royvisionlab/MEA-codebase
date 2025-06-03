package edu.ucsc.neurobiology.vision.gui;

import java.util.*;

import java.awt.*;
import static java.awt.event.KeyEvent.*;
import java.awt.event.*;
import javax.swing.*;


/**
 *
 * Dumitru Petrusca, University of California, Santa Cruz
 */
public class VisionEventQueue
    extends EventQueue {

    private KeyListener listener;

    HashMap<KeyStroke, ActionListener> keyStrokeMap = new HashMap();

    public synchronized void setKeyListener(KeyListener listener) {
        this.listener = listener;
    }


    public void addKeyStroke(KeyStroke keyStroke, ActionListener listener) {
        keyStrokeMap.put(keyStroke, listener);
    }


    public void removeKeyStroke(KeyStroke keyStroke) {
        keyStrokeMap.remove(keyStroke);
    }


    protected void dispatchEvent(AWTEvent event) {
        if (event instanceof KeyEvent) {
            KeyEvent keyEvent = (KeyEvent) event;

            if (keyEvent.getID() == KEY_PRESSED) {
                if (listener != null) {
                    listener.keyPressed(keyEvent);
                }

                for (KeyStroke keyStroke : keyStrokeMap.keySet()) {
                    if (keyEvent.getKeyCode() == keyStroke.getKeyCode() &&

                        (keyEvent.getModifiers() & CTRL_MASK) ==
                        (keyStroke.getModifiers() & CTRL_MASK)
                        &&
                        (keyEvent.getModifiers() & SHIFT_MASK) ==
                        (keyStroke.getModifiers() & SHIFT_MASK)
                        &&
                        (keyEvent.getModifiers() & ALT_MASK) ==
                        (keyStroke.getModifiers() & ALT_MASK)

                        ) {

                        keyStrokeMap.get(keyStroke).actionPerformed(new ActionEvent(this,
                            0, ""));
                    }
                }
            }
        }

        super.dispatchEvent(event);
    }

}
