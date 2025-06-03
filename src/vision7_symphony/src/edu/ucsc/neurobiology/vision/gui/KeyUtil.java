package edu.ucsc.neurobiology.vision.gui;

import java.awt.event.*;


/**
 *
 * Dumitru Petrusca, University of California, Santa Cruz
 */
public class KeyUtil {

    private static final int[] keyCodes = {
                                          KeyEvent.VK_1,
                                          KeyEvent.VK_2,
                                          KeyEvent.VK_3,
                                          KeyEvent.VK_4,
                                          KeyEvent.VK_5,
                                          KeyEvent.VK_6,
                                          KeyEvent.VK_7,
                                          KeyEvent.VK_8,
                                          KeyEvent.VK_9,
                                          KeyEvent.VK_Q,
                                          KeyEvent.VK_W,
                                          KeyEvent.VK_E,
                                          KeyEvent.VK_R,
                                          KeyEvent.VK_T,
                                          KeyEvent.VK_Y,
                                          KeyEvent.VK_U,
                                          KeyEvent.VK_I,
                                          KeyEvent.VK_O,
                                          KeyEvent.VK_P,
    };

    private static final String[] keyTexts = new String[keyCodes.length];


    static {
        for (int i = 0; i < keyTexts.length; i++) {
            keyTexts[i] = KeyEvent.getKeyText(keyCodes[i]);
        }
    }


    public static int getKeyCode(int i) {
        return keyCodes[i];
    }


    public static String getKeyText(int i) {
        return keyTexts[i];
    }
}
