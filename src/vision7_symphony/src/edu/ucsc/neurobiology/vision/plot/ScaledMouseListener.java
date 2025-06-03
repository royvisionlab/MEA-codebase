package edu.ucsc.neurobiology.vision.plot;

import java.awt.*;
import java.awt.event.*;


/**
 * A callback interface to allow click notifications.
 *
 * @author etrusca Dumitru, University of California, Santa Cruz
 */
public interface ScaledMouseListener {

    /**
     * Called by the host (PlotPanel) when a click happened.
     *
     * @param source Component
     * @param event MouseEvent
     * @param x double the real space X coordinte
     * @param y double the real space Y coordinte
     */
    public void clickPerformed(Component source, MouseEvent event, double x, double y);
    
    public void pressPerformed(Component source, MouseEvent event, double x, double y);
    
    public void releasePerformed(Component source, MouseEvent event, double x, double y);

    public void enteredPerformed(Component source, MouseEvent event, double x, double y);
    
    public void exitedPerformed(Component source, MouseEvent event, double x, double y);
}
