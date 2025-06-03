package edu.ucsc.neurobiology.vision.plot;

import java.io.*;

import javax.swing.*;


/**
 * Defines a possible action that can be takes on a selection in PlotPanel.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public abstract class SelectionAction {
    private String description;

    public SelectionAction(String description) {
        this.description = description;
    }

    public String getDescription() {
        return description;
    }

    /**
     * Must do the real action.
     *
     * @param source JComponent
     * @param selection Selection
     */
    public abstract void selectionPerformed(JComponent source, Selection selection);
}
