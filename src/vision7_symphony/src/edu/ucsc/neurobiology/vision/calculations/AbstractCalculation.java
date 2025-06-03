package edu.ucsc.neurobiology.vision.calculations;

import java.util.*;

import javax.swing.*;

import edu.ucsc.neurobiology.vision.*;


/**
 * The syperclass of all Vision Calculations.
 *
 * @author Charles A. Loomis, University of California, Santa Cruz
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public abstract class AbstractCalculation {

    /**
     * Called by Vision to set parameters and configure the calculation.
     * Called before startCalculation().
     *
     * @param parameters HashMap
     * @throws Exception
     */
    public abstract void setParameters(HashMap<String, String> parameters) throws Exception;


    /**
     * Called by Vision to start the calculation.
     * @throws Exception
     */
    public abstract void startCalculation() throws Exception;


    /**
     * The calculation can return a graphical component with diagnostic info.
     * @return JComponent
     */
    public JComponent getDiagnosticPanel() {
        return null;
    }


    /**
     * Must be called by the Calculation to let Vision know that it has ended.
     */
    protected void calculationDone() {
        Vision.getInstance().getCalculationManager().calculationDone();
    }
}
