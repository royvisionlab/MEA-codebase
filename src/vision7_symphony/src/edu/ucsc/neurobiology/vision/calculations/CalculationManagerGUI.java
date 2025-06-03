package edu.ucsc.neurobiology.vision.calculations;

import java.util.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

import org.jdesktop.swingx.*;
import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.parameters.*;


/**
 * This manages a set of AbstractCalculations.  It allows the user to change
 * the calculation parameters and to run calculations interactively.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
@SuppressWarnings("serial")
public class CalculationManagerGUI
    extends JXTaskPane implements ActionListener {

    private JMenu menu;
    CalculationManager calculationManager;


    public CalculationManagerGUI(CalculationManager calculationManager) {
        //              this.setTitle("Data Manager");
        this.setTitle("Calculation");
        this.setScrollOnExpand(true);

        this.calculationManager = calculationManager;

        menu = new JMenu("Calculations");
        menu.setMnemonic(KeyEvent.VK_C);
        for (String cName : calculationManager.calculations.keySet()) {
            if (calculationManager.calculations.get(cName) == null) {
                menu.add(new JSeparator());
            } else {
                JMenuItem item = new JMenuItem(cName);
                item.addActionListener(this);
                menu.add(item);
            }
        }
    }


    public JMenu getMenu() {
        return menu;
    }


    public void setDiagnosticsPanel(AbstractCalculation calculation, String name) {
        this.setTitle("Calculation: " + name);

        JComponent dp = calculation.getDiagnosticPanel();
        removeAll();
        if (dp != null) {
            dp.setBorder(BorderFactory.createEmptyBorder(5, 10, 5, 10));
            add(dp, BorderLayout.CENTER);
        }
    }


    public void removeDiagnosticsPanels() {
        if (SwingUtilities.isEventDispatchThread()) {
            setTitle("Calculation");
            removeAll();
        } else {
            try {
                SwingUtilities.invokeAndWait(new Runnable() {
                    public void run() {
                        setTitle("Calculation");
                        removeAll();
                    }
                });
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }
    }


    public void actionPerformed(ActionEvent event) {
        if (!event.getActionCommand().equals("-")) {
            // do not start the calculation if there is one running
            if (calculationManager.isCurrentlyCalculating()) {
                Vision.getInstance().sendMessage(
                    "Another calculation is already being executed.");
                return;
            }

            String cName = event.getActionCommand();
            ParametersTable table = Vision.getConfig().showDialog(cName, cName, Vision.getInstance().getMainFrame());
            if (table == null) return;

            LinkedHashMap<String, String> params = new LinkedHashMap<String, String>();
            for (int i = 0; i < table.getParametersCount(); i++) {
                Parameter p = (Parameter) table.getParameter(i);
                if (p != null) {
                    params.put(p.getName(), p.valueAsString());
                }
            }

            calculationManager.runCalculation(cName, params, false);
        }
    }

}
