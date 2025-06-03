package edu.ucsc.neurobiology.vision.gui;

import java.awt.*;
import javax.swing.*;


/**
 * An status bar typically displayed at the bottom of the application window.
 * The status bar has an area for text messages. Additional areas can be added
 * to either the left or right sides to contain application specific status
 * areas.
 * And a progress bar which displays
 * the progress of extended operations, as well as a "stop" button that can be
 * used to abort certain prolonged actions.
 *
 * @author tonyj
 * Dumitru Petrusca, University of California, Santa Cruz
 */
public class StatusBar
    extends JPanel {
    private JLabel textLabel;
    private JProgressBar progressBar;


    /**
     * Create a new StatusBar
     */
    public StatusBar() {
        setLayout(new BorderLayout());

        textLabel = new JLabel();
        add(textLabel, BorderLayout.WEST);

        progressBar = new JProgressBar(JProgressBar.HORIZONTAL, 0, 100);
        progressBar.setPreferredSize(new Dimension(400, 5));
        add(progressBar, BorderLayout.EAST);
    }


    /**
     * Set the message to display in the status bar.
     * This message is thread safe and can be called from any thread.
     * @param message The message to display, or null to clear the message
     */
    public void setMessage(final String message) {
        if (!SwingUtilities.isEventDispatchThread()) {
            Runnable run = new Runnable() {
                public void run() {
                    setMessage(message);
                }
            };
            SwingUtilities.invokeLater(run);
        } else {
            if (message == null) {
                textLabel.setText(" ");
            } else {
                textLabel.setText(message);
            }
        }
    }


    public void setProgress(final int progress) {
        progressBar.setValue(progress);
    }
}
