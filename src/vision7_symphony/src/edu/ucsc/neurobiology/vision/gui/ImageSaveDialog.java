package edu.ucsc.neurobiology.vision.gui;

import java.beans.*;

import java.awt.*;
import javax.swing.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.parameters.*;


/**
 * @author not attributable
 */
public class ImageSaveDialog
    extends JFileChooser {

    public ImageSaveDialog(String saveFolder, Config configuration) {
        ParametersTable _table = configuration.getParameterGroup("Image Parameters");
        _table.getColumn("Value").setMinWidth(0);
        _table.getColumn("Value").setMaxWidth(100000);
        _table.getColumn("Value").setWidth(70);
        final JScrollPane epsPanel = new JScrollPane(_table);
        epsPanel.setPreferredSize(new Dimension(200, 80));
        epsPanel.setSize(new Dimension(200, 80));


        final JPanel pngPanel = new JPanel();

        addPropertyChangeListener(new PropertyChangeListener() {
            public void propertyChange(PropertyChangeEvent e) {
                if (e.getPropertyName().equals(JFileChooser.
                                               FILE_FILTER_CHANGED_PROPERTY)) {
                    SimpleFileFilter f = ( (SimpleFileFilter) e.getNewValue());
                    if (f == null) {
                        return;
                    }

                    if (f.getExtension().equals("eps")) {
                        setAccessory(epsPanel);
                    } else if (f.getExtension().equals("png")) {
                        setAccessory(pngPanel);
                    }

                    validate();
                    repaint();
                }
            }
        });
        removeChoosableFileFilter(getAcceptAllFileFilter());
        SimpleFileFilter epsFilter = new SimpleFileFilter("eps",
            "Encapsulated PostScript", false);
        SimpleFileFilter pngFilter = new SimpleFileFilter("png",
            "Portable Network Graphics", false);
        addChoosableFileFilter(epsFilter);
        addChoosableFileFilter(pngFilter);
        setFileFilter(pngFilter);
    }
}
