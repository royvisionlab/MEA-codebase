package edu.ucsc.neurobiology.vision.parameters;

import java.io.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.filechooser.FileFilter;

import edu.ucsc.neurobiology.vision.gui.*;
import static edu.ucsc.neurobiology.vision.parameters.FileParameter.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * This class represents a parameter that represents a String.
 * It can be shown on the screen using the ParametersTable class.
 * The value is displayed in a simple JTextField. When the string gets changed
 * a PropertyChangeEvent is send to all regitered listeners.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class FileListParameter
    extends AbstractParameter {

    private JTextField textField;
    private JPanel p;
    private final String initialValue;
    private String extension;
//    private JFileChooser fileChooser;
    private boolean acceptDirectories = false;


    /**
     * Creates an instance of FileListParameter which will display by default
     * the given string.
     *
     * @param name The internal name of the parameter
     * @param screenName The name which will be displayed on the screen
     *                   (can contain UNICODE characters)
     * @param toolTip The ToolTip associated with this parameter.
     * @param initialValue The initial string displayed by this parameter
     */
    public FileListParameter(
        String name, String screenName, String toolTip, String initialValue,
        final String extension) {

        super(name, screenName, toolTip);
        this.initialValue = initialValue;

//        fileChooser = new JFileChooser(System.getProperty("user.dir"));
//        fileChooser.removeChoosableFileFilter(fileChooser.getAcceptAllFileFilter());
//        fileChooser.setApproveButtonText("Select");
//        fileChooser.setDialogTitle("Select");
        if (extension.startsWith("+")) {
            this.extension = extension.substring(1, extension.length());
//            fileChooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
            acceptDirectories = true;
        } else {
            this.extension = extension;
            acceptDirectories = false;
//            fileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
        }

//        fileChooser.addChoosableFileFilter(new SimpleFileFilter(extension, acceptDirectories));

        initialize();
    }


    /**
     * Retrives the text string of this parameter.
     * a legal double value.
     *
     * @return The double value or NaN
     */
    public String getValue() {
        String s = textField.getText().trim();
        return s;
    }


    protected void initialize() {
        super.initialize();
        p = new JPanel(new BorderLayout());

        textField = new JTextField();
        textField.setEditable(true);
        textField.setToolTipText(toolTip);
        textField.getDocument().addDocumentListener(new DocumentListener() {
            public void insertUpdate(DocumentEvent e) {
                applyInput(); //System.out.println("insert");
            }


            public void removeUpdate(DocumentEvent e) {
                applyInput(); //System.out.println("remove");
            }


            public void changedUpdate(DocumentEvent e) {
                applyInput(); //System.out.println("change");
            }
        });
        p.add(textField, BorderLayout.CENTER);

        JButton b = new JButton(VisionImage.getIcon("open.gif"));
        b.setMargin(new Insets(0, 0, 0, 0));
        b.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                if (acceptDirectories) {
                    fileChooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
                } else {
                    fileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
                }
                FileFilter[] filter = fileChooser.getChoosableFileFilters();
                for (int i = 0; i < filter.length; i++) {
                    fileChooser.removeChoosableFileFilter(filter[i]);
                }
                fileChooser.addChoosableFileFilter(
                    new SimpleFileFilter(extension, "", acceptDirectories));

                Window w = (Window) SwingUtilities.getAncestorOfClass(Window.class,
                    textField);
                if (fileChooser.showOpenDialog(w) == JFileChooser.APPROVE_OPTION) {
                    textField.setText(textField.getText() + ";" +
                                      fileChooser.getSelectedFile().getAbsolutePath());
                }
            }
        });
        p.add(b, BorderLayout.EAST);

        textField.setText(initialValue);
    }


    private void applyInput() {
        if (isValid()) {
            fireParameterChange();
            textField.setForeground(normalForeground);
            textField.setBackground(normalBackground);
        } else {
            textField.setForeground(errorForeground);
            textField.setBackground(errorBackground);
        }
    }


    /**
     * Returns true if the specified number is acceptable for this parameter
     * and false otherwise.
     */
    public boolean isValid() {
        String[] files = StringUtil.decomposeString(textField.getText().trim(), ";");

        boolean areFiles = true, areFolders = true, areAdresses = true;
        for (int i = 0; i < files.length; i++) {
            if (!files[i].startsWith("net://")) {
                areAdresses = false;
            }

            File f = new File(files[i]);
            if (! (f.exists() && f.isFile() && f.canRead() &&
                   f.getName().endsWith(extension))) {
                areFiles = false;
            }

            if (! (f.exists() && f.isDirectory())) {
                areFolders = false;
            }
        }

        if (areFiles) {
            return true;
        } else if (acceptDirectories && areFolders) {
            return true;
        } else if (areAdresses) {
            return true;
        }

        return false;
    }


    public JComponent getValueComponent() {
        return p;
    }


    public String valueAsString() {
        return getValue();
    }


    public Object valueAsObject() {
        return getValue();
    }
}
