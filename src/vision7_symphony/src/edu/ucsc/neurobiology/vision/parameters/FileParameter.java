package edu.ucsc.neurobiology.vision.parameters;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.filechooser.*;

import edu.ucsc.neurobiology.vision.gui.*;
import edu.ucsc.neurobiology.vision.io.*;


/**
 * This class represents a parameter that represents a String.
 * It can be shown on the screen using the ParametersTable class.
 * The value is displayed in a simple JTextField. When the string gets changed
 * a PropertyChangeEvent is send to all regitered listeners.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class FileParameter
extends AbstractParameter {

    private JTextField textField;
    private JPanel p;
    private final String initialValue;
    private String extension;
    private boolean acceptDirectories = false;

    public static final JFileChooser fileChooser;
    static {
        fileChooser = new JFileChooser(System.getProperty("user.dir"));
        fileChooser.removeChoosableFileFilter(fileChooser.getAcceptAllFileFilter());
        fileChooser.setApproveButtonText("Select");
        fileChooser.setDialogTitle("Select");
    }


    /**
     * Creates an instance of FileParameter which will display by default
     * the given string.
     *
     * @param name The internal name of the parameter
     * @param screenName The name which will be displayed on the screen
     *                   (can contain UNICODE characters)
     * @param toolTip The ToolTip associated with this parameter.
     * @param value The initial string displayed by this parameter
     */
    public FileParameter(
            String name, String screenName, String toolTip, String value,
            final String extension) {

        super(name, screenName, toolTip);
        this.initialValue = value;
        this.extension = extension;

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

        if (extension.startsWith("+")) {
            this.extension = extension.substring(1, extension.length());
            acceptDirectories = true;
        } else {
            acceptDirectories = false;
        }

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
                    textField.setText(fileChooser.getSelectedFile().getAbsolutePath());
                }
            }
        });
        p.add(b, BorderLayout.EAST);

        if (IOUtil.isValidFile(initialValue) || IOUtil.isValidFolder(initialValue)) {
            textField.setText(new java.io.File(initialValue).getAbsolutePath());
        } else {
            textField.setText(initialValue);
        }
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
        String name = textField.getText().trim();

        if (IOUtil.isValidNetworkAdress(name)) {
            return true;
        }

        if (IOUtil.isValidFile(name) && name.endsWith(extension)) {
            return true;
        }

        if (acceptDirectories && IOUtil.isValidFolder(name)) {
            return true;
        }
        
        if(IOUtil.isValidDataset(name)) {
            return true;
        }

        return false;
    }


//	public boolean isValid() {
//	String fileName = textField.getText().trim();
//	File file = new File(fileName);

//	if (!file.exists()) {
//	return false;
//	}

//	if (file.isFile() && file.canRead() && file.getName().endsWith(extension)) {
//	return true;
//	} else if (acceptDirectories && file.isDirectory()) {
//	return true;
//	} else {
//	return false;
//	}
//	}


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
