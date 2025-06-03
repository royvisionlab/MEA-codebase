package edu.ucsc.neurobiology.vision.parameters;

import java.io.*;

import javax.swing.filechooser.FileFilter;


/**
 * This is a very simple filter used in the Vision file dialog. It accepts files
 * by extension and can also accept the selection of directories
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class SimpleFileFilter
    extends FileFilter {
    String extension, description;
    boolean acceptDirectories;


    public SimpleFileFilter(String extension, String description,
                            boolean acceptDirectories) {
        this.extension = extension;
        this.acceptDirectories = acceptDirectories;
        this.description = description;
    }


    public boolean accept(File f) {
        return f.isDirectory() || f.getName().endsWith(extension);
    }


    public String getDescription() {
        String s = description + " Files (." + extension + ")";
        if (acceptDirectories) {
            s += " Files & Directories";
        }
        return s;
    }


//    public String getDescription() {
//        return description;
//    }


    public String getExtension() {
        return extension;
    }


    public String toString() {
        return description;
    }
}
