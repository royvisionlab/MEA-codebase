package edu.ucsc.neurobiology.vision.gui;

import java.io.*;
import java.util.*;

import javax.swing.filechooser.FileFilter;


/**
 * This class implements a FileFilter which is a combination of other file
 * filters.  This is useful to create a FileFilter such as "All Java Files"
 * where the *.java and *.class files would be different FileFilters.
 *
 * @author Charles A. Loomis, University of California, Santa Cruz
 */
public class CombinedFileFilter
    extends FileFilter {

    private LinkedList<FileFilter> filters;
    private String description;


    /**
     * Create a new CombinedFileFilter with the given name.  Note that all of
     * the sub-FileFilters must be added with the add() method. */
    public CombinedFileFilter(String description) {
        this.description = description;
        filters = new LinkedList<FileFilter>();
    }


    /**
     * Add a new FileFilter to the list. */
    synchronized public void add(FileFilter filter) {
        if (filter != null && !filters.contains(filter)) {
            filters.add(filter);
        }
    }


    /**
     * Remove a FileFilter from the list. */
    synchronized public void remove(FileFilter filter) {
        if (filter != null) {
            filters.remove(filter);
        }
    }


    /**
     * Return true if the given filter is one of the sub-FileFilters. */
    synchronized public boolean contains(FileFilter filter) {
        return filters.contains(filter);
    }


    /**
     * Return whether or not the given file is accepted by any of the
     * sub-FileFilters. */
    synchronized public boolean accept(File file) {
        for (FileFilter f : filters) {
            if (f.accept(file)) {
                return true;
            }
        }
        return false;
    }


    /**
     * Return the description of this filter. */
    public String getDescription() {
        return description;
    }


    /**
     * Return the first FileFilter which will accept this file. */
    synchronized public FileFilter getAcceptingFileFilter(File file) {
        for (FileFilter f : filters) {
            if (f.accept(file)) {
                return f;
            }
        }
        return null;
    }

}
