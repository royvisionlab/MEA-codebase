package edu.ucsc.neurobiology.vision.io.tags;

// Copyright 2001, FreeHEP.
import java.io.*;


/**
 * Generic Tag to be used by NeuroInputStream/NeuroOutputStream.
 * The tag contains an ID, name and
 * a version. Concrete subclasses should implement the IO Read and Write methods.
 *
 * @author Mark Donszelmann
 * @author Charles Loomis
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public abstract class Tag {
    private int tagID;
    private String name;
    private int version;

    /**
     * This is the tagID for the default tag handler.
     */
    final public static int DEFAULT_TAG = -1;

    /**
     * Creates a new tag with the given ID and version.
     *
     * @param tagID int
     * @param version int
     */
    protected Tag(int tagID, int version) {
        this.tagID = tagID;
        this.version = version;
        this.name = null;
    }


    /**
     * Get the tag number.
     */
    public int getTag() {
        return tagID;
    }


    /**
     * Get the version number.
     */
    public int getVersion() {
        return version;
    }


    /**
     * Get the tag name. */
    public String getName() {
        if (name == null) {
            name = getClass().getName();
            int dot = name.lastIndexOf(".");
            name = (dot >= 0) ? name.substring(dot + 1) : name;
        }
        return name;
    }


    /**
     * This returns the type of block */
    public int getTagType() {
        return 0;
    }


    /**
     * This reads the information from the given input and returns a new Tag*/
    public abstract Tag read(int tagID, NeuroInputStream input, int len) throws
        IOException;


    /**
     * This writes the information to the given output */
    public abstract void write(int tagID, NeuroOutputStream output) throws IOException;


    public abstract String toString();
}
