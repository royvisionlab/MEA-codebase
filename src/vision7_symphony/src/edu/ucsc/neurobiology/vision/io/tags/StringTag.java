package edu.ucsc.neurobiology.vision.io.tags;

import java.io.*;


/**
 * A tag to encapsulate a string.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class StringTag
    extends ObjectEmbeddingTag {

    private String s;


    public StringTag() {
        super(NeuroTagSet.STRING_TAG);
    }


    public StringTag(String s) {
        super(NeuroTagSet.STRING_TAG);

        this.s = s;
    }


    public Tag read(int tagId, NeuroInputStream input, int len) throws IOException {
        return new StringTag(input.readString());
    }


    public void write(int tagId, NeuroOutputStream output) throws IOException {
        output.writeString(s);
    }


    public String toString() {
        return getName() + " " + s;
    }


    public Object getEmbeddedObject() {
        return s;
    }


    public void setEmbeddedObject(Object obj) {
        if (obj == null || obj instanceof String) {
            s = (String) obj;
        } else {
            throw new IllegalAccessError("The object should be String");
        }
    }
}
