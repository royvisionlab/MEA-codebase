package edu.ucsc.neurobiology.vision.io.tags;

import java.io.*;


/**
 * A tag to encapsulate a double.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class DoubleTag
    extends ObjectEmbeddingTag {

    private double value;


    public DoubleTag() {
        super(NeuroTagSet.DOUBLE_TAG);
    }


    public DoubleTag(double value) {
        super(NeuroTagSet.DOUBLE_TAG);
        this.value = value;
    }


    public double doubleValue() {
        return value;
    }


    public int intValue() {
        return (int) value;
    }


    public Tag read(int tagId, NeuroInputStream input, int len) throws IOException {
        return new DoubleTag(input.readDouble());
    }


    public void write(int tagId, NeuroOutputStream output) throws IOException {
        output.writeDouble(value);
    }


    public String toString() {
        return getName() + " " /* + " " + value*/;
    }


    public Object getEmbeddedObject() {
        return new Double(value);
    }


    public void setEmbeddedObject(Object obj) {
        if (obj == null) {
            value = Double.NaN;
        } else if (obj instanceof Double) {
            value = ( (Double) obj).doubleValue();
        } else {
            throw new IllegalAccessError("The object should be Double");
        }
    }
}
