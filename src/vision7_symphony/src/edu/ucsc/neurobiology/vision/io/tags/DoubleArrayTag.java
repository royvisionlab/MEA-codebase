package edu.ucsc.neurobiology.vision.io.tags;

import java.io.*;


/**
 * A tag to encapsulate an array of doubles.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class DoubleArrayTag
    extends ObjectEmbeddingTag {

    private double[] value;


    public DoubleArrayTag() {
        super(NeuroTagSet.DOUBLE_ARRAY_TAG);
    }


    public DoubleArrayTag(double[] value) {
        super(NeuroTagSet.DOUBLE_ARRAY_TAG);
        this.value = cloneArray(value);
    }


    public final double[] cloneArray(double[] a) {
        if (a == null) {
            return null;
        }
        double[] b = new double[a.length];
        System.arraycopy(a, 0, b, 0, a.length);
        return b;
    }


    public double[] getArray() {
        return cloneArray(value);
    }


    public Tag read(int tagId, NeuroInputStream input, int len) throws IOException {
        int n = input.readInt();
        if (n == 0) {
            return new DoubleArrayTag(null);
        } else {
            double[] v = new double[n];
            for (int i = 0; i < n; i++) {
                v[i] = input.readDouble();
            }
            return new DoubleArrayTag(v);
        }
    }


    public void write(int tagId, NeuroOutputStream output) throws IOException {
        if (value == null) {
            output.writeInt(0);
        } else {
            output.writeInt(value.length);
            for (int i = 0; i < value.length; i++) {
                output.writeDouble(value[i]);
            }
        }
    }


    public String toString() {
        return getName() + " ";
    }


    public Object getEmbeddedObject() {
        return cloneArray(value);
    }


    public void setEmbeddedObject(Object obj) {
        if (obj == null || (obj instanceof double[])) {
            value = cloneArray( (double[]) obj);
        } else {
            throw new IllegalAccessError("The object should be double[]");
        }
    }
}
