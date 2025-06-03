package edu.ucsc.neurobiology.vision.util;

import java.io.*;


/**
 * This class represents a growing array of primitive double values.
 * It is used to gain performance over an ArrayList or a LinkedList
 * (don't have to create wrappers for each number).
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class DoubleList
    implements Serializable {

    protected int initialSize;
    protected double incrementFactor;
    protected int counter;
    protected double[] buffer;


    /**
     * Creates a new DoubleList with a initial size of 256 and an increment factor of 1.5.
     */
    public DoubleList() {
        initialize(256, 1.5);
    }


    /**
     * Creates a new DoubleList with a certain initial size and the default increment.
     *
     * @param initialSize The initial size the array will have
     */
    public DoubleList(int initialSize) {
        initialize(initialSize, 1.5);
    }


    /**
     * Creates a new DoubleList with a certain initial size and increment.
     *
     * @param initialSize The initial size the array will have
     * @param incrementFactor The increment used to grow the array when it fills
     */
    public DoubleList(int initialSize, int incrementFactor) {
        initialize(initialSize, incrementFactor);
    }


    protected void initialize(int initialSize, double incrementFactor) {
        this.initialSize = initialSize;
        this.incrementFactor = incrementFactor;
        this.counter = 0;
        buffer = new double[initialSize];
    }


    protected void ensureSize(int size) {
        if (buffer.length >= size) {
            return;
        }

        double[] newBuffer = new double[size];
        System.arraycopy(buffer, 0, newBuffer, 0, buffer.length);
        buffer = newBuffer;
    }


    /**
     * Used to put a new double value into the array.
     *
     * @param d The new value to be appended
     */
    public void add(double d) {
        if (counter >= buffer.length) {
            ensureSize( (int) (buffer.length * incrementFactor));

        }
        buffer[counter] = d;
        counter++;
    }


    public double get(int i) {
        if (i < 0 || i > counter) {
            throw new IndexOutOfBoundsException();
        }

        return buffer[i];
    }


    public void set(int i, double value) {
        if (i < 0 || i > counter) {
            throw new IndexOutOfBoundsException();
        }

        buffer[i] = value;
    }


    public int size() {
        return counter;
    }


    public void clear() {
        counter = 0;
    }


    public double[] toArray() {
        double[] d = new double[counter];
        System.arraycopy(buffer, 0, d, 0, size());
        return d;
    }


    public double getMin() {
        double min = Double.POSITIVE_INFINITY;
        for (int i = 0; i < counter; i++) {
            if (buffer[i] < min) {
                min = buffer[i];
            }
        }
        return min;
    }


    public double getMax() {
        double max = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < counter; i++) {
            if (buffer[i] > max) {
                max = buffer[i];
            }
        }
        return max;
    }


    public String toString() {
        return "DoubleList[" + size() + "]";
    }


    public boolean contains(double value) {
        for (int i = 0; i < counter; i++) {
            if (buffer[i] == value) {
                return true;
            }
        }
        return false;
    }

}
