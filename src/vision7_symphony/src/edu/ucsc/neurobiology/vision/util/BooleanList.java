package edu.ucsc.neurobiology.vision.util;


/**
 * This class represents a growing array of primitive boolean values.
 * It is used to gain performance over an ArrayList or a LinkedList
 * (don't have to create wrappers for each number).
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class BooleanList {
    protected int initialSize;
    protected double incrementFactor;
    protected int counter;
    protected boolean[] buffer;


    /**
     * Creates a new DoubleList with a initial size of 256 and an increment of 256.
     */
    public BooleanList() {
        initialize(256, 1.5);
    }


    /**
     * Creates a new DoubleList with a certain initial size and the default increment.
     *
     * @param initialSize The initial size the array will have
     */
    public BooleanList(int initialSize) {
        initialize(initialSize, 1.5);
    }


    /**
     * Creates a new DoubleList with a certain initial size and increment.
     *
     * @param initialSize The initial size the array will have
     * @param incrementFactor The increment used to grow the array when it fills
     */
    public BooleanList(int initialSize, int incrementFactor) {
        initialize(initialSize, incrementFactor);
    }


    protected void initialize(int initialSize, double incrementFactor) {
        this.initialSize = initialSize;
        this.incrementFactor = incrementFactor;
        this.counter = 0;
        buffer = new boolean[initialSize];
    }


    protected void ensureSize(int size) {
        if (buffer.length >= size) {
            return;
        }

        boolean[] newBuffer = new boolean[size];
        System.arraycopy(buffer, 0, newBuffer, 0, buffer.length);
        buffer = newBuffer;
    }


    /**
     * Used to put a new boolean value into the array.
     *
     * @param d The new value to be appended
     */
    public void add(boolean d) {
        if (counter >= buffer.length) {
            ensureSize( (int) (buffer.length * incrementFactor));
        }

        buffer[counter] = d;
        counter++;
    }


    public boolean get(int i) {
        if (i < 0 || i > counter) {
            throw new IndexOutOfBoundsException();
        }

        return buffer[i];
    }


    public void set(int i, boolean value) {
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


    public boolean[] toArray() {
        boolean[] d = new boolean[counter];
        System.arraycopy(buffer, 0, d, 0, size());
        return d;
    }

}
