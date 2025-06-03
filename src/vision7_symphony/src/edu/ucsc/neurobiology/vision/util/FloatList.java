package edu.ucsc.neurobiology.vision.util;


/**
 * This class represents a growing array of primitive float values.
 * It is used to gain performance over an ArrayList or a LinkedList
 * (don't have to create wrappers for each number).
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class FloatList {
    protected int initialSize;
    protected double incrementFactor;
    protected int counter;
    protected float[] buffer;


    /**
     * Creates a new DoubleList with a initial size of 256 and an increment of 256.
     */
    public FloatList() {
        initialize(256, 1.5);
    }


    /**
     * Creates a new DoubleList with a certain initial size and the default increment.
     *
     * @param initialSize The initial size the array will have
     */
    public FloatList(int initialSize) {
        initialize(initialSize, 1.5);
    }


    /**
     * Creates a new DoubleList with a certain initial size and increment.
     *
     * @param initialSize The initial size the array will have
     * @param incrementFactor The increment used to grow the array when it fills
     */
    public FloatList(int initialSize, int incrementFactor) {
        initialize(initialSize, incrementFactor);
    }


    protected void initialize(int initialSize, double incrementFactor) {
        this.initialSize = initialSize;
        this.incrementFactor = incrementFactor;
        this.counter = 0;
        buffer = new float[initialSize];
    }


    protected void ensureSize(int size) {
        if (buffer.length >= size) {
            return;
        }

        float[] newBuffer = new float[size];
        System.arraycopy(buffer, 0, newBuffer, 0, buffer.length);
        buffer = newBuffer;
    }


    /**
     * Used to put a new int value into the array.
     *
     * @param d The new value to be appended
     */
    public final void add(float d) {
        if (counter >= buffer.length) {
            ensureSize( (int) (buffer.length * incrementFactor));
        }
        buffer[counter] = d;
        counter++;
    }


    public float get(int i) {
        if (i < 0 || i > counter) {
            throw new IndexOutOfBoundsException();
        }

        return buffer[i];
    }


    public void set(int i, int value) {
        if (i < 0 || i > counter) {
            throw new IndexOutOfBoundsException();
        }

        buffer[i] = value;
    }


    public int size() {
        return counter;
    }


    public boolean isEmpty() {
        return counter == 0;
    }


    public void clear() {
        counter = 0;
    }


    public float[] toArray() {
        float[] d = new float[counter];
        System.arraycopy(buffer, 0, d, 0, size());
        return d;
    }


    public float getMin() {
        float min = Integer.MAX_VALUE;
        for (int i = 0; i < counter; i++) {
            if (buffer[i] < min) {
                min = buffer[i];
            }
        }
        return min;
    }


    public float getMax() {
        float max = Integer.MAX_VALUE;
        for (int i = 0; i < counter; i++) {
            if (buffer[i] > max) {
                max = buffer[i];
            }
        }
        return max;
    }


    public boolean contains(int x) {
        for (int i = 0; i < counter; i++) {
            if (buffer[i] == x) {
                return true;
            }
        }

        return false;
    }


    public String toString() {
        return "DoubleList[" + size() + "]";
    }

}
