package edu.ucsc.neurobiology.vision.test;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 *         Matthew Grivich, University of California, Santa Cruz
 */
public class IntegerCircularBuffer {
    public static final int LOSS = 0;
    public static final int NOLOSS = 1;
    int type;
    protected int[] data;
    protected final int n;
    protected int head;
    protected int tail;


    public IntegerCircularBuffer(int type, int capacity) {
        this.type = type;
        this.n = capacity + 1;
        this.data = new int[n];
        for (int i = 0; i < data.length; i++) {
            data[i] = Integer.MIN_VALUE;
        }
        this.head = 0;
        this.tail = 0;
    }


    public final void pushLast(int value) {
        if ( (head - tail + n) % n == n - 1) {
            if (type == LOSS) {
                data[tail] = Integer.MIN_VALUE;
                tail++;
                tail %= n;
            } else {
                throw new IllegalStateException("buffer is full");
            }
        }

        data[head] = value;
        head++;
        head %= n;
    }


    /*
        public final int popLast() {
            //        if (head == tail) throw new IllegalStateException("no data");
            head--;
            if (head < 0) {
                head += n;
            }
            return data[head];
        }
     */

    public final int popFirst() {
        if (head == tail) {
            throw new IllegalStateException("no data");
        }

        int result = data[tail];
//        data[tail] = Integer.MIN_VALUE;
        tail++;
        tail %= n;

        return result;
    }


    public final int inspect(int i) {
        int index = (tail + i) % n;
        return data[index];
    }


    /*
        public final int inspectLast() {
            if (head == tail) {
                throw new IllegalStateException("no data");
            }
            return data[ (head - 1 + n) % n];
        }
     */

    public final int inspectFirst() {
        if (head == tail) {
            throw new IllegalStateException("no data");
        }
        return data[tail];
    }


    public final void dump() {
        System.out.println("===========================");
        System.out.println("tail: " + tail);
        System.out.println("head: " + head);
        for (int i = 0; i < data.length; i++) {
            System.out.print(data[i] + ", ");
        }
        System.out.println();
        System.out.println("===========================");
    }


    public final int size() {
        int size = head - tail;
        if (size < 0) {
            size += n;
        }
        return size;
    }


    public final boolean isFull() {
        return (head - tail + n) % n == (n - 1);
    }


    public final boolean isEmpty() {
        return head == tail;
    }


    public final int capacity() {
        return n - 1;
    }


    public String toString() {
        String s = "head=" + head + ", tail=" + tail + ", size=" + size();
        if (isEmpty()) {
            s += ", empty";
        }
        if (isFull()) {
            s += ", full";
        }

        s += "\n";
        for (int i = 0; i < n; i++) {
            s += data[i] + ", ";
        }
        s += "\n";
        for (int i = 0; i < size(); i++) {
            s += inspect(i) + ", ";
        }

        return s;
    }

}
