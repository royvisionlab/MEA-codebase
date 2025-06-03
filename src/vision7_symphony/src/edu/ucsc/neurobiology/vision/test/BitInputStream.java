package edu.ucsc.neurobiology.vision.test;

import java.io.*;


public class BitInputStream {
    final protected static int MASK_SIZE = 8;
    final protected static int ZERO = 0;
    final protected static int ONES = ~0;
    final protected static int[] BIT_MASK = new int[MASK_SIZE];
    final protected static int[] FIELD_MASK = new int[MASK_SIZE];
    InputStream in;

    /**
     * This is a prefetched byte used to construct signed and unsigned
     * numbers with an arbitrary number of significant bits. */
    private int bits;

    /**
     * The number of valid bits remaining in the bits field. */
    private int validBits;

    // Generate the needed masks for various bit fields and for individual bits.
    static {
        int tempBit = 1;
        int tempField = 1;
        for (int i = 0; i < MASK_SIZE; i++) {
            // Set the masks.
            BIT_MASK[i] = tempBit;
            FIELD_MASK[i] = tempField;

            // Update the temporary values.
            tempBit <<= 1;
            tempField <<= 1;
            tempField++;
        }
    }


    public BitInputStream(InputStream in) {
        this.in = in;

        bits = 0;
        validBits = 0;
    }


    /**
     * A utility method to fetch the next byte in preparation for
     * constructing a bit field.  There is no protection for this
     * method; ensure that it is only called when a byte must be
     * fetched.
     */
    protected void fetchByte() throws IOException {
        bits = in.read();
        if (bits < 0) {
            throw new EOFException();
        }
        validBits = MASK_SIZE;
    }


    /**
     * A utility to force the next read to be byte-aligned. */
    public void byteAlign() {
        validBits = 0;
    }


    /**
     * Read a bit from the input stream and interpret this as a
     * boolean value.  A 1-bit is true; a 0-bit is false. */
    public boolean readBitFlag() throws IOException {
        if (validBits == 0) {
            fetchByte();
        }
        return ( (bits & BIT_MASK[--validBits]) != 0);
    }


    public int readBit() throws IOException {
        if (validBits == 0) {
            fetchByte();
        }
        return bits & BIT_MASK[--validBits];
    }


    /**
     * Read a signed value of n-bits from the input stream.
     */
//    public long readSBits(int n) throws IOException {
//        if (n == 0) {
//            return 0;
//        }
//        int value = (readBitFlag()) ? ONES : ZERO;
//        value <<= (--n);
//        return (value | readUBits(n));
//    }


    /**
     * Read a float value of n-bits from the stream. */
//    public float readFBits(int n) throws IOException {
//        if (n == 0) {
//            return 0.0f;
//        }
//        return ( (float) readSBits(n)) / 0x1000;
//    }


    /**
     * Read an unsigned value of n-bits from the input stream.
     */
//    public long readUBits(int n) throws IOException {
//        long value = ZERO;
//        while (n > 0) {
//            // Take the needed bits or the number which are valid
//            // whichever is less.
//            if (validBits == 0) {
//                fetchByte();
//            }
//            int nbits = (n > validBits) ? validBits : n;
//
//            // Take the bits and update the counters.
//            int temp = ( (bits >> (validBits - nbits)) & FIELD_MASK[nbits - 1]);
//            validBits -= nbits;
//            n -= nbits;
//
//            // Shift the value up to accomodate new bits.
//            value <<= nbits;
//            value |= temp;
//        }
//        return value;
//    }

    public void close() throws IOException {
        in.close();
    }
}
