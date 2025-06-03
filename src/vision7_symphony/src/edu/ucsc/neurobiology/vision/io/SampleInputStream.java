package edu.ucsc.neurobiology.vision.io;

import java.io.*;
import java.util.*;

import javax.swing.*;

import edu.ucsc.neurobiology.vision.*;


/**
    The <tt>SampleInputStream</tt> class reads the "raw" data from an underlying
    input stream (any instance of <tt>java.io.InputStream</tt>) and distributes the
    data to other classes interested in receiving the "sampled" data. It uses buffering
    techniques to make the reading as fast as possible. Any class that
    wants to receive "sampled" data needs to implements the <tt>RawDataListener</tt>
    interface. The interface contains two methods:
    <p>
        <b>processSample()</b> - called by <tt>SampleInputStream</tt> whenever a new
        sample is available.
        <br>
        <b>finishProcessing()</b> - called by <tt>SampleInputStream</tt> when the
        underlying <tt>InputStream</tt> reached it's end. Used to save the results of
        possible calculations and to free resources.
    </p>
    To connect a <tt>RawDataListener</tt> to a <tt>SampleInputStream</tt> one needs
    to call the later's <tt>addDataListener()</tt> method.
    <br>
    The underlying <tt>InputStream</tt> can be bound to a disk file or to a network
    socket (<tt>java.net.Socket</tt>).
    <br>
    <tt>SampleInputStream</tt> is actually a thread (extends <tt>java.lang.Thread
    </tt>) which will run asynchronously. The thread does not start right after the
    class is created. Instead it has to be started by an external class.
    This is because the classes that need "sampled" data may want to do
    preparatory work after <tt>SampleInputStream</tt> creation.
    <br>
    The <tt>run()</tt>
    method of <tt>SampleInputStream</tt> will run (if the thread is started) as long
    as there are bytes in the underlying input stream and it will stop executing
    (return) when there are no more data. While executing the thread does three
    things:
    <ul>
        <li> read byte chunks from the input stream into the internal buffer until
        there are enough bytes to form at least one full sample.
        <li> process the buffer assembling the "raw" bytes into samples (short[]).
        Every time the class has a full sample:
        <li> send the sample to all listeners by calling their <tt>processSample()</tt>
        method.
    </ul>

    These three operations are executed until the end of the input stream is reached.
    Then <tt>SampleInputStream</tt> calls the <tt> finishProcessing()</tt> method
    of all the registered listeners.

    <br><br>
    The data flow is done this way to allow for two or more classes to receive
    independently of each other and at the same time the "samples" of raw data.

    @author Dumitru Petrusca, University of California, Santa Cruz
 */

public final class SampleInputStream
    extends Thread {
    private InputStream input;
    private SampleListener sampleListener1, sampleListener2, sampleListener3;
    private final int nBufferSamples;
    private final int nElectrodes;
    private int samplesRead;

    /**
     * Creates a SampleInputStream which will read data from a specific imput stream.
     *
     * @param input The InputStream which will be used as a source of data
     * @param nElectrodes The number of electrodes the raw data are comming from
     * @param nBufferSamples The size of the internall buffer used for reading.
     * The size should be given in "number of samples"
     */
    public SampleInputStream(
        InputStream input, int nBufferSamples, int nElectrodes) throws IOException {

        super("SampleInputStream Thread");
        this.setPriority(Thread.MAX_PRIORITY);

        this.input = input;
        this.nBufferSamples = nBufferSamples;
        this.nElectrodes = nElectrodes;
    }


    /**
     * Fetches new data from the input stream and makes them available to listeners,
     * do NOT call.
     */
    public void run() {
        // the number of bytes read in the last read operation
        int lastBytesRead;

        // the number of bytes available at the moment in the buffer
        int bytesAvailable = 0;

        // the number of bytes assembled in samples and sent to listeners
        int bytesUsed = 0;

        // the number of samples available for sending to listeners
        int nAvailableSamples = 0;

        /// Not used (0)
        double scale, offset;

        short[] sample = new short[nElectrodes];
        int sampleSize = 2 * nElectrodes;
//        int sampleSize = 2 + (nElectrodes - 1)*3/2;

        // the byte buffer used to read data in
        byte[] buffer = new byte[nBufferSamples * nElectrodes * 2];

        int byteIndex, highByte, lowByte, i, sampleIndex;

        while (true) {
            // shift the remaining bytes to the begining of the buffer
            System.arraycopy(buffer, bytesUsed, buffer, 0, bytesAvailable - bytesUsed);
            bytesAvailable -= bytesUsed;

            // read data from the underlying stream until we have at least one full sample
            do {
                try {
                    lastBytesRead =
                        input.read(buffer, bytesAvailable, buffer.length - bytesAvailable);

                    // if we reached the end of the stream: 1) anounce the end of data;
                    // and 2) renurn from the run() method (the thread will dye).
                    if (lastBytesRead <= 0) {
                        endOfData();
                        return;
                    }

                    bytesAvailable += lastBytesRead;
                    nAvailableSamples = bytesAvailable / sampleSize;
                } catch (IOException e) {
                    // if an error happens: 1) anounce the end of data; 2) inform the
                    // user about the error and 3) return from the run() method
                    // (the thread will dye).
                    endOfData();
                    e.printStackTrace();
                    JOptionPane.showMessageDialog(
                        null, e, "IOException in SampleInputStream",
                        JOptionPane.ERROR_MESSAGE);
                    return;
                }
            } while (nAvailableSamples == 0);

            // calculate the number of bytes we are going to use now
            bytesUsed = nAvailableSamples * sampleSize;

            // send the data to SampleListeners
            byteIndex = 0;
            for (sampleIndex = 0; sampleIndex < nAvailableSamples; sampleIndex++) {
                for (i = 0; i < nElectrodes; i++) {
                    highByte = buffer[byteIndex++];
                    lowByte = buffer[byteIndex++];
                    sample[i] = (short) ( (highByte << 8) + (lowByte & 0xFF));
                }
                if (sampleListener1 != null) {
                    sampleListener1.processSample(sample);
                }
                if (sampleListener2 != null) {
                    sampleListener2.processSample(sample);
                }
                if (sampleListener3 != null) {
                    sampleListener3.processSample(sample);
                }
                samplesRead++;
            }

        } //while(true)
    }


    /*
        int readBytes = 0;
        private final int read(byte[] b, int offset, int length) throws IOException {
            int l = Math.min(length, 16000000*130 - readBytes);
            readBytes += l;

            if (l == 0)
                return -1;
            else
                return l;
        }
     */

    /**
     * Registers a new DataListener to this Input Stream. This function should be
     * called before the stream is started.
     */
    public void addDataListener(SampleListener listener) {
        if (sampleListener1 == null) {
            sampleListener1 = listener;
        } else if (sampleListener2 == null) {
            sampleListener2 = listener;
        } else if (sampleListener3 == null) {
            sampleListener3 = listener;
        } else {
            throw new IllegalArgumentException("Too many SampleListeners (3 allowed)");
        }
    }


    private final void endOfData() {
        if (sampleListener1 != null) {
            try {
                sampleListener1.finishSampleProcessing();
            } catch (IOException e) {
                Vision.reportException("Could not finish processing for " +
                                       sampleListener1.getClass().getName(), e);
            }
        }
        if (sampleListener2 != null) {
            try {
                sampleListener2.finishSampleProcessing();
            } catch (IOException e) {
                Vision.reportException("Could not finish processing for " +
                                       sampleListener2.getClass().getName(), e);
            }
        }
        if (sampleListener3 != null) {
            try {
                sampleListener3.finishSampleProcessing();
            } catch (IOException e) {
                Vision.reportException("Could not finish processing for " +
                                       sampleListener3.getClass().getName(), e);
            }
        }

        // close the stream
        try {
            input.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    public synchronized int getSamplesRead() {
        return samplesRead;
    }


    HashMap<String,String> info = new HashMap<String,String>();
    public synchronized HashMap<String,String> getStateInformation() {
        info.put("Samples Read", String.valueOf(samplesRead));
        return info;
    }
}
