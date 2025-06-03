package edu.ucsc.neurobiology.vision.io;

import java.io.*;

import javax.swing.*;

import edu.ucsc.neurobiology.vision.*;


/**
    This class reads the "raw" data from an underlying
    InputStream and distributes the
    data to classes interested in receiving the "sampled" data. It uses buffering
    techniques to make the reading as fast as possible. Any class that
    wants to receive "sampled" data needs to implements the <tt>SampleListener</tt>
    interface. To connect a <tt>SampleListener</tt> to this class one needs
    to call the <tt>addSampleListener()</tt> method.
    <br>
    The underlying <tt>InputStream</tt> can be bound to a disk file, a disk directory
 * (DirectoryInputStream) or to a network socket (<tt>java.net.Socket</tt>).
    <br>
    <tt>CompressedSampleInputStream</tt> is actually a thread which will run
 * asynchronously. The thread does not start right after the
    class is created. Instead it has to be started by an external call.
    This is because the classes that need "sampled" data may want to do
    preparatory work after the <tt>CompressedSampleInputStream</tt> creation.
    <br>
    The <tt>run()</tt>
    method of <tt>SampleInputStream</tt> will run (if the thread is started) as long
    as there are bytes in the underlying input stream and it will stop executing
    (and return) when there are no more bytes. While executing, the thread does three
    things:
    <ul>
        <li> read byte chunks from the input stream into the internal buffer until
        there are enough bytes to form at least one full sample.
        <li> process the buffer assembling the "raw" bytes into samples (a short[]).
        Every time the class has a full sample:
        <li> send the sample to all listeners by calling their <tt>processSample()</tt>
        method.
    </ul>
    These three operations are executed until the end of the input stream is reached.
 Then <tt>CompressedSampleInputStream</tt> calls the <tt> finishSampleProcessing()</tt>
 * method of all the registered listeners.
    <br><br>
    The data flow is done this way to allow for two or more classes to receive samples
    independently of each other.

 <br>
 <br>
 <b>The listeners should not modify the short[] containing the sample !!!</b>
 <br>
 <br>

    @author Dumitru Petrusca, University of California, Santa Cruz
 */
public final class CompressedSampleInputStream
extends Thread {

    private InputStream input;
    private SampleListener sampleListener1, sampleListener2, sampleListener3;
    private final int nElectrodes;
    private int samplesRead;
    private final int bufferSizeInBytes;

    private boolean initializing = true;
    private final short[][] sampleBuffer;

    private boolean initialize = false;
    private SampleInitializationListener initializationListener;

    private long startSample, stopSample;


    /**
     * Creates a CompressedSampleInputStream which will read data from a specific input stream.
     *
     * @param input The InputStream which will be used as a source of data
     * @param bufferSizeInBytes The size of the internal buffer used for reading.
     * @param nElectrodes The number of electrodes the raw data are coming from
     * @param initialize Whether initialization is needed
     * The size should be given in "number of samples"
     */
    public CompressedSampleInputStream(
            InputStream input, int bufferSizeInBytes, int nElectrodes, boolean initialize, long startSample, long stopSample) throws
            IOException {

        super("CompressedSampleInputStream Thread");
        this.setPriority(Thread.MAX_PRIORITY);

        this.initialize = initialize;
        this.input = input;
        this.bufferSizeInBytes = bufferSizeInBytes;
        this.nElectrodes = nElectrodes;
        this.sampleBuffer = new short[1000 * 20][nElectrodes]; // 1 second buffer
        this.startSample = startSample;
        this.stopSample = stopSample;
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

        short[] sample = new short[nElectrodes];

        //If nElectrodes is odd, treat electrode zero differently.
        final int sampleSize;
        if (nElectrodes % 2 == 1) {
            sampleSize = 2 + (nElectrodes - 1) * 3 / 2;
        } else {
            sampleSize = (nElectrodes) * 3 / 2;
        }

        // the byte buffer used to read data in
        byte[] buffer = new byte[bufferSizeInBytes];

        int byteIndex, sampleIndex;


        while (true) {

            // shift the remaining bytes to the beginning of the buffer
            System.arraycopy(buffer, bytesUsed, buffer, 0, bytesAvailable - bytesUsed);
            bytesAvailable -= bytesUsed;

            // read data from the underlying stream until we have at least one full sample
            do {
                try {

                    lastBytesRead =
                        input.read(buffer, bytesAvailable, buffer.length - bytesAvailable);

                    // if we reached the end of the stream: 1) announce the end of data;
                    // and 2) return from the run() method (the thread will die).
                    // The thread ends cleanly reading up to the last COMPLETE sample.
                    if (lastBytesRead <= 0) {
                        endOfData();
                        return;
                    }

                    bytesAvailable += lastBytesRead;
                    nAvailableSamples = bytesAvailable / sampleSize;
                } catch (IOException e) {
                    // if an error happens: 1) announce the end of data; 2) inform the
                    // user about the error and 3) return from the run() method
                    // (the thread will die).

                    endOfData();
                    e.printStackTrace();
                    JOptionPane.showMessageDialog(
                            null, e, "IOException in CompressedSampleInputStream",
                            JOptionPane.ERROR_MESSAGE);
                    return;
                }
            
            } while (nAvailableSamples == 0);

            // calculate the number of bytes we are going to use now
            bytesUsed = nAvailableSamples * sampleSize;

            // send the data to SampleListeners
            byteIndex = 0;
            int electrode = 0, b1, b2, b3;
            for (sampleIndex = 0; sampleIndex < nAvailableSamples; sampleIndex++) {
                if (nElectrodes % 2 == 1) {
                    // the TTL is not encoded
                    sample[0] = (short) ( (buffer[byteIndex++] << 8) +
                            (buffer[byteIndex++] & 0xff));
                    if (samplesRead >= startSample) {
                        for (electrode = 1; electrode < nElectrodes; ) {
                            b1 = buffer[byteIndex++] & 0xff;
                            b2 = buffer[byteIndex++] & 0xff;
                            b3 = buffer[byteIndex++] & 0xff;

                            sample[electrode++] = (short) ( (b1 << 4) + (b2 >> 4) - 2048);
                            sample[electrode++] = (short) ( ( (b2 & 0x0f) << 8) + b3 - 2048);
                        }
                    } else {
                        //Do not preserve, non-TTL data before start sample
                        for (electrode = 1; electrode < nElectrodes;) {
                            byteIndex+=3;
                            sample[electrode++] = 0;
                            sample[electrode++] = 0;
                        }
                    }
                } else {
                    // the TTL is encoded
                    if (samplesRead >= startSample) {

                        for (electrode = 0; electrode < nElectrodes; ) {

                            b1 = buffer[byteIndex++] & 0xff;
                            b2 = buffer[byteIndex++] & 0xff;
                            b3 = buffer[byteIndex++] & 0xff;

                            sample[electrode++] = (short) ( (b1 << 4) + (b2 >> 4) - 2048);
                            sample[electrode++] = (short) ( ( (b2 & 0x0f) << 8) + b3 - 2048);

                        }
                    } else {
                        for (electrode = 0; electrode < nElectrodes; ) {
                            if (electrode == 0) {
                                b1 = buffer[byteIndex++] & 0xff;
                                b2 = buffer[byteIndex++] & 0xff;
                                b3 = buffer[byteIndex++] & 0xff;

                                sample[electrode++] = (short) ( (b1 << 4) + (b2 >> 4) - 2048);
                                sample[electrode++] = 0;
                            } else {
                                byteIndex+=3;
                                sample[electrode++] = 0;
                                sample[electrode++] = 0;
                            }
                        }

                    }
                }

                if (initialize && initializing) {
                    if (samplesRead < sampleBuffer.length) {
                        // keep a copy of the sample
                        System.arraycopy(
                                sample, 0, sampleBuffer[samplesRead], 0, nElectrodes);
                    } else {
                        // we filled the buffer sample, send it to the listeners
                        if (initializationListener != null) {
                            initializationListener.initialize(sampleBuffer);
                        }

                        // send all the samples in the buffer to the listeners
                        for (int i = 0; i < sampleBuffer.length; i++) {
                            sendSample(sampleBuffer[i]);
                        }

                        // send also the current sample which will be otherwise lost
                        sendSample(sample);

                        initializing = false;
                    }
                } else {
                    sendSample(sample);
                }

                samplesRead++;

                if (samplesRead == stopSample) {
                    //this case occurs when the total number of samples is very small.
                    //functions that require initialization (a full buffer of data)
                    //will fail and throw exceptions.  Functions that do not will succeed.
                    if (samplesRead <= sampleBuffer.length) {
                        for (int i = 0; i < samplesRead; i++) {
                            sendSample(sampleBuffer[i]);
                        }
                    }
            
                    endOfData();
                    
                    return;
                }
            }
        } //while(true)
        
        
    }


    public void sendSample(short[] sample) {
        if (sampleListener1 != null) {
            sampleListener1.processSample(sample);
        }
        if (sampleListener2 != null) {
            sampleListener2.processSample(sample);
        }
        if (sampleListener3 != null) {
            sampleListener3.processSample(sample);
        }
    }


    /**
     * Registers a new DataListener to this Input Stream. This function should be
     * called before the stream is started.
     */
    public void addSampleListener(SampleListener listener) {
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


    public void setInitializationListener(SampleInitializationListener listener) {
        if (initializationListener != null) {
            throw new IllegalArgumentException(
            "The Initialization Listener was already set");
        }

        this.initializationListener = listener;
    }


    private final void endOfData() {
//		System.err.println(samplesRead);
//		System.err.println(samplesSent);

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

}
