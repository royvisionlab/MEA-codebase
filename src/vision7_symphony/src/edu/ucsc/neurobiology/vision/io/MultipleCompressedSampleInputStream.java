package edu.ucsc.neurobiology.vision.io;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import edu.ucsc.neurobiology.vision.Vision;


/**
 * This class can open connections to one or more raw data sources, creates a
 * CompressedSampleInputStream object for each source and then switches between the sources
 * when required to give the impression to the rest of the IO pipeline that it is working 
 * with only one raw data source. The contract of this class is the same as for
 * CompressedSampleInputStream.
 *
 * @see CompressedSampleInputStream
 * @author Dumitru Petrusca, University of California, Santa Cruz
 * @author Matthew Grivich, The Salk Institute
 */
public final class MultipleCompressedSampleInputStream
implements SampleListener, SampleInitializationListener {

    private static final boolean DEBUG_FLOW = false;

    private AsynchronousInputStream asynchronousInput;
    private String[] rawDataSource;
    private final boolean waitForData;

    //Input stream for primary data.
    private InputStream[] inputStream;
    private long[] startSamples, stopSamples;
    //Used so that the ready signal can be sent back to the daq computer
    private OutputStream[] outputStream;
    private RawDataHeader512[] header;

    private int currentInputStream = -1;
    private boolean finished = false;
    private int bufferSizeInBytes;
    private int nBuffers;
    private int nElectrodes;


    private SampleListener
    sampleListener1, sampleListener2, sampleListener3, sampleListener4, sampleListener5;
    private int samplesRead;
    private SampleInitializationListener initializationListener;
    CompressedSampleInputStream sampleInputStream;


    public MultipleCompressedSampleInputStream(String rawDataSources) throws IOException {
        // buffer size: 500 samples
        // nBuffers: 100
        this(rawDataSources, 500 * 770, 100); // 35 Mb
    }

    public MultipleCompressedSampleInputStream(String rawDataSources, int bufferSizeInBytes, int nBuffers) throws IOException {
        this(rawDataSources, bufferSizeInBytes, nBuffers, false);
    }

    /**
     * Creates a CompressedSampleInputStream which will read data from a specific
     * input stream.
     *
     * @param rawDataSources String a string of file/folder names.  See DataFileStringParser for syntax.
     * @param bufferSizeInBytes int The size of the internal buffer used for reading.
     * @param nBuffers int how many such buffers
     * @param waitForData boolean Whether the input stream should take some time to wait for more data 
     *     to come in instead of exiting as soon as EOF is reached.  Only makes sense for file inputs.
     * @throws IOException
     */
    public MultipleCompressedSampleInputStream(String rawDataSources, int bufferSizeInBytes, 
            int nBuffers, boolean waitForData) throws IOException {

        this.waitForData = waitForData;
        
        DataFileStringParser parser = new DataFileStringParser(rawDataSources);
        this.rawDataSource = parser.getDatasets();		
        double[] startTimes = parser.getStartTimes();
        double[] stopTimes = parser.getStopTimes();

        this.bufferSizeInBytes = bufferSizeInBytes;
        this.nBuffers = nBuffers;

        header = new RawDataHeader512[rawDataSource.length];
        startSamples = new long[rawDataSource.length];
        stopSamples = new long[rawDataSource.length];
        inputStream = new InputStream[rawDataSource.length];
        outputStream = new OutputStream[rawDataSource.length];

        for (int i = 0; i < rawDataSource.length; i++) {
            Object[] streams = IOUtil.obtainStreams(rawDataSource[i]);
            inputStream[i] = (InputStream) streams[0];
            outputStream[i] = (OutputStream) streams[1];
            header[i] = new RawDataHeader512(inputStream[i]);
            
            startSamples[i] = ((long) startTimes[i]) * header[i].getSamplingFrequency();
            stopSamples[i] = (long) stopTimes[i];
            if (stopSamples[i] == Long.MAX_VALUE) {
                stopSamples[i] = header[i].getNumberOfSamples();	
            } else {
                stopSamples[i] = ((long) stopTimes[i]) * header[i].getSamplingFrequency();
            }


        }

        nElectrodes = header[0].getNumberOfElectrodes();
    }




    public void commenceWriting() throws IOException {
        for (int i = 0; i < outputStream.length; i++) {
            outputStream[i].write(1);
        }
    }


    public synchronized final void finishSampleProcessing() {
        // move to the next data set
        currentInputStream++;
    
        // maybe it's the last one
        if (currentInputStream == rawDataSource.length) {
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
            if (sampleListener4 != null) {
                try {
                    sampleListener4.finishSampleProcessing();
                } catch (IOException e) {
                    Vision.reportException("Could not finish processing for " +
                            sampleListener4.getClass().getName(), e);
                }
            }
            
            if (sampleListener5 != null) {
                try {
                    sampleListener5.finishSampleProcessing();
                } catch (IOException e) {
                    Vision.reportException("Could not finish processing for " +
                            sampleListener5.getClass().getName(), e);
                }
            }
            finished = true;
            return;
        }


        try {
            InputStream is = inputStream[currentInputStream];
            if (waitForData) is = new WaitingInputStream(is); // Wrapping input in WaitingInputStream so it will wait for data
            asynchronousInput = new AsynchronousInputStream(is, bufferSizeInBytes, nBuffers);
            sampleInputStream = new CompressedSampleInputStream(asynchronousInput, bufferSizeInBytes, nElectrodes,
                    currentInputStream == 0, startSamples[currentInputStream], stopSamples[currentInputStream]);
            sampleInputStream.addSampleListener(this);
            sampleInputStream.setInitializationListener(this);
            sampleInputStream.start();
            asynchronousInput.start();

            if (DEBUG_FLOW) {
                System.out.println(
                        "MCSIS starts reading: " + rawDataSource[currentInputStream]);
            }
        } catch (IOException e) {
            Vision.reportFatalException("Could not read raw data", e);
        }
    }

    public int getCumulativeNSamples() {
        int nSamples = 0;
        for (int i = 0; i < rawDataSource.length; i++) {
            long size;
            if (this.waitForData) size = -1; // We can't trust the file size if we are doing waitForData
            else size = IOUtil.getInputStreamSize(rawDataSource[i]);
            
            long nCurSamples;
            if (size == -1) {
                // network stream or waitForData
                nCurSamples = header[i].getNumberOfSamples();
            } else {
                nCurSamples = (size - header[i].getHeaderSize()) / header[i].getSampleSize();
            }
            
            // Can quit early if stopSamples is set for less than full file
            if (stopSamples[i] < nCurSamples) nCurSamples = stopSamples[i];

            nSamples += nCurSamples;
        }

        return nSamples;
    }

    public RawDataHeader512 getHeader() throws IOException {	
        return new RawDataHeader512(
            header[0].getTimeBase(),
            header[0].getTime(),
            header[0].getNumberOfElectrodes(),
            header[0].getSamplingFrequency(),
            getCumulativeNSamples(),
            header[0].getArrayID(),
            header[0].getFormat(),
            header[0].getDatasetIdentifier(), //FIXME - is that the correct ID ?
            header[0].getComment()            // FIXME: Would be nice to make some kind of concatenated comment?
        );
    }


    public synchronized boolean isFinished() {
        return finished;
    }


    public synchronized int getBufferUsage() {
        if (finished) {
            return 0;
        } else {
            return asynchronousInput.getBufferUsage();
        }
    }


    int nCallsToInitialize = 0;
    public void initialize(short[][] buffer) {
        nCallsToInitialize++;

        if (nCallsToInitialize != 1) {
            throw new Error("initialize was called multiple times");
        }

        if (initializationListener != null) {
            initializationListener.initialize(buffer);
        }
    }


    public final void processSample(short[] sample) {
        samplesRead++;

        if (sampleListener1 != null) {
            sampleListener1.processSample(sample);
        }
        if (sampleListener2 != null) {
            sampleListener2.processSample(sample);
        }
        if (sampleListener3 != null) {
            sampleListener3.processSample(sample);
        }
        if (sampleListener4 != null) {
            sampleListener4.processSample(sample);
        }
        if (sampleListener5 != null) {
            sampleListener5.processSample(sample);
        }
    }


    public void start() {
        finishSampleProcessing();
    }


    public void setInitializationListener(SampleInitializationListener listener) {
        if (initializationListener != null) {
            throw new IllegalArgumentException(
            "The Initialization Listener was already set");
        }

        this.initializationListener = listener;
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
        } else if (sampleListener4 == null) {
            sampleListener4 = listener;
        } else if (sampleListener5 == null) {
            sampleListener5 = listener;
        }else {
            throw new IllegalArgumentException("Too many SampleListeners (5 allowed)");
        }
    }


    public void removeSampleListener(SampleListener listener) {
        if (sampleListener1 == listener) {
            sampleListener1 = null;
        } else if (sampleListener2 == listener) {
            sampleListener2 = null;
        } else if (sampleListener3 == listener) {
            sampleListener3 = null;
        } else if (sampleListener4 == listener) {
            sampleListener4 = null;
        } else if (sampleListener5 == listener) {
            sampleListener5 = null;
        } else {
            throw new IllegalArgumentException("No such listener");
        }
    }


    public synchronized int getSamplesRead() {
        return samplesRead;
    }
    
    public int getNElectrodes() {
        return nElectrodes;
    }

}
