package edu.ucsc.neurobiology.vision.io;

import java.io.*;
import java.net.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.electrodemap.*;


/**
 * Splits the raw data into a certain number of parts depending on how many output
 * locations are provided and sends the data asynchronously to those locations. The
 * output locations can be local files, or TCP/IP network connections. Data can be split
 * in parts only as allowed by the ElectrodeMap class associated with it.
 *
 * @see ElectrodeMap
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class RawDataToVision
    implements SampleListener {

    private AsynchronousOutputStream bOutputStream[];
    private byte[][] buffer;
    private int[] bufferIndex;
    private ElectrodeMap map;
    private int[][] electrodeList;
    private int nOutputFiles;
    // for how long to keep streaming data
    private final int samplesToStream;
    private boolean stopped = false;
    private boolean exceptionOccured = false;
    private boolean[] saveExpandedTTL;


    public RawDataToVision(String[] fileNames, String commonPath, RawDataHeader header,
                        final int nSamplesToBuffer, final int nBuffers,
                        int secondsToStream) throws IOException {

        System.out.println("\n---> Saving Raw Data in " + fileNames.length +
                           " file(s):");
        OutputStream[] outputStreams = new OutputStream[fileNames.length];
        for (int i = 0; i < fileNames.length; i++) {
            if (fileNames[i].startsWith("net://")) { // it's a network location
                String host = null;
                int port = 0;
                try {
                    String address = fileNames[i].substring(6);
                    int indexOfSlash = address.indexOf('/');
                    host = address.substring(0, indexOfSlash);
                    port = Integer.parseInt(address.substring(indexOfSlash + 1));
                } catch (Exception ex) {
                    Vision.reportException(
                        "Wrong network path: " + fileNames[i] + "." +
                        "The format is: net://x.x.x.x/port", ex);
                }
                Socket socket = new Socket(host, port);
                outputStreams[i] = socket.getOutputStream();
            } else { // it's a file
                fileNames[i] = commonPath + File.separator + fileNames[i] +
                               File.separator + header.getExperimentIdentifier();
                new File(fileNames[i]).mkdirs();
                fileNames[i] += File.separator + header.getDatasetName() + ".bin";
                outputStreams[i] = new FileOutputStream(fileNames[i]);
            }
            System.out.println("Location " + (i + 1) + ": " + fileNames[i]);
        }

        this.samplesToStream = 20000 * secondsToStream;
        this.nOutputFiles = outputStreams.length;
        map = ElectrodeMapFactory.getElectrodeMap(header.getArrayID());

        int nSamples = (secondsToStream > 0) ?
                       Math.min(samplesToStream, header.getNumberOfSamples())
                       : header.getNumberOfSamples();

        bOutputStream = new AsynchronousOutputStream[nOutputFiles];
        electrodeList = new int[nOutputFiles][];
        buffer = new byte[nOutputFiles][];
        bufferIndex = new int[nOutputFiles];
        saveExpandedTTL = new boolean[nOutputFiles];

        for (int i = 0; i < nOutputFiles; i++) {
            int[] eList = map.getElectrodeMap(i + 1, nOutputFiles).
                          getParentElectrodeNumbers();
            saveExpandedTTL[i] = (eList.length + 1) % 2 == 1;
            int sampleSize;

            if (saveExpandedTTL[i]) {
                electrodeList[i] = eList;
                sampleSize = 2 + (electrodeList[i].length) * 3 / 2;
            } else { // add the electrode 0 (TTL) in the beginning so it will be compressed
                electrodeList[i] = new int[eList.length + 1];
                System.arraycopy(eList, 0, electrodeList[i], 1, eList.length);
                electrodeList[i][0] = 0;
                sampleSize = (electrodeList[i].length) * 3 / 2;
            }

            buffer[i] = new byte[sampleSize * nSamplesToBuffer];
            /*
                        System.err.println(
                            "Part " + (i + 1) + ": " + eList.length +
                            ", expanded = " + saveExpandedTTL[i] +
                            ", sample size = " + sampleSize
                            );
             */
            // create the output file
            bOutputStream[i] = new AsynchronousOutputStream(
                outputStreams[i], sampleSize * nSamplesToBuffer, nBuffers);
            bOutputStream[i].start();

            System.out.println("Writing header...");
            RawDataHeader newHeader = new RawDataHeader512(
                0,
                eList.length + 1,
                header.getSamplingFrequency(),
                nSamples,
                RawDataHeader.composeArrayID(header.getArrayID(), i + 1, nOutputFiles),
                header.getFormat(),
                header.getDatasetIdentifier(),
                header.getComment());
            // Try to write the spike finding command to the stream.
            int command = 34;
            bOutputStream[i].write(command);
            // byte[] cbyte = new byte[1];
            // cbyte[0] = (byte) command;
            // bOutputStream[i].write(cbyte);
            // Now, write the header to the output stream.
            bOutputStream[i].write(newHeader.getBinaryRepresentation());
        }

        //        System.out.println("n " + nSamples / 20000.);
    }


    public double getBufferUsage(int bufferID) {
        return bOutputStream[bufferID].getBufferUsage();
    }


    public int getOutputStreamsCount() {
        return bOutputStream.length;
    }


    int samples = 0;
    public void processSample(short[] sample) {
        if (stopped || exceptionOccured) {
            return;
        }
        samples++;

        // write the TTL for all files
//        for (int i = 0; i < nOutputFiles; i++) {
//            buffer[i][bufferIndex[i]++] = (byte) (sample[0] >> 8);
//            buffer[i][bufferIndex[i]++] = (byte) (sample[0] & 0x00ff);
//        }

        for (int i = 0; i < nOutputFiles; i++) {
            // write the TTL uncompressed if the total number of electrodes is ODD
            if (saveExpandedTTL[i]) {
                buffer[i][bufferIndex[i]++] = (byte) (sample[0] >> 8);
                buffer[i][bufferIndex[i]++] = (byte) (sample[0] & 0x00ff);
            }

            for (int e = 0; e < electrodeList[i].length; ) {
                int s1 = sample[electrodeList[i][e++]] + 2048;
                int s2 = sample[electrodeList[i][e++]] + 2048;

                buffer[i][bufferIndex[i]++] = (byte) (s1 >> 4);
                buffer[i][bufferIndex[i]++] = (byte) ( ( (s1 & 0x000f) << 4) + (s2 >> 8));
                buffer[i][bufferIndex[i]++] = (byte) (s2 & 0x00ff);
            }
        }

        if (bufferIndex[0] == buffer[0].length) {
            for (int i = 0; i < nOutputFiles; i++) {
                try {
                    bOutputStream[i].write(buffer[i], 0, bufferIndex[i]);
                    bufferIndex[i] = 0;
                } catch (Exception e) {
                    exceptionOccured = true;
                    System.out.println(
                        "DATA SAVING/STREAMING WAS STOPPED BECAUSE OF " +
                        e.getClass().getName() + " : " + e.getMessage());
                    System.out.println(
                        "THE SPIKE DISPLAY AND THE DAQ WILL CONTINUE TO WORK PROPERLY. ");
                }
            }
        }

        // see whether we have to stop streaming
        if (samplesToStream > 0 && samples >= samplesToStream) {
            Vision.getInstance().sendMessage("As demanded the streaming was stopped at " +
                                             samplesToStream / 20000 + " sec.");

            try {
                finishSampleProcessing();
            } catch (IOException e) {
                Vision.reportException(e);
            }
        }
    }


    public void finishSampleProcessing() throws IOException {
        if (stopped || exceptionOccured) {
            return;
        }

        // saving last samples and close the output streams files
        for (int i = 0; i < nOutputFiles; i++) {
            bOutputStream[i].write(buffer[i], 0, bufferIndex[i]);
            bufferIndex[i] = 0;
            bOutputStream[i].close();
        }

        stopped = true;
    }

}
