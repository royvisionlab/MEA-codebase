package edu.ucsc.neurobiology.vision.io;

import java.io.*;
import java.net.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import edu.ucsc.neurobiology.vision.util.VisionParams;


/**
 * Copied from NetworkDaemon, with the addition of breaking the incoming data into multiple files.
 * 
 * Tried to maintain full back-compatibility with NetworkDaemon so this can eventually drop-in replace.
 * 
 * Data must be broken into separate files at sample boundaries, not mid-sample, otherwise later 
 * analyses will get electrode numbers wrong.  Therefore, this daemon accepts an argument for file 
 * length in seconds of data only, and uses the header to determine how many bytes each file should 
 * be clipped to so as to break only at sample boundaries.
 * 
 * In this scheme it seemed useful to have the header saved to a separate file, rather than just to
 * the first data file.  But for backwards compatibility we also still save the header to the first
 * data file too.
 * 
 * @author Peter H. Li, The Salk Institute
 */
public class SplittingNetworkDaemon {
    static ServerSocket myService = null;
    static int port;
    static int bufferSize = 1024;
    static int saveMode = 1;
    static int fileLengthInSeconds = 0;
    
    static String USAGE = "Usage: " +
                      "\n 1. the port number for the server to listen on" +
                      "\n 2. the buffer size in samples, defaults to " + bufferSize +
                      "\n 3. the saving mode: 0 - blocking, 1 - asynchronous, defaults to " + saveMode +
                      "\n 4. the length of each file in seconds (0 = no splitting), defaults to " + fileLengthInSeconds;

    public static void main(final String args[]) throws Exception {
        // Parse inputs
        if (args.length < 1 || args.length > 4) {
            System.out.println(USAGE);
            return;
        }
        try {
            port = Integer.parseInt(args[0]);
            if (args.length > 1) bufferSize          = Integer.parseInt(args[1]) * 1024;
            if (args.length > 2) saveMode            = Integer.parseInt(args[2]);
            if (args.length > 3) fileLengthInSeconds = Integer.parseInt(args[3]);
        } catch (NumberFormatException e) {
            System.out.println(USAGE);
            return;
        }

        
        InetAddress host = InetAddress.getLocalHost();
        System.out.println("Local Host Name: " + host.getHostName());
        
        // Find IP address from local net, but not a loopback.  If there is none, use the first on the list.
        List<InetAddress> inets = new ArrayList<InetAddress>();
        List<NetworkInterface> netifaces = Collections.list(NetworkInterface.getNetworkInterfaces());
        for (NetworkInterface netiface : netifaces)
            inets.addAll(Collections.list(netiface.getInetAddresses()));
        host = inets.get(0);
        for (InetAddress inet : inets) {
            if (inet.isSiteLocalAddress() && !inet.isLoopbackAddress()) { 
                host = inet;
                break;
            }
        }

        System.out.println("Local Host Address (clients should connect to): " + host.getHostAddress());
        System.out.print("Creating the server socket " + port + "...");
        myService = new ServerSocket(port, 0, host);
        System.out.println("done");


        // Wait.  When connection comes, spawn a new thread to handle it and go back to waiting.
        while (true) {
            System.out.println("Waiting for a client to connect...");
            final Socket clientSocket = myService.accept();
            System.out.println("Client connected from " + clientSocket.getInetAddress().getHostAddress());

            Thread streamingThread = new Thread() {
                public void run() {
                    String ID = null;
                    InputStream inputStream = null;
                    OutputStream outputStream = null;
                    
                    try {
                        inputStream = clientSocket.getInputStream();
                        RawDataHeader512 header = new RawDataHeader512(inputStream);                        
                        writeHeaderFile(header);
                        ID = header.getDatasetIdentifier().trim();
                        
                        if (saveMode == 0) {
                            System.out.println(ID + ": Reading and writing data in blocking mode");
                        } else if (saveMode == 1) {
                            System.out.println(ID + ": Reading and writing data in asynchronous mode");
                            AsynchronousInputStream  ais = new AsynchronousInputStream(inputStream, bufferSize, 100);
                            ais.start();
                            inputStream = ais;
                        }

                        int fileNum = -1; // default
                        if (fileLengthInSeconds > 0) fileNum = 0;
                        outputStream = buildOutputStream(header, fileNum);
                        outputStream.write(header.getBinaryRepresentation()); // for backwards compatibility

                        int sampleSize = header.getSampleSize();
                        int frequency = header.getSamplingFrequency();
                        long fileLengthInBytes = ((long) header.getSampleSize()) * header.getSamplingFrequency() * fileLengthInSeconds;
                        long bytesWritten = 0;
                        byte[] buffer = new byte[bufferSize];
                        int lastRead;
                        while ((lastRead = inputStream.read(buffer)) != -1) {
                            if (fileLengthInBytes > 0 && lastRead + bytesWritten > fileLengthInBytes) {
                                int toWriteOld = (int) (fileLengthInBytes - bytesWritten);
                                outputStream.write(buffer, 0, toWriteOld);
                                outputStream.close();
                                
                                // Make new outputStream, put the rest of the bytes in
                                int toWriteNew = lastRead - toWriteOld;
                                outputStream = buildOutputStream(header, ++fileNum);
                                outputStream.write(buffer, toWriteOld, toWriteNew);
                                bytesWritten = toWriteNew;
                            } else {
                                outputStream.write(buffer, 0, lastRead);
                                bytesWritten += lastRead;
                            }
                        }
                        
                        inputStream.close();
                        outputStream.close();
                        System.out.println(ID + ": All data received and saved.");
                    } catch (IOException ex) {
                        try {
                            inputStream.close();
                            outputStream.close();
                        } catch (IOException e) {
                            System.out.println("Could not close the input/output streams.");
                        }
                        System.out.println(ID + " STREAMING THREAD CLOSED UNEXPECTEDLY: " + ex.getMessage());
                    }
                }
            };
            
            streamingThread.start();
        }
    }


    private static String getBaseFilename(RawDataHeader header) {
        return header.getExperimentIdentifier() + File.separator + header.getDatasetName();
    }

    
    private static void writeHeaderFile(RawDataHeader header) throws IOException {
        new File(header.getExperimentIdentifier()).mkdirs();
        String fileName = getBaseFilename(header) + VisionParams.HEADER_FILE_EXTENSION;
        OutputStream os = new FileOutputStream(fileName);
        os.write(header.getBinaryRepresentation());
    }

    private static OutputStream buildOutputStream(RawDataHeader header, int fileNum) throws IOException {
        new File(header.getExperimentIdentifier()).mkdirs();
        
        String fileName = getBaseFilename(header);
        if (fileNum >= 0) {
            new File(header.getExperimentIdentifier() + File.separator + header.getDatasetName()).mkdirs();
            fileName += File.separator + header.getDatasetName();
            fileName += String.format("%03d", fileNum);
        }
        fileName += VisionParams.BIN_FILE_EXTENSION_512;

        System.out.println("Creating output file " + fileName);
        OutputStream outputStream = new FileOutputStream(fileName);
        
        if (saveMode == 1) {
            AsynchronousOutputStream aos = new AsynchronousOutputStream(outputStream, bufferSize, 100);
            aos.start();
            outputStream = aos;
        }
        
        return outputStream;
    }
    
}