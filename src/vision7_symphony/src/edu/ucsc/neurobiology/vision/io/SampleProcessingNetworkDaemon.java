package edu.ucsc.neurobiology.vision.io;

import java.io.*;
import java.net.*;

import edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMap;
import edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory;
import edu.ucsc.neurobiology.vision.util.StreamingSpikeFinder;

/**
 * Copied from NetworkDaemon, with just the simple addition of a sample processing stream to
 * do some analysis on the fly as data streams in.
 * 
 * As of 2010-11, this has not worked; a short analysis like noise finding can be put 
 * inline here successfully, but trying to run spike finding on the fly locks up acquisition.  
 * Could be revisited when streaming to ramdisks or with more powerful stream-to computers.
 *
 * Actually, also worth revisiting also with fewer analyses running during streaming.
 *
 * @author Peter H. Li, The Salk Institute
 */
public class SampleProcessingNetworkDaemon {
    static String usageStr = "Arguments Required: " +
                      "\n 1. the port number for the server to listen on" +
                      "\n 2. the buffer size in Kb (try 64)" +
                      "\n 3. the saving mode (0 - blocking, 1 - asynchronous)" +
                      "\n 4. the spike threshold" +
                      "\n 5. the ttl threshold" +
                      "\n 6. the mean time constant" +
                      "\n 7. the analysis output folder";

    static ServerSocket myService = null;
    static int bufferSize;
    
    public static void main(final String args[]) throws Exception {
        if (args.length != 7) {
            System.out.println(usageStr);
            return;
        }
        int port;
        try {
            port = Integer.parseInt(args[0]);
            bufferSize = Integer.parseInt(args[1]) * 1024;
        } catch (NumberFormatException e) {
            System.out.println(usageStr);
            return;
        }

        InetAddress host = InetAddress.getLocalHost();
        System.out.println("Local Host Name: " + host.getHostName());
        
        InetAddress[] inets = InetAddress.getAllByName(host.getHostName());
        
        //Find IP address from local net, but not a loopback.  If there is none, use the first on the list.
        boolean foundLocal = false;
        for (int i=0; i<inets.length; i++) {
            if (inets[i].isSiteLocalAddress() && !inets[i].isLoopbackAddress()) { 
                foundLocal = true;
                host = inets[i];
                break;
            }
        }
        
        if (!foundLocal) {
            host = inets[0];
        }

        System.out.println("Local Host Adress (clients should connect to): " +
                host.getHostAddress());
        

        System.out.print("Creating the server socket " + port + "...");
        myService = new ServerSocket(port, 0, host);
        System.out.println("done");
  
        while (true) {
            System.out.println("Waiting for a client to connect...");
            final Socket clientSocket = myService.accept();
            System.out.println("client connected from " +
                               clientSocket.getInetAddress().getHostAddress());

            Thread streamingThread = new Thread() {
                
                public void run() {
                    AsynchronousInputStream ais = null;
                    AsynchronousOutputStream aos = null;
                    OutputStream outputStream = null;
                    String ID = null;

                    PipedOutputStream pos = null;
                    PipedInputStream  pis = null;
                    
                    try {
                        InputStream inputStream = clientSocket.getInputStream();
                        RawDataHeader512 header = new RawDataHeader512(inputStream);
                        ID = header.getDatasetIdentifier().trim();

                        new File(header.getExperimentIdentifier()).mkdirs();
                        String fileName = header.getExperimentIdentifier() + File.separator +
                                header.getDatasetName() + ".bin";
                        outputStream = new FileOutputStream(fileName);
                        outputStream.write(header.getBinaryRepresentation());

                        pos = new PipedOutputStream();
                        pis = new PipedInputStream(pos);
                        CompressedSampleInputStream sampleInputStream = new CompressedSampleInputStream(pis, 500 * 770, header.getNumberOfElectrodes(), false, 0, -1);
                        ElectrodeMap electrodeMap = ElectrodeMapFactory.getElectrodeMap(header.getArrayID());
//                        StreamingSpikeFinder streamingSpikeFinder = new StreamingSpikeFinder(electrodeMap, Float.parseFloat(args[3]), Float.parseFloat(args[4]), Float.parseFloat(args[5]), args[6]);
//                        sampleInputStream.addSampleListener(streamingSpikeFinder);
                        sampleInputStream.start();

                        
                        byte[] buffer = new byte[bufferSize];
                        int lastRead;

                        if (args[2].equals("0")) {
                            System.out.println(
                                ID + ": Reading and writing data in blocking mode");
                            while ( (lastRead = inputStream.read(buffer)) != -1) {
                                outputStream.write(buffer, 0, lastRead);
                                pos.write(buffer, 0, lastRead);
                            }
                            inputStream.close();
                            outputStream.close();
                        } else if (args[2].equals("1")) {
                            System.out.println(
                                ID + ": Reading and writing data in asynchronous mode");
                            ais = new AsynchronousInputStream(inputStream, bufferSize,
                                100);
                            aos = new AsynchronousOutputStream(outputStream, bufferSize,
                                100);
                            ais.start();
                            aos.start();
                            while ( (lastRead = ais.read(buffer)) != -1) {
                                aos.write(buffer, 0, lastRead);
                                pos.write(buffer, 0, lastRead);
                            }
                            ais.close();
                            aos.close();
                            pos.close();
                            while (sampleInputStream.isAlive()) {
                                try {
                                    Thread.sleep(100);
                                } catch (InterruptedException e) {}
                            }
                        }

                        System.out.println(ID + ": All data received and saved.");
                    } catch (IOException ex) {
                        // exception occured, close gracefully
                        if (ais != null) {
                            try {
                                ais.close();
                            } catch (IOException e) {
                                System.out.println("Could not close the input stream.");
                            }
                        }

                        if (aos != null) {
                            try {
                                aos.close();
                            } catch (IOException e) {
                                System.out.println("Could not close the output stream.");
                            }
                        }
                        
                        if (pos != null) {
                            try {
                                pos.close();
                            } catch (IOException e) {
                                System.out.println("Could not close the pipe to SampleInputStream.");
                            }
                        }

                        System.out.println(
                            ID + " STREAMING THREAD CLOSED BECAUSE OF " + ex.getMessage());
                            ex.printStackTrace();
                    } // try
                } // Thread.run
            }; // anonymous Thread
            streamingThread.start();
        } // while(true)
    } // main
} // class SampleProcessingNetworkDaemon