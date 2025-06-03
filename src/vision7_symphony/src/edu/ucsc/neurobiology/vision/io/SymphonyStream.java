package edu.ucsc.neurobiology.vision.io;

import java.io.*;
import java.net.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import edu.ucsc.neurobiology.vision.util.*;

public class SymphonyStream {

    static String USAGE = "Arguments Required: " +
        "\n 1. the port number for the server to listen on" +
        "\n 2. the ip address and port of the labview machine (should start with: net://)" +
        "\n 3. the ip address and port of the Vision instance (should start with: net://)";

    static ServerSocket myService = null;
    static String[] fileNames;
    static int port = 9876;
    static int bufferSize = 1024;

    static String outputServerName = "net://192.168.0.100/9000";
    static String rawDataSource = "net://192.168.0.101/7887";

    static String commonPath = "";
    static int bufferSizeInBytes = 1024 * 200;
    static int nBuffers = 200;
    static int secondsToStream = 0; // If set to 0, it will stream the whole epoch
    static int nSamplesToBuffer;
    static boolean waitForData = false;

    static Socket symphonySocket;
    static Socket labviewSocket;

    static MultipleCompressedSampleInputStream inputStream;

    static class PipeWriter extends Thread implements Runnable {
        InputStream in;
        OutputStream out;
        CountDownLatch latch;

        public PipeWriter(InputStream in, OutputStream out, CountDownLatch latch) {
            this.in = in;
            this.out = out;
            this.latch = latch;
        }

        public void run() {
            try {
                byte[] buffer = new byte[bufferSizeInBytes];
                int n = 0;
                while (in.available() != 0 && (n = in.read(buffer)) >= 0) {
//                  System.out.println("PipeWriter Processing: " + new String(Arrays.copyOfRange(buffer,0,n)));
                    out.write(buffer,0,n);
                }
                System.out.println("PipeWriter Terminating");
                in.close();
                out.close();
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }
    }

    public static void main(final String args[]) throws Exception {
        // Parse inputs
        if (args.length < 1 || args.length > 4) {
            System.out.println(USAGE);
            return;
        }
        try {
            if (args.length > 0) port = Integer.parseInt(args[0]);
        } catch (NumberFormatException e) {
            System.out.println(USAGE);
            return;
        }

        InetAddress host = InetAddress.getLocalHost();
        System.out.println("Local Host Name: " + host.getHostName());

        // Find IP address from local net, but not a loopback.  If there is none, use the first on the list.
        //
        // We used to use this to open the port on only one interface, but we switched to opening on all interfaces
        // so now this is just to report the private IP to the user.  The private IP is also not really needed though
        // because DNS should work if you just put the host name (e.g. creampuff-private) into spike finding.
        List<InetAddress> inets = new ArrayList<InetAddress>();
        List<NetworkInterface> netifaces = Collections.list(NetworkInterface.getNetworkInterfaces());
        for (NetworkInterface netiface : netifaces) inets.addAll(Collections.list(netiface.getInetAddresses()));
        host = inets.get(0);
        for (InetAddress inet : inets) {
            if (inet.isSiteLocalAddress() && !inet.isLoopbackAddress()) {
                host = inet;
                break;
            }
        }
        System.out.println("Private Net Address: " + host.getHostAddress());

        final String hostAddress = host.getHostAddress();

        System.out.print("Creating the server socket " + port + "...");
        myService = new ServerSocket(port);
        System.out.println("done");

        final PipedOutputStream pipedOutputStream = new PipedOutputStream();
        final PipedInputStream pipedInputStream = new PipedInputStream(pipedOutputStream);


        // Wait for connections.
        while (true) {
            System.out.println("Waiting for a client to connect...");
            final Socket clientSocket = myService.accept();

            System.out.println("Client connected from " +
                               clientSocket.getInetAddress().getHostAddress());

            String address = clientSocket.getInetAddress().getHostAddress();

            // Get the port of the incoming connection.
            InetSocketAddress sockaddr = (InetSocketAddress) clientSocket.getRemoteSocketAddress();
            int clientPort = sockaddr.getPort();

            // Cast the port to a string.
            String portString = Integer.toString(clientPort);

            if (portString.equals("7887")) {
                System.out.println("Received connection from LabView...");
            } else if (portString.equals("9000")) {
                System.out.println("Received connection from Vision...");
            }
            
            if (portString.equals("9000")) {
                System.out.println("Received connection from Vision...");

                symphonySocket = clientSocket;

                // Start the Symphony thread.
                Thread symphonyThread = new Thread(new Runnable() {
                    @Override
                    public void run() {
                        try {
                            int size = 0;
                            byte[] bytes = null;

                            while (true) {
                                // Read all of the bytes in the input stream.
                                if ((size = pipedInputStream.available()) != 0) {
                                bytes = new byte[size];
                                pipedInputStream.read(bytes, 0, bytes.length);
                                System.out.println(new String(bytes));

                                // Pass the file name to the output stream.
                                ObjectOutputStream oos = new ObjectOutputStream(symphonySocket.getOutputStream());
                                try {
                                    oos.writeObject(new String(bytes));
                                } finally {
                                    oos.flush();
                                }

                                // Pass the file name to Symphony for saving...
                                // ObjectOutputStream oos = new ObjectOutputStream(symphonySocket.getOutputStream());
                                // oos.writeObject(new String(bytes));
                                // oos.flush(); // Flush the stream.
                                // oos.close();
                                }
                            }
                        } catch (IOException e) {
                            e.printStackTrace();
                        } finally {
                            // Closing the streams
                            if (pipedOutputStream != null)
                                try {
                                    pipedOutputStream.close();
                                } catch (IOException e) {
                                    e.printStackTrace();
                                }
                            if (pipedInputStream != null) {
                                try {
                                    pipedInputStream.close();
                                } catch (IOException e) {
                                    e.printStackTrace();
                                }
                            }
                        }
                    }
                });
                symphonyThread.start();
            } else {
                System.out.println("Received connection from LabView...");
                
                Socket labviewSocket = clientSocket;

                Thread labviewThread = new Thread(new Runnable() {
                @Override
                public void run() {
                    try {
                        // InputStream inputStream = labviewSocket.getInputStream();
                        OutputStream outputStream = labviewSocket.getOutputStream();

                        System.out.println("Grabbing input stream from LabView...");
                        inputStream = new MultipleCompressedSampleInputStream(rawDataSource, bufferSizeInBytes, nBuffers, waitForData);

                        // Start the input stream.
                        inputStream.start();

                        if (rawDataSource.startsWith("net://")) {
                            System.out.println("Sending commence signal to DAQ computer");
                            inputStream.commenceWriting();
                        }

                        // Read the command integer at the front of the stream.
                        DataInputStream dis = new DataInputStream(labviewSocket.getInputStream());
                        int command = dis.readInt();

                        System.out.println("Reading data header...");
                        RawDataHeader512 header = inputStream.getHeader();

                        System.out.println("Parsing the file name...");
                        String fileName =
                            header.getExperimentIdentifier() + File.separator +
                            header.getDatasetName() + ".bin";

                        System.out.println("Output file name: " + fileName);

                        // Send the file name to the Symphony thread.
                        if (pipedOutputStream != null) {
                            try {
                                pipedOutputStream.write(fileName.getBytes());
                                System.out.println("Wrote file name to Symphony thread.");
                            } catch (IOException e) {
                                e.printStackTrace();
                            }
                        } else {
                            System.out.println("Piped output stream is null.");
                        }

                    // Set up the handshake with the output server.
                    // if (outputsExist) {
                    //     try {
                    //         System.out.println("Setting up handshake with Vision software...");
                    //         // Pipe the data to the Vision program.
                    //         dataSaver = new RawDataToVision(fileNames, commonPath, header,
                    //                 nSamplesToBuffer, nBuffers, secondsToStream);
                    //         inputStream.addSampleListener(dataSaver);
                    //         System.out.println("Handshake complete!");
                    //     } catch (IOException e) {
                    //         e.printStackTrace();
                    //     }
                    // }

                    } catch (IOException e) {}
                }
                });
                labviewThread.start();
            }



            Thread streamingThread = new Thread() {
                public void run() {
                    OutputStream visionStream = null;
                    OutputStream symphonyStream = null;

                    try {
                        System.out.print("Connecting to Labview...");
                        InputStream inputStream = clientSocket.getInputStream();
                        OutputStream outputStream = clientSocket.getOutputStream();

                        // Create the Vision and Symphony sockets.
                        Socket visionSocket = new Socket(hostAddress, 9000);
                        // Socket symphonySocket = new Socket(hostAddress, 9001);

                        // Get the output streams.
                        visionStream = visionSocket.getOutputStream();
                        // symphonyStream = symphonySocket.getOutputStream();

                        // Send the spike finding command to Vision.
                        System.out.print("Sending spike finding command to Vision...");
                        visionStream.write(java.nio.ByteBuffer.allocate(4).putInt(34).array());

                        // Wait for Vision to connect server side.
                        Socket visionSock = myService.accept();
                        InputStream vin = visionSock.getInputStream();
                        OutputStream vout = visionSock.getOutputStream();

                        // Tell Labview to commence writing data.
                        System.out.println("Sending commence signal to DAQ computer");
                        outputStream.write(1);
                        outputStream.close();

                        byte[] buffer = new byte[bufferSize];
                        int lastRead;
                        int totalRead = 0; // Keep track of the total data read.
                        while ( (lastRead = inputStream.read(buffer)) != -1) {
                            vout.write(buffer, 0, lastRead);

                            // totalRead += lastRead;
                            // if (totalRead < (512 * 1024)) {
                            //     symphonyStream.write(buffer, 0, lastRead);
                            // } else {
                            //     symphonyStream.close();
                            //     symphonySocket.close();
                            // }
                        }
                        inputStream.close();
                        visionStream.close();

                        vin.close();
                        vout.close();

                        // Close the sockets.
                        visionSocket.close();

                        // RawDataHeader512 header = new RawDataHeader512(inputStream);

                        
                        // String fileName =
                        //     header.getExperimentIdentifier() + File.separator +
                        //     header.getDatasetName() + ".bin";
                        // outputStream = new FileOutputStream(fileName);
                        // outputStream.write(header.getBinaryRepresentation());
                        // System.out.println("done");

                        System.out.println("All data received and saved.");
                    } catch (IOException ex) {
                        // exception occurred, close gracefully

                        System.out.println(
                            " STREAMING THREAD CLOSED BECAUSE OF " + ex.getMessage());
                    }
                }
            };
            streamingThread.start();
        }
    }
}
