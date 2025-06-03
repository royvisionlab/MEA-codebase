package edu.ucsc.neurobiology.vision.io;

import java.io.*;
import java.net.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import edu.ucsc.neurobiology.vision.util.*;
// import java.nio.charset.Charset;
// import java.nio.charset.StandardCharsets;
// import java.util.stream.Collectors;


public class DataPipe {
    static String USAGE = "Arguments Required: " +
                      "\n 1. the port number for the server to listen on" +
                      "\n 2. the ip address and port of the labview machine (should start with: net://)" +
                      "\n 3. the ip address and port of the Vision instance (should start with: net://)";


    static ServerSocket myService = null;
    static MultipleCompressedSampleInputStream inputStream;

    static int port = 9876;
    static String outputServerName = "net://192.168.0.100/9000";
    static String rawDataSource = "net://192.168.0.101/7887";

    static String commonPath = "";
    static int bufferSizeInBytes = 1024 * 200;
    static int nBuffers = 200;
    static int secondsToStream = 0; // If set to 0, it will stream the whole epoch
    static int nSamplesToBuffer;
    static boolean waitForData = false;
    static String[] fileNames;
    static RawDataToVision dataSaver;

    static Socket symphonySocket;
    static Socket labviewSocket;

    public static void main(final String args[]) throws Exception {
        // Parse inputs
        if (args.length < 1 || args.length > 4) {
            System.out.println(USAGE);
            return;
        }
        try {
            if (args.length > 0) port = Integer.parseInt(args[0]);
            if (args.length > 1) rawDataSource       = args[1];
            if (args.length > 2) outputServerName    = args[2];
        } catch (NumberFormatException e) {
            System.out.println(USAGE);
            return;
        }

        // Use piped streams to communicate between the LabView and Symphony threads.
        // final PipedInputStream pipedInputStream = new PipedInputStream();
        // final PipedOutputStream pipedOutputStream = new PipedOutputStream(pipedInputStream);

        final PipedOutputStream pipedOutputStream = new PipedOutputStream();
        final PipedInputStream pipedInputStream = new PipedInputStream(pipedOutputStream);

        nSamplesToBuffer = bufferSizeInBytes / 770;
        bufferSizeInBytes = nSamplesToBuffer * 770;

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

        // Parse the output paths.
        final boolean outputsExist = !outputServerName.isEmpty();
        if (outputsExist) {
          fileNames = StringUtil.decomposeString(outputServerName, ";");
        }

        // int i = 570425344;
        // byte[] result = new byte[4];
        // result[0] = (byte) (i >> 24);
        // result[1] = (byte) (i >> 16);
        // result[2] = (byte) (i >> 8);
        // result[3] = (byte) (i /*>> 0*/);
        // System.out.println(result);

        // Wait for connections.
        while (true) {
            System.out.println("Waiting for a client to connect...");
            final Socket clientSocket = myService.accept();

            String address = clientSocket.getInetAddress().getHostAddress();

            // Get the port of the incoming connection.
            InetSocketAddress sockaddr = (InetSocketAddress) clientSocket.getRemoteSocketAddress();
            int clientPort = sockaddr.getPort();

            // Cast the port to a string.
            String portString = Integer.toString(clientPort);

            // Report the connection address and port.
            System.out.println("Client connected from " + address + " and port " + clientPort);

            final boolean foundInput = address.equals(GetInputIPAddress(rawDataSource));

            // Check the client address to determine which computer your getting.
            if (address.equals(GetInputIPAddress(rawDataSource)) || portString.equals("7887")) {
              System.out.println("Received connection from LabView...");

              Thread labviewThread = new Thread(new Runnable() {
                @Override
                public void run() {
                  try {

                    if (foundInput) {
                      System.out.println("Grabbing input stream...");
                      inputStream = new MultipleCompressedSampleInputStream(rawDataSource, bufferSizeInBytes, nBuffers, waitForData);

                      // Start the input stream.
                      inputStream.start();

                      if (rawDataSource.startsWith("net://")) {
                        System.out.println("Sending commence signal to DAQ computer");
                        inputStream.commenceWriting();
                      }

                      // Read the command integer at the front of the stream.
                      DataInputStream dis = new DataInputStream(clientSocket.getInputStream());
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
                      if (outputsExist) {
                        try {
                          System.out.println("Setting up handshake with Vision software...");
                          // Pipe the data to the Vision program.
                          dataSaver = new RawDataToVision(fileNames, commonPath, header,
                                  nSamplesToBuffer, nBuffers, secondsToStream);
                          inputStream.addSampleListener(dataSaver);
                          System.out.println("Handshake complete!");
                        } catch (IOException e) {
                          e.printStackTrace();
                        }
                      }
                    } else { // DAQSimulation connected
                      String fileName = "9999-99-99\\data999.bin";
                      try {
                        pipedOutputStream.write(fileName.getBytes());
                        System.out.println("Wrote file name to Symphony thread.");
                      } catch (IOException e) {
                        e.printStackTrace();
                      }
                    }

                  } catch (IOException e) {}
                }
              });
              labviewThread.start();
            } else if (!(address.equals(GetInputIPAddress(outputServerName)) && portString.equals(GetInputPort(outputServerName)))) {
              System.out.println("Received connection from Symphony...");
              symphonySocket = clientSocket;

              // final ObjectOutputStream oos = new ObjectOutputStream(symphonySocket.getOutputStream());

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
              System.out.println("Received input from " + GetInputIPAddress(outputServerName) + " and port " + GetInputPort(outputServerName));
            }
        }
    }

    // Parse the string inputs.
    private static String[] ParseInputAddress(String input) {
      System.out.println(input);
        input = input.replace("net://","");
        // Get rid of any white space.
        input = input.replaceAll("[\\s]", "");
        // String[] parts = input.split(File.separator);
        String[] parts = input.split("/");
        return parts;
    }

    private static String GetInputIPAddress(String input) {
      String[] parts = ParseInputAddress(input);
      return parts[0];
    }

    private static String GetInputPort(String input) {
      String[] parts = ParseInputAddress(input);
      return parts[1];
    }

}
