package edu.ucsc.neurobiology.vision.io;

import java.io.*;
import java.net.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.util.*;


/**
 * This class is used to simulate the DAQ computer. It read the Neurobiology Raw Data
 * from a file and sends the data over the network to the client.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class DAQ {
    private static String fileName;
    private static int port;
    private static String receiverIP;
    
    static InputStream fileStream;
    static OutputStream outputStream;
    static Socket clientSocket;

    static String commandReceiverIP = "128.114.130.109";
    static int commandReceiverPort = 9876;
    static int RUN_SPIKE_FINDING = 34;

    private static String errorMessage =
        "Please enter the following parameters: " +
        "\n - the file name, required." +
        "\n - the port on wich the DAQ works" +
        "\n - the IP for the receiver (default: " + commandReceiverIP + ")";

    public static void main(String args[]) throws Exception {
        if (args.length < 2 || args.length > 3) {
            System.out.println(errorMessage);
            System.exit(1);
        }

        fileName = args[0];
        port = Integer.parseInt(args[1]);
        if (args.length == 3) receiverIP = args[2];
        else receiverIP = commandReceiverIP;

        ServerSocket serverSocket = null;
        fileStream = (InputStream) (IOUtil.obtainStreams(fileName)[0]);
        RawDataHeader512 header = new RawDataHeader512(fileStream);
        File f = new File(fileName);
        header.setDatasetIdentifier(f.getParentFile().getName() + "-" +
                                    StringUtil.removeExtension(f.getName()));
//        header.setNumberOfSamples(
//            (int) (IOUtil.getInputStreamSize(fileName) - header.getHeaderSize()) /
//            header.getSampleSize());
        System.out.println(header);

        serverSocket = new ServerSocket(port, 5);
        System.out.println(
            "DAQ ready on IP " + InetAddress.getLocalHost().getHostAddress() + " port " +
            port);

        // let Vision know the DAQ is ready
        try {
            Socket s = new Socket(receiverIP, commandReceiverPort);
            new DataOutputStream(s.getOutputStream()).writeInt(RUN_SPIKE_FINDING);
            s.close();
        } catch (Exception ex) {
            System.err.println("\n" + ex.getClass().getName() + " : " + ex.getMessage());
            System.err.println("Could not send SPIKE FINDING command to Vision.");
            serverSocket.close();
            System.exit(1);
        }

        System.out.println("Listening for clients...");
        clientSocket = serverSocket.accept();

        System.out.println("Client connected from: " +
                           clientSocket.getInetAddress().getHostName());
        outputStream = clientSocket.getOutputStream();
        outputStream.write(header.getBinaryRepresentation());

        System.out.println("Sending data...");
        int nPerSecond = 100;
        int bufferSize = (20000 / nPerSecond) * 770;

        final byte[] buffer = new byte[bufferSize];
//        final long fileSize = new File(fileName).length();

        final Timer timer = new Timer(false);
        TimerTask t = new TimerTask() {
            public void run() {
                try {
                    int n = IOUtil.readFully(fileStream, buffer, 0, buffer.length);
                    outputStream.write(buffer, 0, n);
                    if (n != buffer.length) {
                        System.out.println("Closing the streams");
                        // we read all the data
                        timer.cancel();

                        fileStream.close();
                        outputStream.flush();
                        clientSocket.close();
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        };
        timer.scheduleAtFixedRate(t, 0, 10);
    }

}
