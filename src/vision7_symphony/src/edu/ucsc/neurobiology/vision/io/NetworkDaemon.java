package edu.ucsc.neurobiology.vision.io;

import java.io.*;
import java.net.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;


/**
 * Waits indefinitely for raw data connections (acts like a server).
 * When connected asynchronously reads &
 * writes the raw data to an output file. Can handle more than one stream at a time.
 * The output files are usually analyzed with Vision to provide rapid feedback during
 * experiments.
 *
 * Automatically uses ip of current computer and reports it.  Will use the local net if
 * it is available.  Otherwise it will use the public net.
 *
 * Names of the form 10.*, 172.16.* - 172.31.*, 192.168.* are local, by RFC 1918.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class NetworkDaemon {
    static String USAGE = "Arguments Required: " +
                      "\n 1. the port number for the server to listen on" +
                      "\n 2. the buffer size in Kb (try 64)" +
                      "\n 3. the saving mode (0 - blocking, 1 - asynchronous) + ";


    static ServerSocket myService = null;
    static int bufferSize;


    public static void main(final String args[]) throws Exception {
        if (args.length != 3) {
            System.out.println(USAGE);
            return;
        }

        int port;
        try {
            port = Integer.parseInt(args[0]);
            bufferSize = Integer.parseInt(args[1]) * 1024;
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


        System.out.print("Creating the server socket " + port + "...");
        myService = new ServerSocket(port);
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

                    try {
                        System.out.print("Reading data header...");
                        InputStream inputStream = clientSocket.getInputStream();
                        RawDataHeader512 header = new RawDataHeader512(inputStream);
                        ID = header.getDatasetIdentifier().trim();
                        System.out.println("done");

                        System.out.print("Creating output paths...");
                        new File(header.getExperimentIdentifier()).mkdirs();
                        String fileName =
                            header.getExperimentIdentifier() + File.separator +
                            header.getDatasetName() + ".bin";
                        outputStream = new FileOutputStream(fileName);
                        outputStream.write(header.getBinaryRepresentation());
                        System.out.println("done");

                        byte[] buffer = new byte[bufferSize];
                        int lastRead;
                        if (args[2].equals("0")) {
                            System.out.println(
                                ID + ": Reading and writing data in blocking mode");
                            while ( (lastRead = inputStream.read(buffer)) != -1) {
                                outputStream.write(buffer, 0, lastRead);
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
                            }
                            ais.close();
                            aos.close();
                        }

                        System.out.println(ID + ": All data received and saved.");
                    } catch (IOException ex) {
                        // exception occurred, close gracefully
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

                        System.out.println(
                            ID + " STREAMING THREAD CLOSED BECAUSE OF " + ex.getMessage());
                    }
                }
            };
            streamingThread.start();
        }
    }

}
