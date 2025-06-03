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


public class DAQSimulation {
    static String USAGE = "Arguments Required: " +
                      "\n 1. the port number for the server to listen on" +
                      "\n 2. the ip address and port of the labview machine (should start with: net://)" +
                      "\n 3. the ip address and port of the Vision instance (should start with: net://)";


    static Socket clientSocket;

    // Set the port for the DAQSimulation.
    static int clientPort = 7887;
    static int serverPort = 9876;

    public static void main(final String args[]) throws Exception {
        // Parse inputs
        if (args.length > 4) {
            System.out.println(USAGE);
            return;
        }
        // try {
        //     if (args.length > 0) port = Integer.parseInt(args[0]);
        //     if (args.length > 1) rawDataSource       = args[1];
        //     if (args.length > 2) outputServerName    = args[2];
        // } catch (NumberFormatException e) {
        //     System.out.println(USAGE);
        //     return;
        // }

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

        // Create the client socket.
        clientSocket = new Socket(host, serverPort, host, clientPort);
    }

    // Parse the string inputs.
    private static String[] ParseInputAddress(String input) {
      System.out.println(input);
        input = input.replace("net://","");
        // Get rid of any white space.
        input = input.replaceAll("[\\s]", "");
        // String[] parts = input.split(File.separator);
        String[] parts = input.split("/");
        System.out.println(parts[0]);
        return parts;
    }

    private static String GetInputIPAddress(String input) {
      String[] parts = ParseInputAddress(input);
      System.out.println("Parsed client address is: " + parts[0]);
      return parts[0];
    }

    private static String GetInputPort(String input) {
      String[] parts = ParseInputAddress(input);
      return parts[1];
    }

}
