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
public class LabviewSimulator {
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
        // if (args.length < 2 || args.length > 3) {
        //     System.out.println(errorMessage);
        //     System.exit(1);
        // }

        // fileName = args[0];

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
        
        OutputStream visionStream = null;
        Socket visionSocket = new Socket(hostAddress, 9876);
        visionStream = visionSocket.getOutputStream();
        InputStream vin = visionSocket.getInputStream();

        System.out.print("Sending spike finding command to Vision...");
        // visionStream.write(java.nio.ByteBuffer.allocate(4).putInt(34).array());

        // Test. Read a raw data file and stream it to Vision.
        FileInputStream inputStream = new FileInputStream("/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Data/20220406C/data011/data011000.bin");

        byte[] buffer = new byte[1024];
        int length;
        while ((length = inputStream.read(buffer)) != -1) {
            visionStream.write(buffer, 0, length);
        }

        inputStream.close();
        vin.close();
        visionStream.close();
        visionSocket.close();
    }

}
