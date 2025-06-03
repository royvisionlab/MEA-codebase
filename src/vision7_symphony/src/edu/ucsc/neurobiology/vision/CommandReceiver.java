package edu.ucsc.neurobiology.vision;

import java.io.*;
import java.net.*;
import java.util.*;


/**
 * This class waits on port 9876 for the command 34. After receiving this
 * command it start the online spike finding process.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class CommandReceiver
    implements Runnable {

    ServerSocket myService;
    int commandRceiverPort = 9000;
    int RUN_SPIKE_FINDING = 34;


    public CommandReceiver() {}


    public void start() throws IOException {
        System.out.print("Creating Command Server... ");
        myService = new ServerSocket(commandRceiverPort);

        InetAddress host = InetAddress.getLocalHost();
        System.out.println(host.getHostName() + ", " + host.getHostAddress());

        new Thread(this).start();
    }


    public void run() {
        System.out.println("Waiting for commands...");

        while (true) {
            Socket clientSocket;
            DataInputStream dis;
            final int command;

            try {
                clientSocket = myService.accept();
            } catch (IOException ex) {
                System.err.println("Could not establish connection");
                continue;
            }

            try {
                dis = new DataInputStream(clientSocket.getInputStream());
                command = dis.readInt();
            } catch (IOException ex) {
                System.err.println("Could not read command");
                continue;
            }

            System.out.println("Command " + command + " received from " +
                               clientSocket.getInetAddress().getHostAddress());

            Thread commandThread = new Thread() {
                public void run() {
                    if (command == RUN_SPIKE_FINDING) {
                        spikeFinding();
                    } else {
                        System.err.println("Unknown Command " + command);
                    }
                }
            };
            commandThread.start();
        }
    }


    public static void spikeFinding() {
        String cName = "Spike Finding";
        HashMap<String, String> p = Vision.getConfig().getParameterList(cName);
        Vision.getInstance().getCalculationManager().runCalculation(cName, p);
    }

}
