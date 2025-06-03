package edu.ucsc.neurobiology.vision.testing;

import java.io.*;
import java.net.*;

import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * @author nobody, anyone can change
 */
public class TestStreaming {

    public static void test(String outputFolders) throws Exception {
        String[] fileNames = StringUtil.decomposeString(outputFolders, ";");
        int nOutputFiles = fileNames.length;
        Socket[] socket = new Socket[nOutputFiles];

        AsynchronousOutputStream[] bOutputStream = new AsynchronousOutputStream[
            nOutputFiles];
        OutputStream[] outputStreams = new OutputStream[fileNames.length];
        for (int i = 0; i < fileNames.length; i++) {
            if (fileNames[i].startsWith("net://")) { // its a network location
                String address = fileNames[i].substring(6);
                int indexOfSlash = address.indexOf('/');
                String host = address.substring(0, indexOfSlash);
                int port = Integer.parseInt(address.substring(indexOfSlash + 1));
                socket[i] = new Socket(host, port);
                outputStreams[i] = socket[i].getOutputStream();
            }
        }

        for (int i = 0; i < nOutputFiles; i++) {
            bOutputStream[i] = new AsynchronousOutputStream(outputStreams[i], 200, 100);
            bOutputStream[i].start();

            RawDataHeader newHeader = new RawDataHeader512(
                0,
                1 + 128,
                20000,
                0,
                RawDataHeader.composeArrayID(504, i + 1, nOutputFiles),
                RawDataHeader.FORMAT_12BIT_COMPRESSED,
                "test-test",
                "test");

//            outputStreams[i].write(newHeader.getBinaryRepresentation());
            bOutputStream[i].write(newHeader.getBinaryRepresentation());
        }

        byte[] buffer = new byte[1024 * 1024];
        for (int n = 0; n < 100; n++) {
            for (int i = 0; i < nOutputFiles; i++) {
//                outputStreams[i].write(buffer);
                bOutputStream[i].write(buffer);
            }
        }

        for (int i = 0; i < nOutputFiles; i++) {
            socket[i].shutdownOutput();
            //            outputStreams[i].close();
            bOutputStream[i].close();
        }
    }


    public static void main(String[] args) throws Exception {
        TestStreaming.test("net://198.202.70.247/1221");
    }
}
