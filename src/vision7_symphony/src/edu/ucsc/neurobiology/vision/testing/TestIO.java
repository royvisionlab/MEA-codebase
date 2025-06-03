package edu.ucsc.neurobiology.vision.testing;

import java.io.*;

import edu.ucsc.neurobiology.vision.io.*;


/**
 * @author nobody, anyone can change
 */
public class TestIO {

    public static void testSpikeFile(String name) throws IOException {
        SpikeFile f = new SpikeFile(name);
        int[] ttl = f.getTTLTimes();
        System.err.println("TTLs " + ttl.length);
//        System.err.println(f.getSpikesCount());

//        for (int i = 0; i < 65; i++) {
//            System.err.println(i +" : " +f.getSpikesCount(i));
//        }

        SpikeIterator i = f.iterator();
        while (i.hasNext()) {
            Spike s = i.next();

            if (s == null) {
                System.err.println("Null Spike !!!");
            } else {
//                System.err.println(s);
            }
        }
    }


    public static void testVisionHeader(String[] args) throws IOException {
//        ClusteringModelFile m = new ClusteringModelFile(
//            "Y:\\data\\data-new\\2005-04-26-0\\data009\\data009.model");
//            System.out.println(m.getExtractionIds().length);

        VisionHeader m = new VisionHeader();
//        m.electrodeUsage = ElectrodeUsage.SEVEN_ELECTRODES;
        m.magic = 112;
        m.minThreshold = 4;

        RandomAccessFile f = new RandomAccessFile("a", "rw");
        IOUtil.writePublicFields(m, f);
        f.close();

        f = new RandomAccessFile("a", "r");
        IOUtil.readPublicFields(m, f);
        f.close();

        IOUtil.printPublicFields(m);
    }


    public static void testDirectoryInputStream() throws IOException {
        System.out.println(DirectoryInputStream.getDirectorySize(
            "f:\\2004-07-30-0\\data004") / 770.0);
    }


    public static void testSpikeStream1() throws IOException {
        NeuronFile f = new NeuronFile("f:\\data002.neurons");
        NeuronFile.SpikeStream1 ss = f.getSpikeStream();

        int n = 0;
        while (ss.hasNext()) {
            ss.next();
            n++;

            if (n % 1000000 == 0) {
                System.err.println(n);
            }
        }

        f.close();
    }


    public static void profileParamsFileOpen() {
        Thread t = new Thread() {
            public void run() {
                double t1 = System.currentTimeMillis();
                for (int i = 0; i < 10; i++) {
                    try {
                        final ParametersFile p = new ParametersFile(
                            "c:\\data003\\data003.params");
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
                double t2 = System.currentTimeMillis();
                System.out.println( (t2 - t1) / 1000.0);
            }
        };
        t.start();
    }


    public static void main(String[] args) throws Exception {
        testSpikeFile("F:\\data\\data003-12bit\\data003-12bit.spikes");
    }
}
