package edu.ucsc.neurobiology.vision.testing;


/**
 * @author nobody, anyone can change
 */
public class ReadWriteTest
    extends Thread {

    static int kb = 1024;
    static long Mb = 1024 * 1024;

    public static void main(String[] args) throws Exception {
//        ReadTest t1 = new ReadTest(
//            "F:\\2005-01-21-0\\data999\\data999000.bin", 128 * Mb, 16 * kb);
        WriteTest t2 = new WriteTest("d:\\test1", 1280 * Mb, 16 * kb);
//        WriteTest t3 = new WriteTest("e:\\test2", 64 * Mb, 16 * kb);
//
        t2.start();
//        t2.start();
//        t3.start();
    }

}
