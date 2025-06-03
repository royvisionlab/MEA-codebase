package edu.ucsc.neurobiology.vision.testing;

import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;

import org.junit.Test;

import edu.ucsc.neurobiology.vision.io.IOUtil;
import edu.ucsc.neurobiology.vision.io.NeuronFile;
import edu.ucsc.neurobiology.vision.io.ParametersFile;
import edu.ucsc.neurobiology.vision.io.STAFile;


/**
 */
public class JUnitTests {

    public JUnitTests() {

    }


    @Test public void simpleAdd() throws Exception {
        STAFile master = new STAFile("F:\\vision5\\test data\\data002.sta");
        assertTrue(master.isSizeValid());

        String name = "test.sta";

        for (int n = 0; n < 10; n++) {
            new File(name).delete();
            STAFile f = new STAFile(
                name, master.getHeaderCapacity(), master.getWidth(),
                master.getHeight(), master.getSTADepth(), master.getSTAOffset(), master.getStixelWidth(), master.getStixelHeight(),
                master.getRefreshTime());

            for (int i = 1; i <= 2; i++) {
                f.addSTA(i, master.getSTA(3));
            }

            assertTrue(f.isSizeValid());
            assertTrue(f.getSTACount() == 2);
            f.close();

            f = new STAFile(name);
            assertTrue(f.isSizeValid());
            f.close();
        }
    }


    public static junit.framework.Test suite() {
        return new junit.framework.JUnit4TestAdapter(JUnitTests.class);
    }


    public static void main(String[] args) throws IOException {
//        org.junit.runner.JUnitCore.main(JUnitTests.class.getName());

        NeuronFile nf = new NeuronFile(
            "F:\\Good Data\\errors\\OOB\\data004-5\\data004-5.neurons");
        IOUtil.printArray(nf.getIDList());
        System.out.println(nf.getIDList().length);

        STAFile sf = new STAFile(
            "F:\\Good Data\\errors\\OOB\\data004-5\\data004-5.sta");
        IOUtil.printArray(sf.getIDList());
        System.out.println(sf.getIDList().length);

        ParametersFile pf = new ParametersFile(
            "F:\\Good Data\\errors\\OOB\\data004-5\\data004-5.params");
        IOUtil.printArray(pf.getIDList());
        System.out.println(pf.getIDList().length);
    }

}
