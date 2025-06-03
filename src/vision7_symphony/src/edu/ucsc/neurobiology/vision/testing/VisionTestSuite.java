package edu.ucsc.neurobiology.vision.testing;

import junit.framework.*;


/**
 * @author Dumitru
 */
public class VisionTestSuite
    extends TestCase {

    public VisionTestSuite(String s) {
        super(s);
    }


    public static Test suite() {
        TestSuite suite = new TestSuite();
        suite.addTestSuite(edu.ucsc.neurobiology.vision.testing.TestElectrodeMap.class);
        suite.addTestSuite(edu.ucsc.neurobiology.vision.testing.
                           TestMeanVarianceCalculator.class);
        return suite;
    }
}
