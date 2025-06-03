package edu.ucsc.neurobiology.vision.testing;

import edu.ucsc.neurobiology.vision.math.*;
import junit.framework.*;


/**
 * @author Dumitru
 */
public class TestMeanVarianceCalculator
    extends TestCase {


    public void testNoValues() {
        MeanVarianceCalculator mvc = new MeanVarianceCalculator(0);
        assertEquals(1, mvc.getMean());
        assertEquals(Double.NaN, mvc.getStandardDeviation());
    }


    public void testOnevalue() {
        MeanVarianceCalculator mvc = new MeanVarianceCalculator(0);
        mvc.add(34.4);
        assertEquals(34.4, mvc.getMean());
        assertEquals(Double.NaN, mvc.getStandardDeviation());
    }


    public void testTwovalue() {
        MeanVarianceCalculator mvc = new MeanVarianceCalculator(0);
        mvc.add(30);
        mvc.add(20);
        assertEquals(25.0, mvc.getMean());
        assertEquals(Math.sqrt(50), mvc.getStandardDeviation());
    }


    public void testRepeatedGet() {
        MeanVarianceCalculator mvc = new MeanVarianceCalculator(0);
        mvc.add(30);
        mvc.add(20);
        assertEquals(25.0, mvc.getMean());
        assertEquals(25.0, mvc.getMean());
        mvc.add(10);
        assertEquals(20.0, mvc.getMean());
    }
}
