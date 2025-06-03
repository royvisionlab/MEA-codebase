package edu.ucsc.neurobiology.vision.testing;

import java.util.*;

import edu.ucsc.neurobiology.vision.electrodemap.*;
import edu.ucsc.neurobiology.vision.plot.*;
import junit.framework.*;


/**
 * @author Dumitru
 */
public class TestElectrodeMap
    extends TestCase {

    protected void setUp() throws Exception {
        super.setUp();
    }


    protected void tearDown() throws Exception {
        super.tearDown();
    }


    public int displayElectrodeArrayCuts(ElectrodeMap map, int nParts) {
        PlotPanel p = new PlotPanel();
        p.addDataPlotter(new ElectrodeMapPlotter());
        p.addData(map, new ElectrodeMapStyle(""));
        p.setRange(map.getBounds());

        int total = 0, min = 10000, max = -10000;
        for (int part = 1; part <= nParts; part++) {
            Polygon2D region = map.getRegionShape(part, nParts);
            ScatterPlot sp = new ScatterPlot();
            for (int i = 0; i < region.nPoints; i++) {
                sp.add(region.xPoints[i], region.yPoints[i]);
            }
            p.addData(sp, "NONE 0 black, SOLID 1 red");

            ElectrodeMap m1 = map.getElectrodeMap(part, nParts);

            int n = (m1.getNumberOfElectrodes() - 1);
            if (n < min) {
                min = n;
            }
            if (n > max) {
                max = n;
            }
            total += n;
            System.err.println("part " + part + ": " + n);
        }

        System.err.println("Total: " + total);
        PlotUtil.showData("", p);
        return max - min;
    }


    private void testElectrodeArrayCuts(int arrayID, int ...nParts) {
        ElectrodeMap map = ElectrodeMapFactory.getElectrodeMap(arrayID);

        for (int n = 0; n < nParts.length; n++) {
//            System.out.println("Testing nParts = " + nParts[n]);
            int nElectodes = 0;
            ArrayList<Integer> electodes = new ArrayList();

            for (int pIndex = 1; pIndex <= nParts[n]; pIndex++) {
                ElectrodeMap m = map.getElectrodeMap(pIndex, nParts[n]);
                nElectodes += m.getNumberOfElectrodes();

                int[] parentElectodes = m.getParentElectrodeNumbers();
                for (int i = 0; i < parentElectodes.length; i++) {
                    if (!electodes.contains(parentElectodes[i])) {
                        electodes.add(parentElectodes[i]);
                    } else {
                        fail(" ===> Electrode " + parentElectodes[i] + " is duplicated");
                    }
                }
            }

            for (int i = 1; i < map.getNumberOfElectrodes(); i++) {
                if (!electodes.contains(i)) {
                    fail(" ===> Electrode " + i + " is not included in any subarray.");
                }
            }

            if (nElectodes - nParts[n] + 1 != map.getNumberOfElectrodes()) {
                fail("Number of electrodes - failed");
            }
        }
    }


    public void test512map() {
        int[] nParts = {1, 2, 4, 6};
        for (int i = 0; i < nParts.length; i++) {
//            int n = m.displayElectrodeArrayCuts(nParts[i]);
            testElectrodeArrayCuts(504, nParts[i]);
        }
    }


    public void test519map() {
        int[] nParts = {1, 6};
        for (int i = 0; i < nParts.length; i++) {
//            int n = m.displayElectrodeArrayCuts(nParts[i]);
            testElectrodeArrayCuts(1501, nParts[i]);
        }
    }


//    public static void main(String[] args) throws Exception {
//        Rectangular512ElectrodeMap electrodeMap = new Rectangular512ElectrodeMap();
//        TestElectrodeMap tester = new TestElectrodeMap();
//        tester.test512map();
//        tester.displayElectrodeArrayCuts(electrodeMap, 8);
//    }

    public static void main(String[] args) throws Exception {
    Rectangular512ElectrodeMap electrodeMap = new Rectangular512ElectrodeMap();
    TestElectrodeMap tester = new TestElectrodeMap();
    tester.test512map();
    tester.displayElectrodeArrayCuts(electrodeMap, 8);
}


}


