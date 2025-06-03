package edu.ucsc.neurobiology.vision.test;

import java.io.*;

import edu.ucsc.neurobiology.vision.analysis.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.plot.*;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ReversingGratingsDumitru {
    static String[] name = {
                           "f:\\data\\2005-04-26-0\\data005-mapped-data009",
                           "f:\\data\\2005-09-09-1\\data003",
                           "f:\\data\\2005-04-14-0\\data004"
    };
    static String[] pName = {
                            "f:\\data\\2005-04-26-0\\data009\\data009.params",
                            "f:\\data\\2005-09-09-1\\data002\\data002.params",
                            "f:\\data\\2005-04-14-0\\data002\\data002.params"
    };
    static String[] cName = {"All/OFF/Midget", "All/OFF/YFast", "All/OFF/Parasol"};


    static void test() throws Exception {
        DoubleHistogram nHist = new DoubleHistogram("", 0, 15, 0.5);
        int N = 0;

        for (int c = 1; c <= 2; c++) {
            for (int pIndex = 0; pIndex < name.length; pIndex++) {
                ReversingGratings m = new ReversingGratings(name[pIndex], 70);
                ParametersFile pFile = new ParametersFile(pName[pIndex]);

                // make n-index hist
                int[] id = pFile.getNeuronsInClass(cName[c]);
                for (int k = 0; k < id.length; k++) {
                    double[] f1 = pFile.getArrayCell(id[k], "T1reversingF1");
                    double[] f2 = pFile.getArrayCell(id[k], "T1reversingF2");
                    double maxN = Double.NEGATIVE_INFINITY;

                    for (int i = 0; i < f1.length - 1; i++) {
                        if (f1[i] > 0 || f2[i] > 0) {
                            double n = f2[i] / f1[i];
                            if (n > maxN) {
                                maxN = n;
                            }
                        }
                    }

                    N++;
                    nHist.fill(maxN, 1);
                }
            }
        }

        System.err.println(N);
        HistogramStyle h = new HistogramStyle();
        h.setFillingTowers(true);
        PlotPanel p = PlotUtil.showData("", nHist, h, 500, 300);
        p.setXRange( -0.05, 15);
    }


    static void simpleTest() throws IOException {
        ReversingGratings m = new ReversingGratings(
            "F:\\Good Data\\2005-04-26-0\\data005-mapped-data009", 30);

        int id = 31;
        m.setCurrentNeuron(id);

        m.showOnePeriodPSTH();
        m.showSpectra();
        m.showHarmonic(1);

//        int i = m.getRunID(768, 30, 0.25);
//        DoubleHistogram[] hh = m.getAverageResponseHistograms();
//        PlotPanel p = PlotUtil.showData(hh[i], new HistogramStyle(m.getRunInfo(i)));
//        p.addToLegend(m.getRunInfo(i));
    }


    public static void main(String[] args) throws IOException {
        simpleTest();
    }
}
