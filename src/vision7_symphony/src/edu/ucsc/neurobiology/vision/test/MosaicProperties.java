package edu.ucsc.neurobiology.vision.test;

import static java.lang.Math.*;

import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.plot.*;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class MosaicProperties {

    public MosaicProperties() {
    }


    public static void main(String[] args) throws Exception {
//        String pname = "Y:\\data\\2005-09-09-1\\data002\\data002.params";
//        int[][] pairs = { {2341, 306}, {2341, 2111}, {2111, 1475}, {1475, 1740}, {1475,
//                        1740}, {1475, 1329}, {1475, 1282}, {1329, 1282}, {1329, 1740},
//                        {1329, 1500}, {1329, 703}, {1329, 1046}, {1329, 306}, {1500, 306},
//                        {1500, 1740}, {1500, 703}, {306, 1740}, {703, 1046}, {1282, 1046},
//        };

//        String pname = "Y:\\data\\2005-04-14-0\\data002\\data002.params";
//        int[][] pairs = { {242, 1065}, {242, 422}, {242, 238}, {242, 5}, {422, 238}, {422,
//                        703}, {422, 1065}, {703, 1065}, {5, 1484}, {5, 1646}, {1484, 1646},
//        };

        String pname = "Y:\\data\\2005-04-26-0\\data009\\data009.params";
        int[][] pairs = { {258, 675}, {258, 1945}, {258, 1670}, {1945, 1670}, {922, 675},
        };

        ParametersFile pf = new ParametersFile(pname);
        MeanVarianceCalculator mvc = new MeanVarianceCalculator();
        DoubleHistogram h = new DoubleHistogram("", 0, 10, 0.2);

        /*
                String cName = "All/OFF/YFast";

                int[] id = pf.getNeuronsInClass(cName);

                // get average sigma
                for (int i = 0; i < id.length; i++) {
                    double sigmaX = pf.getDoubleCell(id[i], "SigmaX");
                    double sigmaY = pf.getDoubleCell(id[i], "SigmaY");
                    mvc.add(sqrt(sigmaX * sigmaY));
                }
                double s = mvc.getMean();

                DoubleHistogram h = new DoubleHistogram("", 0, 10, 0.2);
                mvc.reset();

                        for (int i = 0; i < id.length; i++) {
                            double x0i = pf.getDoubleCell(id[i], "x0");
                            double y0i = pf.getDoubleCell(id[i], "y0");
                            double sigmaXi = pf.getDoubleCell(id[i], "SigmaX");
                            double sigmaYi = pf.getDoubleCell(id[i], "SigmaY");
                            double thetai = pf.getDoubleCell(id[i], "Theta");
                            double si = sqrt(sigmaXi * sigmaYi);

                            for (int j = i + 1; j < id.length; j++) {
                                double x0j = pf.getDoubleCell(id[j], "x0");
                                double y0j = pf.getDoubleCell(id[j], "y0");
                                double sigmaXj = pf.getDoubleCell(id[j], "SigmaX");
                                double sigmaYj = pf.getDoubleCell(id[j], "SigmaY");
                                double thetaj = pf.getDoubleCell(id[j], "Theta");
                                double sj = sqrt(sigmaXj * sigmaYj);

                 double d = sqrt( (x0i - x0j) * (x0i - x0j) + (y0i - y0j) * (y0i - y0j));
//                if (2*d / (si + sj) < 4) {
                                System.err.println("pair " + id[i] + " - " + id[j]);
                                h.fill(2 * d / (si + sj), 1);
                                mvc.add(2 * d / (si + sj));
//                }
                            }
                        }
         */

        for (int i = 0; i < pairs.length; i++) {
            double x0i = pf.getDoubleCell(pairs[i][0], "x0");
            double y0i = pf.getDoubleCell(pairs[i][0], "y0");
            double sigmaXi = pf.getDoubleCell(pairs[i][0], "SigmaX");
            double sigmaYi = pf.getDoubleCell(pairs[i][0], "SigmaY");
            double thetai = pf.getDoubleCell(pairs[i][0], "Theta");
            double si = sqrt(sigmaXi * sigmaYi);

            double x0j = pf.getDoubleCell(pairs[i][1], "x0");
            double y0j = pf.getDoubleCell(pairs[i][1], "y0");
            double sigmaXj = pf.getDoubleCell(pairs[i][1], "SigmaX");
            double sigmaYj = pf.getDoubleCell(pairs[i][1], "SigmaY");
            double thetaj = pf.getDoubleCell(pairs[i][1], "Theta");
            double sj = sqrt(sigmaXj * sigmaYj);

            double d = sqrt( (x0i - x0j) * (x0i - x0j) + (y0i - y0j) * (y0i - y0j));
            h.fill(2 * d / (si + sj), 1);
            mvc.add(2 * d / (si + sj));
        }

        System.err.println(mvc);
        PlotUtil.showData("", h, new HistogramStyle());
    }
}
