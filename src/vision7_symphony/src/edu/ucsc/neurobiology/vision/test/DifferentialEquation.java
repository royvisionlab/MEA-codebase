package edu.ucsc.neurobiology.vision.test;

import edu.ucsc.neurobiology.vision.plot.*;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class DifferentialEquation {

    public DifferentialEquation() {
    }


    public static void main(String[] args) {
        ScatterPlotStyle style = new ScatterPlotStyle();
        style.setConnectingPoints(true);
        style.setSymbolType(SymbolType.NONE);

        double t0 = 0;
        double tMax = 50;
        double dt = 0.01;

        double A0 = 0.5, B = 1, S0 = 10, dS = 0.1, C = 0.2, D = 0.00047, E = 18.836;
//        double[] z = new double[ (int) Math.ceil( (tMax - t0) / dt)];
        int n = (int) Math.ceil( (tMax - t0) / dt);
        DoubleHistogram h = new DoubleHistogram("", t0, tMax, dt);
        for (int i = 0; i < h.getBinCount(); i++) {
            if (i > 500 && i < 510) {
                h.setBin(i, S0 + dS);
            } else {
                h.setBin(i, S0);
            }
        }

        double[] s = h.toArray();
        double a = A0 + D * E * S0 / (C + D * S0);
        double z = a * B / (a + S0);
        ;
        ScatterPlot spA = new ScatterPlot();
        ScatterPlot spZ = new ScatterPlot();
        for (int i = 0; i < n; i++) {
            double t = t0 + (i + 0.5) * dt;
//            a += dt * ( -C * (a - A0) + D * (E - (a - A0)) * s[i]);

            z += dt * (a * (B - z) - s[i] * z);

            spA.add(t, a);
            spZ.add(t, s[i] * z);
        }
        PlotUtil.showData("A", spA, style);
        PlotUtil.showData("Z", spZ, style);
    }
}
