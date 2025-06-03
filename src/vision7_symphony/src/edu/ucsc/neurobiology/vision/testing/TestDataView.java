package edu.ucsc.neurobiology.vision.testing;

import java.util.*;

import javax.swing.*;

import edu.ucsc.neurobiology.vision.analysis.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.plot.*;


/**
 * @author nobody, anyone can change
 */
public class TestDataView {
    static double[] xValues;


    public static void main(String[] args) throws Exception {
        final HashMap vMap = new HashMap();

        int deltaT = 100 * 20; // in samples
        int dt = 1 * 20;
        NeuronFile neuronFile = new NeuronFile(
            "D:\\new\\2005-04-26-0\\data026\\data026.neurons");
        ParametersFile paramsFile = new ParametersFile(
            "D:\\new\\2005-04-26-0\\data022\\data022.params");
        System.out.println("Files opened");

        String classPath1 = "All/ON/Parasol";
        String classPath2 = "All/OFF/Parasol";
        HashMap<Integer, Integer> idMap1 =
            paramsFile.evaluate("ID", "classID==\"" + classPath1 + "\"");
        HashMap<Integer, Integer> idMap2 =
            paramsFile.evaluate("ID", "classID==\"" + classPath2 + "\"");

        int[] idList1 = new int[idMap1.size()];
        int[] idList2 = new int[idMap2.size()];
        int i = 0;
        for (Integer id : idMap1.keySet()) {
            idList1[i++] = id.intValue();
        }
        i = 0;
        for (Integer id : idMap2.keySet()) {
            idList2[i++] = id.intValue();
        }

        for (int i1 = 0, k = 0; i1 < idList1.length; i1++) {
            System.out.println(i1);
            for (int i2 = i1 + 1; i2 < idList2.length; i2++) {
                int[] t1 = neuronFile.getSpikeTimes(idList1[i1]);
                int[] t2 = neuronFile.getSpikeTimes(idList2[i2]);

                DoubleHistogram ccH = CrossCorrelationCalculator.
                                      getCrossCorrelationHistogram(
                                          t1, t2, dt, deltaT, null);
                if (k == 0) {
                    xValues = ccH.getXValues();
                }
                ccH.scale(1 / ccH.getBinSum());

                vMap.put(idList1[i1] + "-" + idList2[i2], ccH.toArray());
                k++;
            }
        }

        PCACalculator pca = new PCACalculator(vMap);
        final HashMap<Integer, Double> xMap = pca.getPCAComponent(0);
        final HashMap<Integer, Double> yMap = pca.getPCAComponent(1);
        ScatterPlot sp = new ScatterPlot();
        for (Object key : xMap.keySet()) {
            double x = xMap.get(key).doubleValue();
            double y = yMap.get(key).doubleValue();
            sp.add(x, y);
        }
        PlotPanel p = new PlotPanel();

        p.addSelectionAction(new SelectionAction("Print IDs") {
            public void selectionPerformed(JComponent source, Selection selection) {
                SelectionRegion r = selection.getSelection();
                for (Object key : xMap.keySet()) {
                    double x = xMap.get(key).doubleValue();
                    double y = yMap.get(key).doubleValue();
                    if (r.contains(x, y)) {
                        System.out.println(key);
                    }
                }
            }
        });
        p.addSelectionAction(new SelectionAction("Show") {
            public void selectionPerformed(JComponent source, Selection selection) {
                SelectionRegion r = selection.getSelection();
                PlotPanel p = new PlotPanel();
                ScatterPlotStyle style = new ScatterPlotStyle();
                style.setConnectingPoints(true);
                style.setSymbolType(SymbolType.NONE);

                for (Object key : xMap.keySet()) {
                    double x = xMap.get(key).doubleValue();
                    double y = yMap.get(key).doubleValue();
                    if (r.contains(x, y)) {
                        p.addData(
                            new ScatterPlot(xValues, (double[]) vMap.get(key), null),
                            style);
                    }
                }
                p.autoscale();
                PlotUtil.showData("", p);
            }
        });
//            p.setLabels(expression + " : PCA 1", expression + " : PCA 2");

        ScatterPlotStyle style = new ScatterPlotStyle();
        style.setSymbolSize(3);
        p.addData(sp, style);
        p.autoscale();
        p.pad();
        PlotUtil.showData("PCA ", p);

    }

}
