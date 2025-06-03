package edu.ucsc.neurobiology.vision.neuronviewer;

import static java.awt.Color.*;
import java.awt.*;
import javax.swing.tree.*;

import static edu.ucsc.neurobiology.vision.math.MathUtil.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.math.fitting.STATimeFunction1;
import static edu.ucsc.neurobiology.vision.plot.PlotUtil.*;
import static edu.ucsc.neurobiology.vision.plot.SymbolType.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class STVTimeCoursePlotMaker
extends PlotMaker {

    public boolean showVarianceLine = true;
    public boolean normalize;

    int staDepth;
    double refreshInterval;
    private ScatterPlotStyle varianceLineStyle = new ScatterPlotStyle(
            "Stimulus Variance Line", NONE, 0, black, true, black, 1);

    public double stimulusVariance; //Needed for the normalization to be correct


    public STVTimeCoursePlotMaker() {
        super("STV Time Course", CLASS_AND_NEURON_PLOT);
    }


    public void initialize(NeuronViewer viewer) {
        super.initialize(viewer);

        staDepth = viewer.staCollection.getSTADepth();
        refreshInterval = viewer.staCollection.getRefreshTime();
    }


    public Component makePlot(IntegerList list, int plotType, TreePath classPath) {
        int[] neurons = list.toArray();

        PlotPanel timeCoursePanel = null;

        timeCoursePanel = new PlotPanel("tcV");



        // calculate the average TC
        double[][] averageTC = new double[3][staDepth];
        int n = 0;
        for (int i = 0; i < neurons.length; i++) {
            double[][] timeFilter = new double[3][];

            timeFilter[0] = paramsFile.getArrayCell(neurons[i], "RedVTimeCourse");
            timeFilter[1] = paramsFile.getArrayCell(neurons[i], "GreenVTimeCourse");
            timeFilter[2] = paramsFile.getArrayCell(neurons[i], "BlueVTimeCourse");

            if (timeFilter[0] != null) {
                // calculate average tc
                for (int k = 0; k < 3; k++) {
                    MathUtil.add(averageTC[k], timeFilter[k]);
                }
                n++;
            }
        }

        // normalize and draw time courses
        for (int i = 0; i < neurons.length; i++) {
            double[][] timeFilter = new double[3][];

            timeFilter[0] = paramsFile.getArrayCell(neurons[i], "RedVTimeCourse");
            timeFilter[1] = paramsFile.getArrayCell(neurons[i], "GreenVTimeCourse");
            timeFilter[2] = paramsFile.getArrayCell(neurons[i], "BlueVTimeCourse");
            if (timeFilter[1] != null) {
                if (normalize) {

                    MathUtil.add(timeFilter[0], -stimulusVariance);
                    MathUtil.add(timeFilter[1], -stimulusVariance);
                    MathUtil.add(timeFilter[2], -stimulusVariance);

                    double norm = Math.sqrt(sumSquare(timeFilter[0]) + sumSquare(timeFilter[1]) + 
                            sumSquare(timeFilter[2]));
                    if (normalize) {

                        MathUtil.divide(timeFilter[0], norm);
                        MathUtil.divide(timeFilter[1], norm);
                        MathUtil.divide(timeFilter[2], norm);

                    }



                    MathUtil.add(timeFilter[0], stimulusVariance);
                    MathUtil.add(timeFilter[1], stimulusVariance);
                    MathUtil.add(timeFilter[2], stimulusVariance);
                }
                for (int c = 0; c < 3; c++) {
                    ScatterPlot sp = new ScatterPlot("" + neurons[i]);
                    for (int j = 0; j < timeFilter[1].length; j++) {
                        sp.add( (j - staDepth + 1) * refreshInterval, timeFilter[c][j]);
                    }
                    timeCoursePanel.addData(sp, timeCoursesStyle2[c]);
                }
            }
        }

        if (showVarianceLine) {
            ScatterPlot zero = new ScatterPlot("");
            for (int i = 0; i < staDepth; i++) {
                zero.add( (i - staDepth + 1) * refreshInterval, stimulusVariance);
            }
            timeCoursePanel.addData(zero, varianceLineStyle);
        }

        if (plotType == NEURON_PLOT) {
            int neuronID = list.get(0);
            try {
                STATimeFunction1 f = new STATimeFunction1(
                        paramsFile.getDoubleCell(neuronID, "a1V"),
                        paramsFile.getDoubleCell(neuronID, "t1V"),
                        paramsFile.getDoubleCell(neuronID, "a2V"),
                        paramsFile.getDoubleCell(neuronID, "t2V"),
                        paramsFile.getDoubleCell(neuronID, "a3V"),
                        paramsFile.getDoubleCell(neuronID, "t3V"),
                        paramsFile.getDoubleCell(neuronID, "n1V"),
                        paramsFile.getDoubleCell(neuronID, "n2V"),
                        paramsFile.getDoubleCell(neuronID, "n3V"),
                        paramsFile.getDoubleCell(neuronID, "tVOffset"));
                timeCoursePanel.addData(f, new FunctionStyle("Fit"));

//				double[] amps = edu.ucsc.neurobiology.vision.analysis.TimeCourseCalculator.findPeakAmps(f, -650, -40, .01, 10);
//				for(int i=amps.length-1; i>(amps.length-4); i--) {
//				System.out.println(amps[i]);
//				}
//				System.out.println();

            } catch (Exception e) {
            }
        }


        timeCoursePanel.setLabels("Time to spike (ms)", "STA (arb. units)");
        timeCoursePanel.autoscale();
        timeCoursePanel.padY();

        return timeCoursePanel;
    }
}
