package edu.ucsc.neurobiology.vision.neuronviewer;

import static java.awt.Color.*;
import java.awt.*;
import javax.swing.tree.*;

import static edu.ucsc.neurobiology.vision.math.MathUtil.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.math.fitting.*;
import static edu.ucsc.neurobiology.vision.plot.PlotUtil.*;
import static edu.ucsc.neurobiology.vision.plot.SymbolType.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class TimeCoursePlotMaker
    extends PlotMaker {
    
    public boolean showZeroLine = true;
    public boolean normalize;
    public int backupSTAOffset;

    int staDepth, staOffset;
    double refreshInterval;
    private ScatterPlotStyle zeroLineStyle = new ScatterPlotStyle(
        "Zero Line", NONE, 0, black, true, black, 1);


    public TimeCoursePlotMaker() {
        super("Time Course", CLASS_AND_NEURON_PLOT);
    }


    public void initialize(NeuronViewer viewer) {
        super.initialize(viewer);

        if (staFile != null) {
            staDepth = staFile.getSTADepth();
            staOffset = staFile.getSTAOffset();
            refreshInterval = staFile.getRefreshTime();
        } else {
            int id = paramsFile.getIDList()[0];
            staDepth = paramsFile.getArrayCell(id, "GreenTimeCourse").length;
            staOffset = backupSTAOffset;
            refreshInterval = 8.33;
        }
    }


    public Component makePlot(IntegerList list, int plotType, TreePath classPath) {
        int[] neurons = list.toArray();

        PlotPanel timeCoursePanel = new PlotPanel("tc");



        // calculate the average TC
        double[][] averageTC = new double[3][staDepth];
        int n = 0;
        
        for (int i = 0; i < neurons.length; i++) {
            double[][] timeFilter = new double[3][];
            timeFilter[0] = paramsFile.getArrayCell(neurons[i], "RedTimeCourse");
            timeFilter[1] = paramsFile.getArrayCell(neurons[i], "GreenTimeCourse");
            timeFilter[2] = paramsFile.getArrayCell(neurons[i], "BlueTimeCourse");
            if (timeFilter[0] != null) {
                // calculate average tc
                for (int k = 0; k < 3; k++) {
                    MathUtil.add(averageTC[k], timeFilter[k]);
                }
                n++;
            }
        }


        // normalize and draw time courses
        double norm = 1.0;
        for (int i = 0; i < neurons.length; i++) {
            double[][] timeFilter = new double[3][];
            timeFilter[0] = paramsFile.getArrayCell(neurons[i], "RedTimeCourse");
            timeFilter[1] = paramsFile.getArrayCell(neurons[i], "GreenTimeCourse");
            timeFilter[2] = paramsFile.getArrayCell(neurons[i], "BlueTimeCourse");

            if (timeFilter[1] != null) {

                norm = Math.sqrt(sumSquare(timeFilter[0]) + sumSquare(timeFilter[1]) + 
                sumSquare(timeFilter[2]));
                if (normalize) {
                    
                    MathUtil.divide(timeFilter[0], norm);
                    MathUtil.divide(timeFilter[1], norm);
                    MathUtil.divide(timeFilter[2], norm);
                    
                }
                for (int c = 0; c < 3; c++) {
                    ScatterPlot sp = new ScatterPlot("" + neurons[i]);
                    for (int j = 0; j < timeFilter[1].length; j++) {
                        sp.add( (j - staDepth + 2 + staOffset) * refreshInterval, timeFilter[c][j]);
                    }
                    timeCoursePanel.addData(sp, timeCoursesStyle2[c]);
                }
            }
        }

        if (plotType == NEURON_PLOT) {
            int neuronID = list.get(0);
            try {
                STATimeFunction1 function = new STATimeFunction1(
                    paramsFile.getDoubleCell(neuronID, "a1"),
                    paramsFile.getDoubleCell(neuronID, "t1"),
                    paramsFile.getDoubleCell(neuronID, "a2"),
                    paramsFile.getDoubleCell(neuronID, "t2"),
                    paramsFile.getDoubleCell(neuronID, "a3"),
                    paramsFile.getDoubleCell(neuronID, "t3"),
                    paramsFile.getDoubleCell(neuronID, "n1"),
                    paramsFile.getDoubleCell(neuronID, "n2"),
                    paramsFile.getDoubleCell(neuronID, "n3"));
                
                if(normalize) { 	
                    function.setNormalization(norm);
                }
                    
                    
                timeCoursePanel.addData(function, new FunctionStyle("Fit"));
                
                
                
                
            } catch (Exception e) {
            }
        }
       // final int firstFrame =
        //    (int) ( (time - firstImageChangeTime) * conversionFactor) - staDepth +
        //    2 + staOffset;
        //Adds zero axis to scatter plot
       if (showZeroLine) {
           ScatterPlot zero = new ScatterPlot("");
           for (int i = 0; i < staDepth; i++) {
               zero.add( (i - staDepth + 2 + staOffset) * refreshInterval, 0);
           }
           timeCoursePanel.addData(zero, zeroLineStyle);
       }


        timeCoursePanel.setLabels("Time to spike (ms)", "STA (arb. units)");
        timeCoursePanel.autoscale();
        timeCoursePanel.padY();

        return timeCoursePanel;
    }

}
