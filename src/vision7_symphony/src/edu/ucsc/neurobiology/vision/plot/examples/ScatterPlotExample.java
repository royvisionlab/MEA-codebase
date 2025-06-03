package edu.ucsc.neurobiology.vision.plot.examples;

import java.awt.*;
import javax.swing.*;

import edu.ucsc.neurobiology.vision.plot.*;


public class ScatterPlotExample {

    public static void main(String args[]) throws Exception {

        // Create a scatter plot
        ScatterPlot scatter = new ScatterPlot("A scatter plot");
        // fill it with some points
//        for (int i = 0; i < 100; i++) {
//            scatter.add(i, 2 * i + 50 + 25 * Math.random());
//        }
        scatter.add(100, 100);
        scatter.add(150, 150);
        scatter.add(200, 200);
        scatter.add(250, 250);

        // create the PlotPanel
        PlotPanel panel = new PlotPanel();
        // add the scatter plot with a standard ScatterPlotStyle (right clich on the
        // scatter plot to customize the style)
        panel.addData(scatter,
                      new ScatterPlotStyle(SymbolType.FILLED_SQUARE, 5, Color.black, true,
                                           Color.black, 1));
//        panel.addData(ElectrodeMapFactory.getElectrodeMap(0), new ElectrodeMapStyle("Electrode Map"));
//        panel.addData(ElectrodeMapFactory.create64Map(), null);
        // set the ranges
        panel.setRange(0, 300, 0, 300);
        // set the axis lables
        panel.setLabels("time (s)", "distance (m)");

        // display a Frame with the Histogram
        JFrame f = new JFrame("Scatter Plot Example");
        f.setDefaultCloseOperation(f.EXIT_ON_CLOSE);
        f.setBounds(50, 50, 500, 500);
        f.add(panel);
        f.setVisible(true);

    }

}
