package edu.ucsc.neurobiology.vision.plot.examples;

import javax.swing.*;

import edu.ucsc.neurobiology.vision.plot.*;


public class HistogramExample {

    public static void main(String args[]) throws Exception {
        // create the PlotPanel
        PlotPanel panel = new PlotPanel();
        // add the histogram with a standard HistogramStyle (right clich on the histogram
        // to customize the style)
        panel.addData(new RandomHistogram(), new HistogramStyle());
        // set the ranges
        panel.setRange( -20, +20, 0, 1);
        // set the axis lables
        panel.setLabels("the bin number", "Math.random()");

        // display a Frame with the Histogram
        JFrame f = new JFrame("Histogram Example");
        f.setBounds(50, 50, 900, 300);
        f.add(panel);
        f.setVisible(true);
    }

}


// Look at the documentation of HistogramData for a description of all the methods
class RandomHistogram
    implements HistogramData {
    private double[] bins;

    public RandomHistogram() {
        // 100 bin hitogram
        bins = new double[100];
        // filled with random numbers
        for (int i = 0; i < bins.length; i++) {
            bins[i] = Math.random();
        }
    }


    public String getDescription() {
        return "Example Random Histogram";
    }


    public int getBinCount() {
        return bins.length;
    }


    public double getBin(int bin) {
        return bins[bin];
    }


    public double getMin() {
        return -20;
    }


    public double getMax() {
        return +20;
    }


    public double getBinInterval() {
        return 1;
    }
}
