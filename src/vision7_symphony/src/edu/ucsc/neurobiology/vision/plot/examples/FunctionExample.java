package edu.ucsc.neurobiology.vision.plot.examples;

import java.awt.*;
import javax.swing.*;

import edu.ucsc.neurobiology.vision.math.fitting.*;
import edu.ucsc.neurobiology.vision.plot.*;


public class FunctionExample {

    public FunctionExample() {
    }


    public static void main(String args[]) throws Exception {
        // create the PlotPanel
        PlotPanel panel = new PlotPanel();
        // add the function with a standard FunctionStyle (right clich on the function
        // to customize the style). Look at the Gauss1DFunction and FunctionData classes
        // in the "plot" package

        Gaussian1DFunction f1 = new Gaussian1DFunction(0, 50, -5, 5);
        Gaussian1DFunction f2 = new Gaussian1DFunction(0, 50, +5, 5);
        Gaussian1DFunction f3 = new Gaussian1DFunction(0, -50, 0, 15);

        panel.addData(f1, new FunctionStyle("f1"));
        panel.addData(f2, new FunctionStyle("f", Color.red, 2));
        panel.addData(f3, new FunctionStyle("f3", Color.blue, 2));
//         FunctionSum fs = new FunctionSum(new FittableFunction[] {f1, f2, f3});
//         panel.addData(fs, new FunctionStyle("fs", Color.green, 2));

        // set the ranges
        panel.setRange( -20, 20, -100, 100);
        // set the axis lables
        panel.setLabels("x (m)", "y (m)");

        // display a Frame with the Function
        JFrame ff = new JFrame("Function Example");
        ff.setBounds(50, 50, 400, 400);
        ff.add(panel);
        ff.setVisible(true);
    }


}
