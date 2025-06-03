package edu.ucsc.neurobiology.vision.math.fitting;

import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.plot.*;


/**
 * A function that encapsulates another function with the goal of taking the absolute
 * value of it.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class AbsAdapterFunction
    implements FunctionData, PlotData {

    FunctionData f;


    public AbsAdapterFunction(FunctionData f) {
        this.f = f;
    }


    public double getValueAt(double t) throws CannotEvaluateException {
        return Math.abs(f.getValueAt(t));
    }


    public String getDescription() {
        return "Abs";
    }
}
