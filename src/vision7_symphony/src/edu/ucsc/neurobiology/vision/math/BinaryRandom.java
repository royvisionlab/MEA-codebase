package edu.ucsc.neurobiology.vision.math;

import java.util.*;


/**
 * A boolean random number generator that returns "false" with the given probability.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class BinaryRandom {
    private Random random;
    private double p;


    public BinaryRandom(double p) {
        random = new Random(11111);
        this.p = p;
    }


    public boolean random() {
        double v = random.nextDouble();
        if (v > 1 - p) {
            return true;
        } else {
            return false;
        }
    }
}
