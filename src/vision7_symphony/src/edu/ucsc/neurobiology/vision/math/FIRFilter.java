package edu.ucsc.neurobiology.vision.math;


/**
 * An implementation ofa Finite Impulse Response filter. Used to filter the raw data.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class FIRFilter {
    private double[] taps, delay;
    private int i;
    private double accum;
    private int nTaps;


    public FIRFilter(double[] taps) {
        this.taps = taps;
        this.nTaps = taps.length;
        this.delay = new double[nTaps];
    }


    /*
        public final double filterBasic(double input) {
            // store input at the beginning of the delay line
            delay[0] = input;
            // calculate FIR
            accum = 0;
            for (i = 0; i < coef.length; i++) {
                accum += coef[i] * delay[i];
            }
            //shift delay line
            for (i = coef.length - 2; i >= 0; i--) {
                delay[i + 1] = delay[i];
            }
            return accum;
        }
     */


    /*  circular buffer implementation */
    int state;
    public final double filter(double input) {
        // store input at the beginning of the delay line
        delay[state] = input;
        if (++state >= nTaps) { // incr state and check for wrap
            state = 0;
        }

        /* calc FIR and shift data */
        accum = 0;
        for (i = nTaps - 1; i >= 0; i--) {
            accum += taps[i] * delay[state];
            if (++state >= nTaps) { // incr state and check for wrap
                state = 0;
            }
        }

        return accum;
    }


    public int getTapsCount() {
        return nTaps;
    }
}
