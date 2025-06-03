package edu.ucsc.neurobiology.vision.io;


/**
 * This is the base class for all spikes.  It contains the minimum
 * information about a spike, namely the electrode ID, the amplitude,
 * and the time of the spike.
 *
 * @author Charles A. Loomis, University of California, Santa Cruz<br>
 *         Dumitru Petrusca, University of California, Santa Cruz
 */
public class Spike
    implements Cloneable {

    public int time;
    public int electrode;
    public float amplitude;


    /**
     * Create a spike with the given electrode ID, time, amplitude and width.
     */
    public Spike(int time, int electrode, double amplitude) {
        this.time = time;
        this.electrode = electrode;
        this.amplitude = (float) amplitude;
    }


    public final void setValues(int time, int electrode, double amplitude) {
        this.time = time;
        this.electrode = electrode;
        this.amplitude = (float) amplitude;
    }


    /**
     * Performs a check according to the spike time. Used by SpikeBuffer.
     */
    public final boolean equals(Object obj) {
        return (obj == null) ? false : time == ( (Spike) obj).time;
    }


    /**
     * Need to implement this because equals was redefined.
     * @return int
     */
    public int hashCode() {
        return 42; // any arbitrary constant will do
    }


    public String toString() {
        return " T:" + time + " El:" + electrode;
    }


    public Object clone() {
        try {
            return super.clone();
        } catch (CloneNotSupportedException e) {
            return null;
        }
    }

}
