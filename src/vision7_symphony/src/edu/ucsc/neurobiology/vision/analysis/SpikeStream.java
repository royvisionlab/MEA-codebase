package edu.ucsc.neurobiology.vision.analysis;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class SpikeStream {
    int[] times, ttl;
    int tIndex, ttlIndex;


    public SpikeStream(int[] times, int[] ttl) {
        if (times == null) {
            throw new IllegalArgumentException("the times array cannot be null");
        }
        if (ttl == null) {
            throw new IllegalArgumentException("the ttl array cannot be null");
        }

        this.times = times;
        this.ttl = ttl;
    }


    public int getNext() {
        if (tIndex >= times.length && ttlIndex >= ttl.length) {
            return Integer.MAX_VALUE; // reached the end of the stream
        } else if (tIndex >= times.length) {
            return -ttl[ttlIndex++]; // return a TTL
        } else if (ttlIndex >= ttl.length) {
            return times[tIndex++]; // return a spike
        } else {
            if (times[tIndex] < ttl[ttlIndex]) {
                return times[tIndex++]; // return a spike
            } else {
                return -ttl[ttlIndex++]; // return a TTL
            }
        }
    }
}
