package edu.ucsc.neurobiology.vision.test;

import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.plot.*;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class AmplitudeWidth {

    public static void main(String args[]) throws Exception {
        final SpikeFile spikes = new SpikeFile("data051.spikes");
//        final WhiteNoiseMovie movie = WhiteNoiseMovie.load("data051.movie");

//        final SpikeFile spikes = new SpikeFile(new File("data051cf.spikes"));
//        final WhiteNoiseMovie movie = WhiteNoiseMovie.load("data051.movie");

        int el1 = 24, el2 = 33;

        //        DoubleHistogram2D h = new DoubleHistogram2D("", 0, 40, 0, 2000, 1, 20);

        ScatterPlot[] sp = new ScatterPlot[65];
        for (int i = el1; i < el2; i++) {
            sp[i] = new ScatterPlot("" + i);
        }

        SpikeIterator iter = spikes.iterator();
        int nSpikes = spikes.getSpikesCount(), sIndex = 0, oldPercent = -1;
        while (iter.hasNext()) {
            Spike s = iter.next();

            if (s.electrode >= el1 && s.electrode < el2) {
//                sp[s.electrode].add(s.width, s.amplitude);
            }

            sIndex++;
            int percent = (int) (100l * sIndex / nSpikes);
            if (percent != oldPercent) {
                System.out.println(percent);
                oldPercent = percent;
            }
        }

        for (int i = el1; i < el2; i++) {
            PlotUtil
                .showData("" + i, sp[i], new ScatterPlotStyle())
                .setLabels("Width (samples)", "Ampltitude (ADC counts)");
        }
    }

}
