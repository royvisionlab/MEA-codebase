package edu.ucsc.neurobiology.vision.analysis;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.calculations.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.stimulus.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 * @owner Matthew Grivich, The Salk Institute
 */
public class STACalculation
    extends AbstractCalculation {

    private String filePath;
    private String rawMoviePath;
    private BufferedMovie movie1, movie2;
    private NeuronFile neuronFile;
    private int staDepth;
    private int staOffset;
    private int atOnce;
    private int spikesToCalculate;
    boolean threaded;
    boolean calculateSTV;
    boolean resumeCalculation;

    public void startCalculation() throws Exception {
        Vision app = Vision.getInstance();
        long startTime = System.currentTimeMillis();

        String fileName = new File(filePath).getName();
        String neuronFileName = filePath + File.separator + fileName +
                                VisionParams.NEURON_FILE_EXTENSION;
        String globalsFileName = filePath + File.separator + fileName + ".globals";
        if (rawMoviePath == null || rawMoviePath.length() == 0) {
            String movieFileName = filePath + File.separator + fileName +
                                   VisionParams.MOVIE_FILE_EXTENSION;
            WhiteNoiseMovie m1 = new WhiteNoiseMovie(movieFileName, globalsFileName);
            movie1 = new BufferedMovie(m1);
            movie2 = new BufferedMovie(m1.duplicate());
        } else {
            String movieFileName = rawMoviePath;

            RawMovie rawMovie = new RawMovie(movieFileName, globalsFileName);
            movie1 = new BufferedMovie(rawMovie);
            movie2 = new BufferedMovie(rawMovie);
        }

        app.sendMessage("Loading NeuronFile");
        neuronFile = new NeuronFile(neuronFileName);

        app.sendMessage("Creating STAFile");
        STAFile staFile = null;
        STAFile steFile = null;
        String staFileName = StringUtil.removeExtension(neuronFileName) +
                             VisionParams.STA_FILE_EXTENSION;
        String steFileName = StringUtil.removeExtension(neuronFileName) + ".stv";
        if (new File(staFileName).exists() && resumeCalculation) {
            staFile = new STAFile(staFileName);
            if (new File(steFileName).exists()) {
                steFile = new STAFile(steFileName);
            }
        } else {
            staFile = new STAFile(
                StringUtil.removeExtension(neuronFileName) +
                VisionParams.STA_FILE_EXTENSION, 10000,
                movie1.getWidth(), movie1.getHeight(), staDepth,staOffset,
                movie1.getFrame(0).getStixelWidth(), movie1.getFrame(0).getStixelHeight(), movie1.getRefreshTime());

            if (calculateSTV) {

                steFile = new STAFile(
                    StringUtil.removeExtension(neuronFileName) + ".stv", 10000,
                    movie1.getWidth(), movie1.getHeight(), staDepth, staOffset,
                    movie1.getFrame(0).getStixelWidth(), movie1.getFrame(0).getStixelHeight()
                    , movie1.getRefreshTime());
            }
        }
        int[] ttl = neuronFile.getTTLTimes();
//        System.out.println("First TTL : " + ttl[0]);
        final int[] idList = neuronFile.getIDList();
        double tps = 1000.0 / neuronFile.getSamplingFrequency();

        // choose a reasonable "atOnce". make sure is not 1
        if (idList.length % atOnce == 1) {
            atOnce++;
        }

        int startIndex =  staFile.getIDList().length;  //0 if new file, larger if partial file

        while (true) {
            final int n = Math.min(idList.length - startIndex, atOnce);
            System.out.println("Calculating STAs for " + n + " neurons");

            app.startProgressBar();
            STACalculator[] sta = calculateSTAsThreaded(
                startIndex, n, neuronFile, movie1, movie2, staDepth, staOffset, ttl[0], tps,
                threaded, spikesToCalculate);
            app.endProgressBar();

            // save the STAs
            app.sendMessage("Saving " + n + " STAs");
            System.out.println("Saving " + n + " STAs");
            for (int i = 0; i < n; i++) {
                STA s = sta[i].getSTA();
                staFile.addSTA(idList[startIndex + i], s);
                if (calculateSTV) {
                    steFile.addSTA(idList[startIndex + i], sta[i].getSTV());
                }
            }

            if (n != atOnce) {
                break;
            }
            startIndex += n;

//            // FIXME
//            if (true) {
//                break;
//            }
        }

        neuronFile.close();
        staFile.close();
        if (steFile != null) {
            steFile.close();
        }

        app.sendMessage(
            "Done in: " + (System.currentTimeMillis() - startTime) / 1000. + " s.");

        Vision.getInstance().getCalculationManager().calculationDone();
    }


    static class STAThread
        extends Thread {

        public STACalculator[] sta;
        public int progress = 0;

        private final int startIndex, n;
        private BufferedMovie movie;
        private NeuronFile neuronFile;
        private final int staDepth, staOffset;
        private final int firstTTL;
        private final double timePerSample;
        private final int spikesToCalculate;


        public STAThread(
            int startIndex, int n, NeuronFile neuronFile, BufferedMovie movie,
            int staDepth, int staOffset, int firstTTL, double timePerSample, int spikesToCalculate) {

            this.staDepth = staDepth;
            this.staOffset = staOffset;
            this.startIndex = startIndex;
            this.n = n;
            this.movie = movie;
            this.neuronFile = neuronFile;
            this.firstTTL = firstTTL;
            this.timePerSample = timePerSample;
            this.spikesToCalculate = spikesToCalculate;

            sta = new STACalculator[n];
        }


        public void run() {
            try {
                movie.reset();
                Vision.getInstance().sendMessage("Loading spike times");

                int[] idList = neuronFile.getIDList();
                int[][] t = new int[n][]; //time of spike [neuron][spike]
                int totalNumberOfSpikes = 0;
                for (int i = 0; i < n; i++) {
                    try {
                        t[i] = neuronFile.getSpikeTimes(idList[startIndex + i]);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                    totalNumberOfSpikes += t[i].length;
                }

                Vision.getInstance().sendMessage("Creating STACalculators");
                for (int i = 0; i < n; i++) {
                    sta[i] = new STACalculator(staDepth, staOffset, firstTTL, timePerSample, movie);
                }

                Vision.getInstance().sendMessage("Processing Spikes");
                int min = 0;
                int spikeIndex = 0;
                int[] spikesProcessed = new int[n]; //number of spikes processed [neuron]
                while (true) {
                    min = 0; //neuron with the minimum number of processed spikes
                    //find first neuron with spikes remaining.
                    //if no neurons remain, break
                    boolean done = false;
                    while (!done) {
                        if (min == n) {
                            done = true;
                        } else if (spikesProcessed[min] < t[min].length) {
                            done = true;
                        } else {
                            min++;
                        }
                    }
                    if (min == n) {
                        break;
                    }

                    //for each neuron, find the one with the earliest spike.
                    //call this min.
                    for (int i = 0; i < n; i++) {
                        //if i neuron has spikes to process
                        if (spikesProcessed[i] < t[i].length) {
                            //if i neuron current spike time <
                            //min neuron current spike time
                            //set i neuron to min neuron
                            if (t[i][spikesProcessed[i]] < t[min][spikesProcessed[min]]) {
                                min = i;
                            }
                        }
                    }

                    //if spikes processed for min neuron < spikesToCalculate,
                    //process it.
                    if (spikesProcessed[min] < spikesToCalculate) {
                        sta[min].addSpike(t[min][spikesProcessed[min]]);
                    }

                    //increment spikes processed for min neuron
                    spikesProcessed[min]++;

                    //increment total spikes processed
                    spikeIndex++;
                    progress = (int) (100L * spikeIndex / totalNumberOfSpikes);
                }

                // very important
                for (int i = 0; i < n; i++) {
                    sta[i].finish();
                }
            } catch (IOException e) {
                e.printStackTrace();
                return;
            }
        }
    };

    public static STACalculator[] calculateSTAsThreaded(
        final int startIndex, final int n, final NeuronFile neuronFile,
        final BufferedMovie movie1, final BufferedMovie movie2, final int staDepth, int staOffset,
        final int firstTTL, final double timePerSample, boolean threaded,
        int spikesToCalculate) {

        Vision app = Vision.getInstance();
        if (threaded) {
            final int n1 = n / 2;
            STAThread t1 = new STAThread(startIndex, n1, neuronFile, movie1, staDepth, staOffset,
                                         firstTTL, timePerSample, spikesToCalculate);
            STAThread t2 = new STAThread(startIndex + n1, n - n1, neuronFile, movie2,
                                         staDepth, staOffset, firstTTL, timePerSample,
                                         spikesToCalculate);
            t1.start();
            t2.start();

            while (t1.isAlive() || t2.isAlive()) {
                int progress = (t1.progress + t2.progress) / 2;
                app.setProgress(progress);
                try {
                    Thread.sleep(100);
                } catch (InterruptedException e) {}
            }

            STACalculator[] sta = new STACalculator[n];
            System.arraycopy(t1.sta, 0, sta, 0, t1.sta.length);
            System.arraycopy(t2.sta, 0, sta, n1, t2.sta.length);
            return sta;
        } else {
            STAThread t = new STAThread(
                startIndex, n, neuronFile, movie1, staDepth, staOffset, firstTTL, timePerSample,
                spikesToCalculate);
            t.start();

            while (t.isAlive()) {
                app.setProgress(t.progress);
                try {
                    Thread.sleep(100);
                } catch (InterruptedException e) {}
            }

            return t.sta;
        }
    }

    public void setParameters(HashMap<String, String> parameters) {
        filePath = new File(parameters.get("File_Path")).getAbsolutePath();
        staDepth = Integer.parseInt(parameters.get("STA Depth"));
        staOffset = Integer.parseInt(parameters.get("STA Offset"));
        spikesToCalculate = Integer.parseInt(parameters.get("Spikes To Calculate"));
        if (spikesToCalculate <= 0)
            spikesToCalculate = Integer.MAX_VALUE;
        atOnce = Integer.parseInt(parameters.get("STAs Calculated At Once"));
        rawMoviePath = parameters.get("Raw_Movie_Path");
        threaded = Boolean.valueOf(parameters.get("Double Threaded")).
                   booleanValue();
        calculateSTV = Boolean.valueOf(parameters.get("Calculate STV")).booleanValue();
        resumeCalculation = Boolean.valueOf(parameters.get("Resume Calculation")).booleanValue();
    }
}