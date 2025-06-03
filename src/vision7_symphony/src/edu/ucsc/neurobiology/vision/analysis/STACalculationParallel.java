package edu.ucsc.neurobiology.vision.analysis;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.calculations.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.io.chunk.ChunkFile;
import edu.ucsc.neurobiology.vision.io.chunk.GlobalsFile;
import edu.ucsc.neurobiology.vision.io.chunk.WhiteNoiseMovieFile;
import edu.ucsc.neurobiology.vision.stimulus.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * Copied from STACalculation.java, but with multithreading implemented differently to scale better to many threads.
 *
 * STACalculation multithreads by giving different threads subsets of the total set of cell IDs.  In testing, this 
 * gives only modest speed-up.  I believe this is because every thread is still looping over almost the entire 
 * WhiteNoiseMovie, which is inefficient.
 *
 * STACalculationParallel multithreads by giving each thread a subregion of the full WhiteNoiseMovie to calculate.
 * In tests, this seems to scale better.
 *
 * @author Peter H. Li, The Salk Institute
 * @owner  Peter H. Li, The Salk Institute
 */
public class STACalculationParallel extends AbstractCalculation {
    private String filePath;
    private String rawMoviePath;
    private NeuronFile neuronFile;
    private int staDepth, staOffset;
    private int atOnce;		// How many to do at once on this machine, to avoid memory overflow
    private int spikesToCalculate;
    private int numThreads;			// How many subregions to break the movie into and process separately
    private int numParts, partNum;	// For farming the STA calculation out to multiple machines
    private float[] proportions; // For farming the STA calculation out to multiple machines
    private boolean calculateSTV;
    private boolean resumeCalculation;
    
    public void startCalculation() throws Exception {
        Vision app = Vision.getInstance();
        
        long startTime = System.currentTimeMillis();
        
        String fileName = new File(filePath).getName();
        String neuronFileName  = filePath + File.separator + fileName + VisionParams.NEURON_FILE_EXTENSION;
        String globalsFileName = filePath + File.separator + fileName + VisionParams.GLOBALS_FILE_EXTENSION;
        String movieFileName   = filePath + File.separator + fileName + VisionParams.MOVIE_FILE_EXTENSION;
        
        // Override default movie path?
        if (rawMoviePath != null && rawMoviePath.length() > 0) {
            movieFileName = rawMoviePath;
        }
        
        app.sendMessage("Loading Movie Params"); // Only loaded at this point to get parameters and frames for STAFile
        WhiteNoiseMovieFile movieFile = new WhiteNoiseMovieFile(movieFileName, ChunkFile.READ);
        WhiteNoiseMovieFile.WhiteMovieParams movieParams = movieFile.getWhiteMovieParams();
        GlobalsFile globalsFile = new GlobalsFile(globalsFileName, ChunkFile.READ);
        GlobalsFile.RunTimeMovieParams runTime = globalsFile.getRunTimeMovieParams();
        
        
        app.sendMessage("Creating STAFile");
        
        String partSuffix = "";
        if (numParts > 1) partSuffix = "-" + partNum;
        
        STAFile staFile = null;
        STAFile steFile = null;
        String staFileName = StringUtil.removeExtension(neuronFileName) + VisionParams.STA_FILE_EXTENSION + partSuffix;
        String steFileName = StringUtil.removeExtension(neuronFileName) + VisionParams.STE_FILE_EXTENSION + partSuffix;
        if (new File(staFileName).exists() && resumeCalculation) {
            staFile = new STAFile(staFileName);
            if (new File(steFileName).exists()) {
                steFile = new STAFile(steFileName);
            }
        } else {
            staFile = new STAFile(
                staFileName, 10000,
                movieParams.width, movieParams.height, staDepth, staOffset,
                runTime.micronsPerStixelX, runTime.micronsPerStixelY, runTime.refreshPeriod);

            if (calculateSTV) {
                steFile = new STAFile(
                    steFileName, 10000,
                    movieParams.width, movieParams.height, staDepth,staOffset,
                    runTime.micronsPerStixelX, runTime.micronsPerStixelY,
                    runTime.refreshPeriod);
            }
        }
        
        
        app.sendMessage("Loading NeuronFile");
        neuronFile = new NeuronFile(neuronFileName);        
        int[] ttl = neuronFile.getTTLTimes();
        double tps = 1000.0 / neuronFile.getSamplingFrequency();

        final int[] fullIDList = neuronFile.getIDList();
        final int[] idList;
        if (proportions == null) idList = ParallelUtil.getMyPart(partNum, numParts, fullIDList);
        else idList = ParallelUtil.getMyPart(partNum, numParts, fullIDList, proportions);
        
        // choose a reasonable "atOnce". make sure is not 1
        if (idList.length % atOnce == 1) atOnce++;
        
        int startIndex = staFile.getIDList().length;  // 0 if new file, larger if partial file
        
        app.sendMessage("Seeding white noise movie frames");
        BufferedMovie[] movies = buildMoviesThreaded(movieFileName, globalsFileName, numThreads, staDepth);
        
        while (startIndex < idList.length) {
            final int n = Math.min(atOnce, idList.length - startIndex);
            final int endIndex = startIndex + n - 1;
            
            app.sendMessage("Loading spike times");
            int[][] spikeTimes = new int[n][]; //time of spike [neuron][spike]
            int totalNumberOfSpikes = 0;
            for (int i = 0; i < n; i++) {
                try {
                    spikeTimes[i] = neuronFile.getSpikeTimes(idList[startIndex + i]);
                } catch (IOException e) {
                    e.printStackTrace();
                }
                totalNumberOfSpikes += spikeTimes[i].length;
            }
            
            app.sendMessage("Calculating frame indices");
            STACalculator dummyConverter = new STACalculator(staDepth, staOffset, ttl[0], tps, movies[0]);
            int[][] allFirstFrameIndices = new int[n][];
            for (int i = 0; i < n; i++)
                allFirstFrameIndices[i] = dummyConverter.spikeTimesToFirstFrameIndices(spikeTimes[i]);
            
            // Show user which neurons we're working on.
            // Only show N/M if there will be multiple rounds.
            String nOutOfN = "" + n;
            if (endIndex < idList.length-1) nOutOfN += "/" + (idList.length-startIndex);
            String theseIDs = "" + idList[startIndex] + "-" + idList[endIndex];
            String calcMsg = "Calculating STAs for " + nOutOfN + " remaining neurons (IDs " + theseIDs + ")";

            app.sendMessage(calcMsg);
            System.out.println(calcMsg);
            
            app.startProgressBar();
            STAThread[] staThreads = calculateSTAsThreaded(n, neuronFile, allFirstFrameIndices, 
                    totalNumberOfSpikes, staDepth, staOffset, ttl[0], tps, numThreads, spikesToCalculate, movies);
            app.endProgressBar();

            // Save the STAs
            String saveMsg = "Saving " + n + " STAs";
            app.sendMessage(saveMsg);
            System.out.println(saveMsg);
            for (int i = 0; i < n; i++) {
                STA[] staSubregions = new STA[numThreads];
                for (int j = 0; j < numThreads; j++) staSubregions[j] = staThreads[j].stacs[i].getSTA();
                STA sta = STA.concatenatingBuild(staSubregions, movieParams.width, movieParams.height);
                staFile.addSTA(idList[startIndex + i], sta);
                
                if (calculateSTV) {
                    STA[] stvSubregions = new STA[numThreads];
                    for (int j = 0; j < numThreads; j++) stvSubregions[j] = staThreads[j].stacs[i].getSTV();
                    STA stv = STA.concatenatingBuild(stvSubregions, movieParams.width, movieParams.height);
                    steFile.addSTA(idList[startIndex + i], stv);
                }
            }

            startIndex += n;
        }

        neuronFile.close();
        staFile.close();
        if (steFile != null) steFile.close();

        app.sendMessage("Done in: " + (System.currentTimeMillis() - startTime) / 1000. + " s.");
        Vision.getInstance().getCalculationManager().calculationDone();
    }

    
    static class SubmovieThread extends Thread {
        private final String movieFileName;
        private final String globalsFileName;
        private final int thisThread;
        private final int numThreads;
        private WhiteNoiseMovieSubregion wnms;
        
        public SubmovieThread(String movieFileName, String globalsFileName, int thisThread, int numThreads) {
            this.movieFileName   = movieFileName;
            this.globalsFileName = globalsFileName;
            this.thisThread      = thisThread;
            this.numThreads   	 = numThreads;
        }
        
        public void run() {
            try {
                wnms = WhiteNoiseMovieSubregion.threadWhiteNoiseMovieSubregion(movieFileName, globalsFileName, thisThread, numThreads);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        
        public WhiteNoiseMovieSubregion getWNMS() {
            return wnms;
        }
        
        public static boolean anySubmovieThreadsAlive(SubmovieThread[] submovieThreads) {
            for (SubmovieThread s : submovieThreads) {
                if (s.isAlive()) return true;
            }
            return false;
        }
    }
    
    
    static class STAThread
        extends Thread {

        public STACalculator[] stacs;
        public int progress = 0;

        private final int n;
        private BufferedMovie submovie;
        private final int staDepth, staOffset;
        private final int[][] allFirstFrameIndices;
        private final int totalNumberOfSpikes;
        private final int firstTTL;
        private final double timePerSample;
        private final int spikesToCalculate;


        public STAThread(int n, int[][] allFirstFrameIndices, int totalNumberOfSpikes, int staDepth, int staOffset, int firstTTL, 
                double timePerSample, int spikesToCalculate, BufferedMovie submovie) throws IOException {

            this.staDepth = staDepth;
            this.staOffset = staOffset;
            this.n = n;
            this.allFirstFrameIndices = allFirstFrameIndices;
            this.totalNumberOfSpikes = totalNumberOfSpikes;
            this.firstTTL = firstTTL;
            this.timePerSample = timePerSample;
            this.spikesToCalculate = spikesToCalculate;
            this.submovie = submovie;
        }

        public void run() {
            stacs = new STACalculator[n];

            try {
                submovie.reset();

                Vision.getInstance().sendMessage("Creating STACalculators");
                for (int i = 0; i < n; i++)
                    stacs[i] = new STACalculator(staDepth, staOffset, firstTTL, timePerSample, submovie);

                Vision.getInstance().sendMessage("Processing Spikes");
                int[] spikesProcessed = new int[n];
                Arrays.fill(spikesProcessed, 0);
                int spikeIndex = 0;
                
                // Advance spikeIndex past any spikes that are too early, i.e. those that would put
                // FirstFrameIndex for the current STA Depth earlier than the first frame of the movie.
                for (int i = 0; i < n; i++) {
                    while (spikesProcessed[i] < allFirstFrameIndices[i].length && allFirstFrameIndices[i][spikesProcessed[i]] < 0) {
                        spikesProcessed[i]++;
                        spikeIndex++;
                    }
                }

                for (int f = 0; f < submovie.size() - staDepth; f++) {
                    for (int i = 0; i < n; i++) {
                        while (spikesProcessed[i] < spikesToCalculate &&
                               spikesProcessed[i] < allFirstFrameIndices[i].length && 
                               allFirstFrameIndices[i][spikesProcessed[i]] == f) {

                            stacs[i].addFramesAt(f);
                            spikesProcessed[i]++;
                            spikeIndex++;
                            progress = (int) (100L * spikeIndex / totalNumberOfSpikes);
                        }
                    }
                }

                // very important
                for (int i = 0; i < n; i++) stacs[i].finish();

            } catch (IOException e) {
                e.printStackTrace();
                return;
            }
        }
        
        public BufferedMovie getSubmovie() {
            return submovie;
        }
        
        
        public static boolean anySTAThreadsAlive(STAThread[] staThreads) {
            for (STAThread s : staThreads) {
                if (s.isAlive()) return true;
            }
            return false;
        }
        
        public static int averageProgress(STAThread[] staThreads) {
            int cumProg = 0;
            for (STAThread s : staThreads) cumProg += s.progress;
            return cumProg / staThreads.length;
        }
    };

    
    public static BufferedMovie[] buildMoviesThreaded(String movieFileName, String globalsFileName, 
            int numThreads, int staDepth) throws IOException {

        SubmovieThread[] submovieThreads = new SubmovieThread[numThreads];
        for (int i = 0; i < numThreads; i++) {
            submovieThreads[i] = new SubmovieThread(movieFileName, globalsFileName, i, numThreads);
            submovieThreads[i].start();
        }
        while (SubmovieThread.anySubmovieThreadsAlive(submovieThreads)) {
            try {
                Thread.sleep(100);
            } catch (InterruptedException e) {}
        }
        BufferedMovie[] movies = new BufferedMovie[numThreads];
        for (int i = 0; i < numThreads; i++) movies[i] = new BufferedMovie(submovieThreads[i].getWNMS(), staDepth);
        return movies;
    }
    

    public static STAThread[] calculateSTAsThreaded(
        final int n, final NeuronFile neuronFile, final int[][] allFirstFrameIndices, final int totalNumberOfSpikes, final int staDepth,
        final int staOffset, final int firstTTL, final double timePerSample, final int numThreads,
        int spikesToCalculate, BufferedMovie[] movies) throws IOException {
        
        Vision app = Vision.getInstance();
        
        STAThread[] staThreads = new STAThread[numThreads];
        for (int i = 0; i < numThreads; i++) {
            staThreads[i] = new STAThread(n, allFirstFrameIndices, totalNumberOfSpikes, staDepth, staOffset, firstTTL, timePerSample, 
                    spikesToCalculate, movies[i]);
        }
        for (int i = 0; i < numThreads; i++) staThreads[i].start();

        while (STAThread.anySTAThreadsAlive(staThreads)) {
            app.setProgress(STAThread.averageProgress(staThreads));
            try {
                Thread.sleep(100);
            } catch (InterruptedException e) {}
        }

        return staThreads;
    }

    
    public void setParameters(HashMap<String,String> parameters) {
        filePath = new File( (String) parameters.get("File_Path")).getAbsolutePath();
        rawMoviePath = parameters.get("Raw_Movie_Path");

        staDepth = Integer.parseInt( (String) parameters.get("STA Depth"));
        staOffset = Integer.parseInt((String) parameters.get("STA Offset"));
        atOnce = Integer.parseInt(parameters.get("STAs Calculated At Once"));
        spikesToCalculate = Integer.parseInt(parameters.get("Spikes To Calculate"));
        if (spikesToCalculate <= 0) spikesToCalculate = Integer.MAX_VALUE;
        
        numThreads = Integer.parseInt(parameters.get("nThreads"));
        numParts = parameters.containsKey("nParts")  ? Integer.parseInt(parameters.get("nParts"))  : 1;
        partNum  = parameters.containsKey("partNum") ? Integer.parseInt(parameters.get("partNum")) : 0;
        if (parameters.containsKey("proportions")) proportions = parseProportions(parameters.get("proportions"));
        
        calculateSTV = Boolean.valueOf( (String) parameters.get("Calculate STV")).booleanValue();
        resumeCalculation = Boolean.valueOf((String) parameters.get("Resume Calculation")).booleanValue();
    }
    
    private float[] parseProportions(String s) {
        if (s.length() == 0) return null;	

        String[] ss = s.split(",");
        float[] r = new float[ss.length];
        for (int i = 0; i < r.length; i++) r[i] = Float.parseFloat(ss[i]);
        
        float tot = 0;
        for (int i = 0; i < r.length; i++) tot += r[i];
        if (tot != 1.0f) System.err.printf("Warning: proportions add up to %f instead of 1\n", tot);
        
        return r;
    }
}