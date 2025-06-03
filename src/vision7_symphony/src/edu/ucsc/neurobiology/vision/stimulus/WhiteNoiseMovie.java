package edu.ucsc.neurobiology.vision.stimulus;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.io.chunk.*;
import edu.ucsc.neurobiology.vision.stimulus.FrameGenerator.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * This class encapsulates the writing and reading of a white noise movie file.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class WhiteNoiseMovie
/*extends Tag*/ implements AdvancedMovie {

    public long startSeed;
    private int nFrames;
    private int width;
    private int height;
    private double stixelWidth;
    private double stixelHeight;
    private double refreshTime;
    public ColorType colorType;
    public MovieType noiseType;
    public RandomNumberGenerator rng;
    public float contrastValue;
    public int[] droppedFrameIndex;

    public AbstractRandomFrameGenerator frameGenerator;
    public long[] seeds;
    private STAFrame generalFrame;
    private float[] commonFrame;
    private FrameEncoding frameEncoding = FrameEncoding.PACKED_ENCODING;
    private boolean calculateSTV = true;
    private boolean sparse = false;
    private float probability;
    

    public WhiteNoiseMovie(
            int width, int height, double stixelWidth, double stixelHeight, double refreshTime,
            int nFrames, long startSeed, ColorType colorType, MovieType noiseType, boolean sparse, float probability, 
            RandomNumberGenerator rng, float contrastValue, int[] droppedFrameIndex,
            long[] seeds, boolean generateSeeds) throws IOException {

        this.width = width;
        this.height = height;
        this.stixelWidth = stixelHeight;
        this.stixelHeight = stixelHeight;
        this.refreshTime = refreshTime;
        this.nFrames = nFrames;
        this.startSeed = startSeed;
        this.colorType = colorType;
        this.noiseType = noiseType;
        this.rng = rng;
        this.contrastValue = contrastValue;
        this.droppedFrameIndex = droppedFrameIndex;
        this.sparse = sparse;
        this.probability = probability;
        
        init();

        // obtain the seeds
        if (seeds != null) {
            this.seeds = seeds;
        } else if (generateSeeds) {
            generateSeeds();
        }
    }
    
    public WhiteNoiseMovie(
            int width, int height, double stixelWidth, double stixelHeight, double refreshTime,
            int nFrames, long startSeed, ColorType colorType, MovieType noiseType,
            RandomNumberGenerator rng, float contrastValue, int[] droppedFrameIndex,
            long[] seeds, boolean generateSeeds) throws IOException {
        this(width, height, stixelWidth, stixelHeight, refreshTime, nFrames, startSeed, colorType, 
                noiseType, false, 0, rng, contrastValue, droppedFrameIndex, seeds, generateSeeds);
    }
    
    
    public WhiteNoiseMovie(
            int width, int height, double stixelWidth, double stixelHeight, double refreshTime,
            int nFrames, long startSeed, ColorType colorType, MovieType noiseType,
            RandomNumberGenerator rng, float contrastValue, int[] droppedFrameIndex,
            long[] seeds) throws IOException {

        this(width, height, stixelWidth, stixelHeight, refreshTime, nFrames, startSeed, colorType, 
                noiseType, rng, contrastValue, droppedFrameIndex, seeds, true);
        // true means generateSeeds, but if the seeds passed isn't null then that will override and generateSeeds won't be called.
    }

    
    public WhiteNoiseMovie(String movieFileName, String globalsFileName, boolean generateSeeds) {
        loadRawMovieFileFields(movieFileName);		
        loadGlobalsFileFields(globalsFileName);
        init();
        if (generateSeeds) generateSeeds();
    }

    
    public WhiteNoiseMovie(String movieFileName, String globalsFileName) {
        this(movieFileName, globalsFileName, true);
    }

    
    
    void init() {
        if ((colorType != ColorType.DEPENDENT) &&
            (colorType != ColorType.INDEPENDENT) &&
            (colorType != ColorType.SEPARATED))
            throw new IllegalArgumentException("Illegal color type");
        
        switch (noiseType) {
        case BINARY_MOVIE:
            frameGenerator = new BinaryFrameGenerator(width, height, contrastValue,
                    rng, startSeed, colorType);
            break;

        case GAUSSIAN_MOVIE:
        	System.out.println("here");
            frameGenerator = new GaussFrameGenerator(width, height, 
                    new float[] {contrastValue, contrastValue, contrastValue},
                    rng, startSeed, colorType);
            break;
            
        default:
            throw new IllegalArgumentException("Unknown movie type " + noiseType);
        }
        if (sparse) convertToSparse();
        
        commonFrame = new float[width * height * 3];
        float[] colorFrame = new float[width * height * 3];
        generalFrame = new STAFrame(width, height, stixelWidth, stixelHeight, colorFrame);
    }

    
    /**
     * Wrap existing frameGenerator into a sparse generator
     */
    private void convertToSparse() {
        SparseFrameGenerator sfg = new SparseFrameGenerator(probability, (AbstractConstantSeedFrameGenerator) frameGenerator);
        frameGenerator = sfg;
    }
    
    
    void generateSeeds() {
        seeds = new long[nFrames];
        long currentSeed = frameGenerator.getSeed();
        int dropIndex = 0;
        int frameIndex = 0;
        
        try {
            for (int i = 0; i < nFrames; i++) {
                seeds[i] = currentSeed;

                // Change the seed to the next frame, unless we have dropped frames to check for
                if (droppedFrameIndex.length > 0 && dropIndex < droppedFrameIndex.length) {
                    // Change the seed to the next frame, unless this frame was dropped
                    if (droppedFrameIndex[dropIndex] == frameIndex) {
                        dropIndex++;
                    } else {
                        // Frame wasn't dropped: we can advance the frame index
                        frameGenerator.advanceFrameSeed(generalFrame.getBuffer());
                        currentSeed = frameGenerator.getSeed();
                        frameIndex++;
                    }
                } else {
                    // No need to worry about the frame index if there were no dropped frames
                    frameGenerator.advanceFrameSeed(generalFrame.getBuffer());
                    currentSeed = frameGenerator.getSeed();
                }
            }
        } catch(Exception e) {
            Vision.reportException("Unable to generate seeds.", e);
        }
    }
    
    
    void loadRawMovieFileFields(String movieFileName) {
        try {
            WhiteNoiseMovieFile movieFile = new WhiteNoiseMovieFile(movieFileName, ChunkFile.READ);
            WhiteNoiseMovieFile.WhiteMovieParams movieParams = movieFile.getWhiteMovieParams();
            width         = movieParams.width;
            height        = movieParams.height;
            startSeed     = movieParams.seed;
            colorType     = movieParams.colorType;
            noiseType     = movieParams.movieType;
            rng           = movieParams.rng;
            contrastValue = (float) movieParams.contrastValue;		
            sparse		  = movieParams.sparse;
            probability   = movieParams.probability;
        } catch (IOException e) {
            Vision.reportException(".movie file cannot be loaded.", e);
        }
    }
    
    
    void loadGlobalsFileFields(String globalsFileName) {
        try {
            GlobalsFile globalsFile = new GlobalsFile(globalsFileName, ChunkFile.READ);
            GlobalsFile.RunTimeMovieParams runTime = globalsFile.getRunTimeMovieParams();
            refreshTime       = runTime.refreshPeriod;
            droppedFrameIndex = runTime.droppedFrames;
            nFrames           = runTime.nFramesRequired;
            stixelWidth       = runTime.micronsPerStixelX;
            stixelHeight      = runTime.micronsPerStixelY;
        } catch (IOException e) {
            Vision.reportException(".globals file cannot be loaded.", e);
        }
    }


    void resize(int width, int height) {
        this.width = width;
        this.height = height;
        init(); // Regenerate frameGenerator and STA frames.
    }
    
    
    public UpdatingFrame createCompatibleFrame() {
        switch (noiseType) {
            case BINARY_MOVIE:
                return new BinaryUpdatingFrame(width, height, colorType);

            case GAUSSIAN_MOVIE:
                return new GaussUpdatingFrame(width, height, calculateSTV);
            
            default:
                throw new IllegalArgumentException("Unknown movie type.");
        }
    }


    public void setFrameEncoding(FrameEncoding frameEncoding) {
        this.frameEncoding = frameEncoding;
    }


    public Object getEncodedFrame(int frameIndex, Object frame) throws IOException {
        if (frame == null) {
            switch (frameEncoding) {
            case PACKED_ENCODING:
                frame = frameGenerator.createFrameBuffer();
                break;

            case FLOAT_ARRAY_ENCODING:
                frame = new float[width * height * 3];
                break;
            }
        }
        
        frameGenerator.nextFrame(seeds[frameIndex], frame);
        return frame;
    }


    public ImageFrame getFrame(int frameIndex) throws IOException {
        frameGenerator.nextFrame(seeds[frameIndex], commonFrame);
        generalFrame.setBuffer(commonFrame);
        return generalFrame;
    }


    public String getDescription() {
        return "White Noise Movie(" + size() + ")";
    }


    public double getRefreshTime() {
        return refreshTime;
    }


    public int getWidth() {
        return width;
    }


    public int getHeight() {
        return height;
    }


    @Override
    public int size() {
        return nFrames;
    }


    public String toString() {
        return "startSeed = " + startSeed + "\n" + 
        "width = " + width + "\n" + 
        "height = " + height + "\n" + 
        "refreshTime = " + refreshTime + "\n" + 
        "nFrames = " + nFrames + "\n" + 
        "colorType = " + colorType;
    }


    /**
     * Creates a new instance and then reads the relevant information
     * from the stream.
     */
//	public Tag read(int tagId, NeuroInputStream input, int len) throws IOException {
//	long seed = input.readLong();
//	int n = input.readInt();
//	int w = input.readInt();
//	int h = input.readInt();
//	double refresh = input.readDouble();
//	double pixSize = input.readDouble();
//	ColorType _colorType = ColorType.values()[input.readInt()];
////	System.out.println(_colorType);
//	MovieType _noiseType = MovieType.values()[input.readInt()];
////	System.out.println(_noiseType);
//	RandomNumberGenerator _rng = RandomNumberGenerator.values()[input.readInt()];
////	System.out.println(_rng);
//	float _constrastValue = input.readFloat();

//	int nDropped = input.readInt();
//	int[] droppedFrameInd = new int[nDropped];
//	for (int i = 0; i < droppedFrameInd.length; i++) {
//	droppedFrameInd[i] = input.readInt();
//	}

//	int nSeeds = input.readInt();
//	long[] seeds = new long[nSeeds];
//	for (int i = 0; i < seeds.length; i++) {
//	seeds[i] = input.readLong();
//	}

//	return new WhiteNoiseMovie(
//	w, h, pixSize, refresh, n, seed, _colorType, _noiseType, _rng,
//	_constrastValue,
//	droppedFrameInd, seeds);
//	}


    /**
     * Write out the electrode ID and the amplitude in addition to
     * first writing out the parent's data.
     */
//	public void write(int tagId, NeuroOutputStream output) throws IOException {
//	output.writeLong(startSeed);
//	output.writeInt(nFrames);
//	output.writeInt(width);
//	output.writeInt(height);
//	output.writeDouble(refreshTime);
//	output.writeDouble(pixelSize);
//	output.writeInt(colorType.ordinal());
//	output.writeInt(noiseType.ordinal());
//	output.writeInt(rng.ordinal());
//	output.writeFloat(contrastValue);

//	output.writeInt(droppedFrameIndex.length);
//	for (int i = 0; i < droppedFrameIndex.length; i++) {
//	output.writeInt(droppedFrameIndex[i]);
//	}

//	output.writeInt(seeds.length);
//	for (int i = 0; i < seeds.length; i++) {
//	output.writeLong(seeds[i]);
//	}
//	}


//	private static WhiteNoiseMovie load(String fileName) throws IOException {
        //   NeuroInputStream nis = new NeuroInputStream(new FileInputStream(fileName));
        //   WhiteNoiseMovie m = (WhiteNoiseMovie) nis.readTag();
        //   nis.close();
        //  return m;
//		return null;
//	}


    public WhiteNoiseMovie duplicate() throws IOException {
        return new WhiteNoiseMovie(width, height, stixelWidth, stixelHeight, refreshTime, nFrames,
                startSeed, colorType, noiseType, sparse, probability, rng, contrastValue,
                droppedFrameIndex, seeds, true); 
        // true means generateSeeds, but if the seeds passed isn't null then that will override and generateSeeds won't be called.
    }


//	public static final void createWhiteNoiseMovie(
//	String neuronFileName, int refreshInterval, int width, int height, int seed,
//	ColorType colorType, MovieType noiseType, RandomNumberGenerator rng,
//	double contrastSigma, double pixelSize) throws IOException {

////	System.err.println(colorType);

//	if (pixelSize == -1) {
//	pixelSize = 640 / width;
//	}

//	NeuronFile nf = new NeuronFile(neuronFileName);
//	int[] ttlTimes = nf.getTTLTimes();
//	ArrayList<Integer> droppedFramesList = new ArrayList<Integer> ();
//	double refreshPeriod = WhiteNoiseMovie.calculateRefreshPeriod(droppedFramesList,
//	ttlTimes, refreshInterval);

//	//The frame number of each dropped frame
//	int[] droppedFrameNumbers = new int[droppedFramesList.size()];

//	for (int i = 0; i < droppedFramesList.size(); i++) {
//	droppedFrameNumbers[i] = (int) (droppedFramesList.get(i));
//	}

//	// create the movie
//	WhiteNoiseMovie movie = new WhiteNoiseMovie(
//	width,
//	height,
//	pixelSize,
//	refreshPeriod,
//	( (ttlTimes.length + 1) * 100) / refreshInterval,
//	seed,
//	colorType,
//	noiseType,
//	rng,
//	(float) contrastSigma,
//	droppedFrameNumbers,
//	null
//	);

//	String fName = StringUtil.removeExtension(neuronFileName) +
//	VisionParams.MOVIE_FILE_EXTENSION;
//	NeuroOutputStream s = new NeuroOutputStream(new
//	FileOutputStream(fName));
//	s.writeTag(movie);
//	s.close();
//	}


    /**
     * Calculates the refresh period (in ms) and the dropped frames list (in movie frame number).
     * droppedFramesList must be newly initialized before this function is called.
     * 
     * @author Matthew Grivich, UCSC
     */
    public static double calculateRefreshPeriod(ArrayList<Integer> droppedFramesList,
            int[] ttlTimes, int interval, double monitorFrequency, int setFramesPerTTL) {
        double sum = 0;
        int nGoodTTL = 0;
        int partialDroppedFramesCounter = 0;
        
        for (int i = 1; i < ttlTimes.length; i++) {
            //Gives the number of frames, between two ttls.
            int nFramesPerTTL = (int) Math.round( (ttlTimes[i] - ttlTimes[i - 1]) /
                (1000.0 * VisionParams.SAMPLES_PER_MILLISECOND / monitorFrequency));
            
            if (nFramesPerTTL == setFramesPerTTL) {
                // if measurement of frames per ttl is correct (100), no frames
                // were dropped.  Use in average calculation
                sum += ttlTimes[i] - ttlTimes[i - 1];
                nGoodTTL++;

            } else if (nFramesPerTTL > setFramesPerTTL) {	
                //else, record approximate time of each dropped frame, in samples
                int currentDroppedFrames = nFramesPerTTL - setFramesPerTTL;
                for (int j = 0; j < currentDroppedFrames; j++) {
                    System.out.println("A stimulus machine frame was dropped.");
                    //prevents multiple dropped frames from having the same index.
                    int numberDropped = 0;  // In reality, appears to do nothing useful. Why is it reset every time?
                                            // Declaration in the wrong scope maybe?
                    partialDroppedFramesCounter++;

                    //Dropped frames are checked for in stimulus machine frames,
                    //While STAs are calculated in movie frames.
                    //If enough stimulus machine dropped frames have accumulated,
                    //Create a movie dropped frame.
                    if (partialDroppedFramesCounter == interval) {
                        System.out.println("A movie frame was dropped.");
                        partialDroppedFramesCounter = 0;
                        droppedFramesList.add(new Integer( (ttlTimes[i] + ttlTimes[i - 1]) /
                                2) + numberDropped++);
                    }
                }
                
            } else { 
                // this means nFramesPerTTL < VisionParams.nFramesPerTTL, a weird result
                System.out.println("Warning: only " + nFramesPerTTL + " frames between TTLs " + (i-1) + " and " + i + 
                     "; is monitorFrequency set too low?");
                // Since this is not handled more, this situation breaks the movie, 
                // So perhaps this should throw an Error?
            }
        }
        
        double refreshPeriod = sum * interval /
            (nGoodTTL * VisionParams.SAMPLES_PER_MILLISECOND * setFramesPerTTL);

        //Changes dropped frame list to return frame number of each dropped frame
        for (int i = 0; i < droppedFramesList.size(); i++) {
            droppedFramesList.set(i,
                    new Integer( (int) (droppedFramesList.get(i) /
                            VisionParams.SAMPLES_PER_MILLISECOND /
                            refreshPeriod)));
            System.out.println("Dropped frame at movie frame " + droppedFramesList.get(i) +
            ".");
        }
        
        return refreshPeriod;
    }
}
