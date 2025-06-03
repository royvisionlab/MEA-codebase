package edu.ucsc.neurobiology.vision.util;

import java.util.HashMap;
import java.util.Map;

import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.io.chunk.*;

/**
 * This class contains some global-scope parameters.
 *
 * @author Jon Shlens, Salk Institute, La Jolla<br>
 *         Dumitru Petrusca, University of California, Santa Cruz
 */
public final class VisionParams {
    
    // no objects of this type may be created
    private VisionParams() {}
    
    public static final int VERSION = 7003051;  //short for X.XXX.XXX.
    
    /**
     * Maximum number of clusters allowed on an electrode. This number is crucial
     * for determining the neuronID for every single cluster. Change this number with care.
     */
    public static final int maxClusters = 15;
    
    public static final String VISION_OUTPUT_FILE = "vision-output.txt";
    
    // maximum number of slots to store neurons in the model file
    public static final int MAX_MODEL_SLOTS = 5000;
    
    // maximum number of slots to store neurons in the neurons file
    public static final int NEURONS_HEADER_CAPACITY = 50000;
    
    public final static String BIN_FILE_EXTENSION_512 = ".bin";
    public final static String BZ2_FILE_EXTENSION_512 = ".bz2";
    public final static String HEADER_FILE_EXTENSION =  ".header";
    public final static String REM_FILE_EXTENSION_512 = ".rem";
    public final static String NOISE_FILE_EXTENSION =   ".noise";
    public final static String NEURON_FILE_EXTENSION =  ".neurons";
    public final static String MODEL_FILE_EXTENSION =   ".model";
    public final static String PROJ_FILE_EXTENSION =    ".prj";
    public final static String SPIKES_FILE_EXTENSION =  ".spikes";
    public final static String STA_FILE_EXTENSION =     ".sta";
    public final static String STE_FILE_EXTENSION =     ".ste";
    public final static String GLOBALS_FILE_EXTENSION = ".globals";
    public final static String PARAMS_FILE_EXTENSION =  ".params";
    public final static String EI_FILE_EXTENSION =      ".ei";
    public final static String MOVIE_FILE_EXTENSION =   ".movie";
    
    public final static Map<String,Class<?>> EXTENSION_MAP = new HashMap<String,Class<?>>();
    static {
        EXTENSION_MAP.put(NEURON_FILE_EXTENSION,  NeuronFile.class);
        EXTENSION_MAP.put(BIN_FILE_EXTENSION_512, RawDataFile.class);
        EXTENSION_MAP.put(EI_FILE_EXTENSION,      PhysiologicalImagingFile.class);
        EXTENSION_MAP.put(GLOBALS_FILE_EXTENSION, GlobalsFile.class);
        EXTENSION_MAP.put(MOVIE_FILE_EXTENSION,   WhiteNoiseMovieFile.class);
    }
    
    // In milliseconds
    public static final double ACFT1 = 0.5;
    public static final double ACFT2 = 1.0;
    
    public static final double SAMPLES_PER_MILLISECOND = 20;
    
    // TTL Pulse defaults
    public final static double DEFAULT_MONITOR_FREQUENCY = 120.0;
    public final static int DEFAULT_FRAMES_PER_TTL = 100;
    
    // bonus PCA Clustering parameter.  Takes the number given
    // by EM spikes and further reduces it during the convergence determination.
    // extracted to make it more visible.
    // larger numbers make the convergence determination less accurate
    // and faster.
    public final static int LIKELIHOOD_CALC_STRIDE = 3;
    
    
    public static String versionString(int version) {
        return "" + version/1000000 + "." + (version/1000) % 1000 + "." + version % 100;
    }
    
    public static String versionString() {
        return versionString(VERSION);
    }
    
}
