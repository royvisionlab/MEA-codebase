package edu.ucsc.neurobiology.vision.io;

import java.lang.reflect.Field;
import java.lang.reflect.Modifier;

import edu.ucsc.neurobiology.vision.anf.ElectrodeUsage;


/**
 * This class implements the header for a large number of high level IO classes in Vision.
 * Only int, float, double and Enum fields are supported.
 * See CVS/formats/VisionHeaderFormat.txt for a description of the fields.
 * 
 * NeuronFile uses a bastardization of this header.  It does not have enough room
 * so only some of the fields are set.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class VisionHeader implements Cloneable {

    //Public final variables are not saved to the file.
    public final static int NOISE_UNWHITENED = -2; 
    public final static int NOISE_WHITENED = 1;

    //Private final variables are saved.
    private static final int SIZE;

    static {
        int size = 0;

        Field[] f = VisionHeader.class.getFields();
        for (int i = 0; i < f.length; i++) {
            if ((f[i].getModifiers() & (Modifier.FINAL + Modifier.PUBLIC)) 
                    == (Modifier.FINAL + Modifier.PUBLIC)) {
                // Do nothing
            } else if (f[i].getType().isEnum()) {
                size += 4;
            } else if (f[i].getType().equals(int.class)) {
                size += 4;
            } else if (f[i].getType().equals(float.class)) {
                size += 4;
            } else if (f[i].getType().equals(double.class)) {
                size += 8;
            }
        }

        SIZE = size;
    }


    // general params
    public int magic;
    public int headerVersion; // -- shlens
    public int version;
    public int time;

    // raw data params
    public int arrayID;
    public int nSamples;
    public int samplingFrequency;

    // spike finding params
    public double meanTimeConstant;
    public double threshold;

    // covariance params
    public ElectrodeUsage electrodeUsage;
    public int minCovarianceSpikes;
    public int maxCovarianceSpikes;
    public int nlPoints;
    public int nrPoints;
    public double minimizationError;

    // projections params
    public int nDimensions;

    // clustering params
    public int binsPerDimension;
    public double clusteringSignificance;
    public int densityClusteringMaxSpikeLoss;
    public int minClusters;
    public int maxClusters;
    public int nEMSpikes;
    public int minEMIterations;
    public int maxEMIterations;
    public double emLikelihoodDelta;

    // neuron cleaning params
    public int minNeuronSpikes;
    public double acfT1;
    public double acfT2;
    public double maxContamination;
    public double coincidenceTime;
    public double maxCorrelation;
    public int removeDuplicates; // boolean if removed duplicates

    // EI imaging parameters
    public int nlPointsEI;
    public int nrPointsEI;
    public int maxEISpikes;

    // SNF parameters
    public double minThreshold;
    public int maxElectrodePatternSize;

    public int covarianceType;  //either NOISE_UNWHITENED or NOISE_WHITENED
    
    public int visionVersion;

    // unused space, 32 bytes
    public int unused3;
    public int unused4;
    public int unused5;
    public int unused6;
    public int unused7;
    public int unused8;


    /**
     * Initializes the header with -1 or null fields, meaning unset values.
     */
    public VisionHeader() {
        Field[] f = this.getClass().getFields();
        for (int i = 0; i < f.length; i++) {
            try {

                if((f[i].getModifiers() & (Modifier.FINAL + Modifier.PUBLIC)) 
                        == (Modifier.FINAL + Modifier.PUBLIC)) {
                    // Do nothing
                } else if (f[i].getType().isEnum()) {
                    f[i].set(this, null);
                } else if (f[i].getType().equals(int.class)) {
                    f[i].setInt(this, -2);
                } else if (f[i].getType().equals(float.class)) {
                    f[i].setFloat(this, -2);
                } else if (f[i].getType().equals(double.class)) {
                    f[i].setDouble(this, -2);
                }  
            } catch (IllegalAccessException ex) {
                throw new Error("Cannot access field " + f[i]);
            } catch (IllegalArgumentException ex) {
                throw new Error("Cannot access field " + f[i]);
            }
        }
    }


    public void initializeWithGarbage() {
        Field[] f = this.getClass().getFields();
        for (int i = 0; i < f.length; i++) {
            try {
                if ((f[i].getModifiers() & (Modifier.FINAL + Modifier.PUBLIC)) 
                        == (Modifier.FINAL + Modifier.PUBLIC)) {
                    // Do nothing
                } else if (f[i].getType().isEnum()) {
                    f[i].set(this, f[i].getType().getEnumConstants()[0]);
                } else if (f[i].getType().equals(int.class)) {
                    f[i].setInt(this, 0);
                } else if (f[i].getType().equals(float.class)) {
                    f[i].setFloat(this, 0);
                } else if (f[i].getType().equals(double.class)) {
                    f[i].setDouble(this, 0);
                }
            } catch (IllegalAccessException ex) {
                throw new Error("Cannot access field " + f[i]);
            } catch (IllegalArgumentException ex) {
                throw new Error("Cannot access field " + f[i]);
            }
        }
    }


    public int getFieldIndex(String fieldName) {
        Field[] f = this.getClass().getFields();
        for (int i = 0; i < f.length; i++) {
            if (f[i].getName().equals(fieldName)) {
                return i;
            }
        }
        return -1;
    }





    public VisionHeader(VisionHeader h) {
        if (h == null) {
            throw new NullPointerException("The given header is null");
        }

        Field[] f1 = this.getClass().getFields();
        Field[] f2 = h.getClass().getFields();
        Field[] f = f1.length < f2.length ? f1 : f2;
        for (int i = 0; i < f.length; i++) {
            try {
                Object o = f[i].get(h);
                if(!((f[i].getModifiers() & (Modifier.FINAL + Modifier.PUBLIC)) 
                        == (Modifier.FINAL + Modifier.PUBLIC))) {	
                    f[i].set(this, o);
                }
            } catch (IllegalAccessException ex) {
                throw new Error("Cannot access field " + f[i]);
            } catch (IllegalArgumentException ex) {
                //ex.printStackTrace();
                throw new Error("Cannot access field " + f[i] + ": " + ex.getMessage());
            }
        }
    }


    public VisionHeader(VisionHeader h, int junk) {
        // Shlens's implementation to work with CellFinder
        // Copy all fields in THIS (instance) that are available in h over to THIS (instance)
        if (h == null) throw new NullPointerException("The given header is null");

        Field[] theseFields = this.getClass().getFields();

        for (int i = 0; i < theseFields.length; i++) {
            try {
                Object o = theseFields[i].get(h);
                if(!((theseFields[i].getModifiers() & (Modifier.FINAL + Modifier.PUBLIC)) 
                        == (Modifier.FINAL + Modifier.PUBLIC))) {
                    if (o != null) {
                        theseFields[i].set(this, o);
                    }
                }
            } catch (IllegalAccessException ex) {
                throw new Error("Cannot access field " + theseFields[i]);
            } catch (IllegalArgumentException ex) {
                //ex.printStackTrace();
                //System.out.println("unable to find" + theseFields[i]);
                //throw new Error("Cannot access field " + theseFields[i] + ": " + ex.getMessage());
            }
        }
    }


    public Object clone() {
        try {
            return super.clone();
        } catch (CloneNotSupportedException ex) {
            return null;
        }
    }


    public static boolean equals(VisionHeader s1, VisionHeader s2) {
        if (s1.arrayID != s2.arrayID) return false;
        if (s1.electrodeUsage != s2.electrodeUsage) return false;
        if (s1.nDimensions != s2.nDimensions) return false;
        if (s1.nlPoints != s2.nlPoints) return false;
        if (s1.nrPoints != s2.nrPoints) return false;

        // compatible model !
        return true;
    }


    public static final int size() {
        return SIZE;
    }

}