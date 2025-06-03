package edu.ucsc.neurobiology.vision.io;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.util.*;


/**
 * An interface to the model file used in PCANFClustering.
 * The first structure in the file is the VisionHeader (look at this class).
 * Then for every neuron extraction (i.e. electrode) there is a "cell" (see CellFile)
 * containing an instance of the ClusteringModelFile.Model class
 * (see ClusteringModelFile.Model for details).
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ClusteringModelFile extends CellFile {
    private static final boolean DEBUG = false;
    private static final int MAGIC = 0xFAFAFA;
    private static final int EM_MODEL = 0;
    private static final int MANUAL_MODEL = 1;
    
    public ClusteringModelFile(String fileName, int nSlots, VisionHeader header) throws IOException {
        super(fileName, nSlots, createUserHeader(header));
    }
    
    public ClusteringModelFile(String fileName) throws IOException {
        super(fileName);
    }
    
    
    private static byte[] createUserHeader(VisionHeader header) throws IOException {
        header.magic = MAGIC;
        header.version = 1;
        
        ByteArrayOutputStream bos = new ByteArrayOutputStream();
        DataOutputStream dos = new DataOutputStream(bos);
        IOUtil.writePublicFields(header, dos);
        dos.close();
        bos.close();
        return bos.toByteArray();
    }


    synchronized public VisionHeader getUserHeader() throws IOException {
        byte[] byteHeader = super.readUserHeader();
        ByteArrayInputStream bos = new ByteArrayInputStream(byteHeader);
        DataInputStream dis = new DataInputStream(bos);
        VisionHeader header = new VisionHeader();
        IOUtil.readPublicFields(header, dis);
        dis.close();
        bos.close();
        return header;
    }
    
    
    synchronized public int getNextExtractionID() {
        int[] idList = getIDList();
        if (idList.length == 0) return 1;
        return idList[idList.length - 1] + 1;
    }
    
    
    synchronized public int[] getNeuronsForCleanning(int cleaningLevel) throws IOException {
        IntegerList id = new IntegerList();
        int[] extractionID = getIDList();
        for (int i = 0; i < extractionID.length; i++) {
            Model m = getNeuronExtraction(extractionID[i]);
            if (m.cleaningLevel == cleaningLevel) {
                for (int j = 0; j < m.neuronID.length; j++) {
                    id.add(m.neuronID[j]);
                }
            }
        }
        return id.toArray();
    }


    /**
     * This class contains a record of all the information relevant to a neuron extraction
     * (done on a specific electrode). The information stored includes the neuron
     * indices and IDs, the number of dimensions and number of clusters used, the PCA
     * eigenvectors used, the list of electrodes used, etc. Subclasses add more specific
     * information depending on the exact neuron identification scheme used
     * (automatic, EM or manual, a user selection). See subclasses for details.
     *
     * @author Dumitru Petrusca, University of California, Santa Cruz
     */
    public static abstract class Model implements Cloneable {
        public int extractionID;
        public int[] neuronIndex;
        public int[] neuronID;
        public int cleaningLevel;
        public int[] electrodes;
        public double threshold;
        public int nDimensions;
        public double[][] eigenvectors;
        public int nClusters;
        
        synchronized protected void write(DataOutput dis) throws IOException {
            dis.writeInt(neuronIndex.length);
            for (int i = 0; i < neuronIndex.length; i++) {
                dis.writeInt(neuronIndex[i]);
                dis.writeInt(neuronID[i]);
            }
            dis.writeInt(cleaningLevel);

            // write the electrode pattern
            dis.writeInt(electrodes.length);
            for (int i = 0; i < electrodes.length; i++) {
                dis.writeInt(electrodes[i]);
            }

            dis.writeDouble(threshold);

            // write the eigenvectors
            dis.writeInt(eigenvectors.length);
            final int vectorLength = eigenvectors[0].length;
            dis.writeInt(vectorLength);
            for (int i = 0; i < eigenvectors.length; i++) {
                for (int j = 0; j < vectorLength; j++) {
                    dis.writeDouble(eigenvectors[i][j]);
                }
            }
            dis.writeInt(nDimensions);
            dis.writeInt(nClusters);
        }


        synchronized protected void read(DataInput dis) throws IOException {
            int nNeurons = dis.readInt();
            neuronIndex = new int[nNeurons];
            neuronID = new int[nNeurons];
            for (int i = 0; i < nNeurons; i++) {
                neuronIndex[i] = dis.readInt();
                neuronID[i] = dis.readInt();
            }
            cleaningLevel = dis.readInt();

            // write the electrode pattern
            final int nElectrodes = dis.readInt();
            electrodes = new int[nElectrodes];
            for (int i = 0; i < electrodes.length; i++) {
                electrodes[i] = dis.readInt();
            }

            threshold = dis.readDouble();

            // write the eigenvectors
            final int nEigenvectrors = dis.readInt();
            final int vectorLength = dis.readInt();
            eigenvectors = new double[nEigenvectrors][vectorLength];
            for (int i = 0; i < nEigenvectrors; i++) {
                for (int j = 0; j < vectorLength; j++) {
                    eigenvectors[i][j] = dis.readDouble();
                }
            }

            // write the EM model parameters
            nDimensions = dis.readInt();
            nClusters = dis.readInt();
        }


        public void dump() {
            System.out.println("------- " + getClass().getName() + " for " + neuronID);
            System.out.println("Neuron Index " + neuronIndex);
            System.out.println("Cleanning Level " + cleaningLevel);
            IOUtil.printArray(electrodes);
            System.out.println("Threshold " + threshold);
            System.out.println("nDimensions " + nDimensions);
            System.out.println("Eigenvectors");
            for (int i = 0; i < nDimensions; i++) {
                IOUtil.printArray(eigenvectors[i]);
            }
            System.out.println("nClusters " + nClusters);
        }
    }


    /**
     * Sublass that adds more information to the model in case of automatic clustering (EM).
     * The specific info is the number of Gaussian components, their probabilities, their
     * locations and their sizes.
     *
     * @author Dumitru Petrusca, University of California, Santa Cruz
     */
    public static class EMModel
        extends Model {
        public int nGaussians;
        public double[] probability;
        public double[][] means;
        public double[][] covariances;


        synchronized public void write(DataOutput dis) throws IOException {
            super.write(dis);

            // write the EM model parameters
            dis.writeInt(nGaussians);

            for (int neuron = 0; neuron < nGaussians; neuron++) {
                dis.writeDouble(probability[neuron]);

                for (int d = 0; d < nDimensions; d++)
                    dis.writeDouble(means[neuron][d]);

                for (int d = 0; d < nDimensions; d++)
                    dis.writeDouble(covariances[neuron][d]);
            }
        }


        synchronized protected void read(DataInput dis) throws IOException {
            super.read(dis);

            nGaussians = dis.readInt();

            probability = new double[nGaussians];
            means = new double[nGaussians][nDimensions];
            covariances = new double[nGaussians][nDimensions];

            for (int neuron = 0; neuron < probability.length; neuron++) {
                probability[neuron] = dis.readDouble();

                for (int d = 0; d < nDimensions; d++)
                    means[neuron][d] = dis.readDouble();

                for (int d = 0; d < nDimensions; d++)
                    covariances[neuron][d] = dis.readDouble();
            }
        }


        public void dump() {
            super.dump();

            System.out.println("nGaussians " + nGaussians);
            System.out.print("Probabilities: ");
            IOUtil.printArray(probability);
            for (int i = 0; i < probability.length; i++) {
                System.out.print("Means: ");
                IOUtil.printArray(means[i]);
                System.out.print("Sigmas: ");
                IOUtil.printArray(covariances[i]);
            }
            System.out.println("-------");
        }
    }


    /**
     * Subclass that adds more information to the model in case of manual clustering.
     * The specific info is the selection region entered by the user for every cluster and
     * the dimensions in which these selections were made.
     *
     * @author Dumitru Petrusca, University of California, Santa Cruz
     */
    public static class ManualModel
        extends Model {

        public ArrayList<double[]> selections;
        public int[] dimension1;
        public int[] dimension2;


        synchronized public void write(DataOutput dis) throws IOException {
            super.write(dis);

            // write the Box model parameters

            for (int neuron = 0; neuron < nClusters - 1; neuron++) {
                // write the dimensions
                dis.writeInt(dimension1[neuron]);
                dis.writeInt(dimension2[neuron]);

                // write the selection
                double[] s = (double[]) selections.get(neuron);
                dis.writeInt(s.length);
                for (int i = 0; i < s.length; i++) {
                    dis.writeDouble(s[i]);
                }
            }
        }


        synchronized protected void read(DataInput dis) throws IOException {
            super.read(dis);

            dimension1 = new int[nClusters - 1];
            dimension2 = new int[nClusters - 1];
            selections = new ArrayList<double[]>();

            for (int neuron = 0; neuron < nClusters - 1; neuron++) {
                dimension1[neuron] = dis.readInt();
                dimension2[neuron] = dis.readInt();

                int n = dis.readInt();
                double[] s = new double[n];
                for (int i = 0; i < n; i++) {
                    s[i] = dis.readDouble();
                }
                selections.add(s);
            }
        }


        public void dump() {
            super.dump();

            System.out.println("Dimensions: " + dimension1 + ", " + dimension2);
            for (int i = 0; i < selections.size(); i++) {
                double[] box = (double[]) selections.get(i);
                System.out.println("Location: " + box[0] + ", " + box[1]);
                System.out.println("Size: " + box[2] + ", " + box[3]);
            }
            System.out.println("-------");
        }
    }


    synchronized public void addExtraction(Model model) throws IOException {
//        System.out.println("Add Neuron " + model.extractionID);
        ByteArrayOutputStream bos = new ByteArrayOutputStream();
        DataOutputStream dis = new DataOutputStream(bos);

        if (model instanceof EMModel) dis.writeInt(EM_MODEL);
        else dis.writeInt(MANUAL_MODEL);
        model.write(dis);

        dis.close();
        super.addCell(model.extractionID, bos.toByteArray());
        if (DEBUG) model.dump();
    }


    synchronized public Model getNeuronExtraction(int extractionID) throws IOException {
        byte[] cell = super.getCell(extractionID);
        if (cell == null) return null;
        
        ByteArrayInputStream bos = new ByteArrayInputStream(cell);
        DataInputStream dis = new DataInputStream(bos);
        
        Model model;
        if (dis.readInt() == 0) model = new EMModel();
        else model = new ManualModel();
        model.read(dis);
        model.extractionID = extractionID;
        
        dis.close();
        if (DEBUG) model.dump();
        return model;
    }


    synchronized private static ClusteringModelFile.EMModel getModel() throws IOException {
        int nNeurons = 5;

        ClusteringModelFile.EMModel m = new ClusteringModelFile.EMModel();
        m.extractionID = 1;
        m.neuronIndex = new int[nNeurons];
        m.neuronID = new int[nNeurons];
        for (int neuron = 0; neuron < nNeurons; neuron++) {
            m.neuronIndex[neuron] = neuron;
            m.neuronID[neuron] = neuron;
        }
        m.cleaningLevel = 1;
        m.electrodes = new int[7];
        m.threshold = 0;
        m.nDimensions = 5;
        m.eigenvectors = new double[5][24];
        m.nClusters = 5;
        m.nGaussians = 5;
        m.probability = new double[m.nGaussians];
        m.means = new double[m.nGaussians][];
        m.covariances = new double[m.nGaussians][];
        for (int i = 0; i < m.nGaussians; i++) {
            m.probability[i] = i;
            m.means[i] = new double[5];
            m.covariances[i] = new double[5];
        }

        return m;
    }


    public static void profile() throws IOException {
        VisionHeader header = new VisionHeader();
        header.nlPoints = 0;
        header.nrPoints = 0;
        header.nlPointsEI = 10;
        header.nrPointsEI = 20;
        header.acfT1 = 0.5;
        header.acfT2 = 1.0;
        header.nDimensions = 5;
        header.maxClusters = 10;
        //header.minimizationInterval = -1;
        header.minimizationError = -1;
        header.coincidenceTime = -1;
        header.maxElectrodePatternSize = -1;
        header.minThreshold = -1;

        new File("1.model").delete();
        final ClusteringModelFile modelFile = new ClusteringModelFile("1.model", 50000, header);
        final Model m = getModel();

        Thread t = new Thread() {
            public void run() {
                try {
                    for (int i = 0; i < 50000; i++) {
                        modelFile.addExtraction(m);
                    }
                } catch (IOException ex) {
                    ex.printStackTrace();
                }
            }
        };
        t.start();
    }
}