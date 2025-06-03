package edu.ucsc.neurobiology.vision.io;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.Vision;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * An interface to the neuron file. For every neuron only the list of spike times is stored.
 * Each neuron is associated with a unique ID. The TTL times can also be retrieved.
 * One can also use this class to create neuron files. Look in the "anf" package
 * for uses of this class.
 * 
 * 
 * Note that this file does not have an actual VisionHeader.
 * There is some glue to make it looks like it does, but many 
 * of the fields are not set, because there is no room.
 * <br>
 * <br>
 * The class is thread safe.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 * @owner Matthew Grivich, The Salk Institute
 */
public class NeuronFile {


    /* =========================================================== */
    /* === Static Declarations =================================== */

    private static final boolean DEBUG = false;
    
    public static final int INT_VERSION = 32;
    public static final int SALK_VERSION = 33;
    public static final int DOUBLE_VERSION = 100;
    
    private static final int HEADER_LENGTH = 152;
    private static final int SLOT_LENGTH = 4 + 4 + 8;
    private static final int UNUSED_SLOT = Integer.MIN_VALUE;
    
    /* =========================================================== */
    /* === Static Functions ======================================= */

    // Note that the neuronID depends on the maximum clusters assigned to each electrode
    public static int getNeuronID(int electrode, int cluster) {
        // safety check
        if (cluster >= VisionParams.maxClusters) {
            System.out.println("NeuronID for cluster " + cluster +
                               " is invalid.");
            return -1;
        }

        int neuronID = ( (electrode - 1) * VisionParams.maxClusters) + cluster + 1;

        return neuronID;
    }


    public static int getNeuronIDCluster(int neuronID) {

        int cluster = (int) Math.IEEEremainder( ( (double) (neuronID - 1)),
                                               ( (double) VisionParams.
                                                maxClusters));

        // Math.IEEEremainder has some funny behavior. Please see the documentation
        // on this function to understand my work around.
        if (cluster < 0) {
            cluster = cluster + VisionParams.maxClusters;
        }
        return cluster;
    }


    public static int getNeuronIDElectrode(int neuronID) {

        int cluster = (int) Math.IEEEremainder( ( (double) (neuronID - 1)),
                                               ( (double) VisionParams.maxClusters));

        // Math.IEEEremainder has some funny behavior. Please see the documentation
        // on this function to understand my work around.
        if (cluster < 0) {
            cluster = cluster + VisionParams.maxClusters;
        }

        int electrode = (int) ( (neuronID - 1 - cluster) /
                               VisionParams.maxClusters) + 1;
        return electrode;
    }


    /* =========================================================== */
    /* =========================================================== */




    private final String fileName;
    private RandomAccessFile file;
    private final int valueLength;
    private byte[] buffer = new byte[10 * 1024 * 1024];

    private final int version;
    private final int headerCapacity;
    private final int nSamples;
    private final int samplingFrequency;
    private final long time;
    
    private final double maxContamination;
    private final int minNeuronSpikes, removeDuplicates, covarianceType, visionVersion;
    
    
    private final int[] neuronIDList;
    private final int[] electrodeList;
    private final long[] startLocationList;
    private int automaticNeuronID;
    private VisionHeader header;


    public NeuronFile(String fileName, VisionHeader header, int headerCapacity, int[] ttl) throws
        IOException {
        this(fileName, header, headerCapacity, ttl, ttl.length);
    }




    public NeuronFile(String fileName, VisionHeader header, int headerCapacity,
                      Object ttl, int nTTL) throws IOException {

        if (header.version == INT_VERSION) {
            valueLength = 4;
        } else if (header.version == DOUBLE_VERSION) {
            valueLength = 8;
        } else if (header.version == SALK_VERSION) {
            throw new Error(
                "The neuron file " + fileName + " has version " + SALK_VERSION +
                " used by EJ. Please modify the version to " + DOUBLE_VERSION);
        } else {
            throw new Error("Wrong version, ");
        }

        this.header = header;

        this.fileName = fileName;
        File f = new File(fileName);
        if (f.exists()) {
            f.delete();
        }
        file = new RandomAccessFile(fileName, "rw");
        this.headerCapacity = headerCapacity;

        this.version = header.version;
        file.writeInt(version);
        file.writeInt(headerCapacity);
        this.nSamples = header.nSamples;
        file.writeInt(nSamples);
        this.samplingFrequency = 20000;
        file.writeInt(samplingFrequency);
        this.time = 0;
        file.writeLong(time);

       
        // - shlens additions to neurons file
        this.maxContamination = header.maxContamination;
        file.writeDouble(maxContamination); // 8 bytes
        
        this.minNeuronSpikes = header.minNeuronSpikes;
        file.writeInt(minNeuronSpikes); // 4 bytes
        this.removeDuplicates = header.removeDuplicates;
        file.writeInt(removeDuplicates); // 4 bytes
        
        //mgrivich addition
        this.covarianceType = header.covarianceType;
        file.writeInt(covarianceType);
   
        this.visionVersion = header.visionVersion;
        file.writeInt(visionVersion);
    
        //empty space
        for (int i = 0; i < 104; i++) {
            file.writeByte(0);
        }
      


  //      for (int i = 0; i < 128; i++) {
  //          file.writeByte(0);
  //      }

        electrodeList = new int[headerCapacity];
        Arrays.fill(electrodeList, UNUSED_SLOT);
        neuronIDList = new int[headerCapacity];
        Arrays.fill(neuronIDList, -1);
        startLocationList = new long[headerCapacity];
        flushAll();

        if (version == INT_VERSION && ttl instanceof int[]) {
            addNeuron(0, -1, (int[]) ttl, nTTL);
        } else if (version == DOUBLE_VERSION && ttl instanceof double[]) {
            addExactNeuron(0, -1, (double[]) ttl, nTTL);
        } else {
            throw new Error("Wrong TTL object");
        }

        automaticNeuronID = 1;
    }


    public NeuronFile(String fileName) throws IOException {
        this.fileName = fileName;
        File f = new File(fileName);
        if (!f.exists() || f.isDirectory()) {
            throw new IllegalArgumentException(
                "The neuron file " + fileName + " does not exist or cannot be read");
        }
        file = new RandomAccessFile(fileName, "rw");

        this.version = file.readInt();
        if (DEBUG) System.out.println("version : " + version);
        if (version == INT_VERSION) {
            valueLength = 4;
        } else if (version == DOUBLE_VERSION) {
            valueLength = 8;
        } else if (version == SALK_VERSION) {
            throw new Error(
                "The neuron file " + fileName + " has version " + SALK_VERSION +
                " used by EJ. Please modify the version to " + DOUBLE_VERSION);
        } else {
            throw new Error("Wrong version, ");
        }

        this.headerCapacity = file.readInt();
        if (DEBUG) System.out.println("headerCapacity : " + headerCapacity);
        
        this.nSamples = file.readInt();
        if (DEBUG) System.out.println("nSamples : " + nSamples);
        
        this.samplingFrequency = file.readInt();
        if (DEBUG) System.out.println("samplingFrequency : " + samplingFrequency);
        
        this.time = file.readLong();
        if (DEBUG) System.out.println("time : " + time);
        
        this.maxContamination = file.readDouble();
        if (DEBUG) System.out.println("maxContamination: " + maxContamination);
        
        this.minNeuronSpikes = file.readInt();
        if (DEBUG) System.out.println("minNeuronSpikes: " + minNeuronSpikes);
        
        this.removeDuplicates = file.readInt();
        if (DEBUG) System.out.println("removeDuplicates: " + removeDuplicates);
        
        this.covarianceType = file.readInt();
        if (DEBUG) System.out.println("covarianceType: " + covarianceType); 
        
        this.visionVersion = file.readInt();
        if (DEBUG) System.out.println("visionVersion: " + visionVersion);
        
        for (int i = 0; i < 104; i++) file.readByte();

        this.header = new VisionHeader();
        header.initializeWithGarbage();

        /*****************************/
        // Additions by shlens to fix bugs in NeuronsFile.getHeader()
        
        //NeuronsFile doesn't actually use the vision header.  This is a crude hack to make it look 
        //like it does.  -- mgrivich
        header.version = version;

        /*****************************/
        header.nSamples = nSamples;
        header.samplingFrequency = samplingFrequency;
        header.maxContamination = maxContamination;
        header.minNeuronSpikes = minNeuronSpikes;
        header.removeDuplicates = removeDuplicates;
        header.covarianceType = covarianceType;
        header.visionVersion = visionVersion;

        byte[] buffer = new byte[headerCapacity * SLOT_LENGTH];
        file.readFully(buffer);
        DataInputStream di = new DataInputStream(new ByteArrayInputStream(buffer));
        
        electrodeList = new int[headerCapacity];
        neuronIDList = new int[headerCapacity];
        startLocationList = new long[headerCapacity];
        int maxID = -1;
        for (int i = 0; i < headerCapacity; i++) {
            neuronIDList[i] = di.readInt();
            if (neuronIDList[i] > maxID) {
                maxID = neuronIDList[i];
            }
            electrodeList[i] = di.readInt();
            if (DEBUG) {
           //     System.out.println("ID: " + neuronIDList[i] + " on " + electrodeList[i]);
            }
            startLocationList[i] = di.readLong();
        }

        di.close();

        automaticNeuronID = maxID + 1;

//        int[] idList = getIDList();
//        int nn = 0;
//        for (int i = 0; i < idList.length; i++) {
//            nn += getSpikeCount(idList[i]);
//        }
//        System.out.println("Total spikes: " + nn);
    }

    
    synchronized public int[] getIDList() {
        IntegerList list = new IntegerList();

        for (int i = 1; i < headerCapacity; i++) {
            if (!isEmpty(i)) {
//                System.out.println(i);
                list.add(neuronIDList[i]);
            }
        }

        int[] idList = list.toArray();
        Arrays.sort(idList);
        return idList;
    }


    /**
     * Use only if you really need it !
     *
     * @return int[]
     */
    synchronized public int[] getFullIDList() {
        IntegerList list = new IntegerList();

        for (int i = 1; i < headerCapacity; i++) {
            if (electrodeList[i] != UNUSED_SLOT) {
                list.add(neuronIDList[i]);
            }
        }

        int[] idList = list.toArray();
        Arrays.sort(idList);
        return idList;
    }


    synchronized private void flush(int i) throws IOException {
        file.seek(HEADER_LENGTH + i * SLOT_LENGTH);
        file.writeInt(neuronIDList[i]);
        file.writeInt(electrodeList[i]);
        file.writeLong(startLocationList[i]);
    }


    synchronized private void flushAll() throws IOException {
        file.seek(HEADER_LENGTH);

        for (int i = 0; i < headerCapacity; i++) {
            file.writeInt(neuronIDList[i]);
            file.writeInt(electrodeList[i]);
            file.writeLong(startLocationList[i]);
        }
    }


    synchronized public int getNumberOfNeurons() {
        int n = 0;

        // consider all slots except 0 which is the TTL slot
        for (int i = 1; i < headerCapacity; i++) {
            if (!isEmpty(i)) {
                n++;
            }
        }

        return n;
    }


    synchronized public boolean containsID(int neuronID) {
        for (int i = 0; i < headerCapacity; i++)
            if (neuronIDList[i] == neuronID && !isEmpty(i)) return true;
        return false;
    }


    synchronized private boolean isEmpty(int slotIndex) {
        if (electrodeList[slotIndex] < 0) {
            return true;
        } else {
            return false;
        }
    }


    synchronized public int getHeaderCapacity() {
        return headerCapacity;
    }


    synchronized private int getSomeEmptySlot() {
        for (int i = 0; i < headerCapacity; i++) {
            if (isEmpty(i)) {
                return i;
            }
        }

        return -1;
    }


    /**
     * Adds a new neuron to the file.
     * The ID is automatically assigned (without collisions).
     *
     * @param electrode
     * @param times
     * @param nSpikes
     * @return the associated ID number
     * @throws IOException
     */

    // Removed by shlens. All neurons can only be appended by first determining
    // the neuronID using the static functions above.

    /*    synchronized public int addNeuron(
            int electrode, int[] times, int nSpikes) throws IOException {

            addNeuron(electrode, automaticNeuronID, times, nSpikes);
            automaticNeuronID++;

            return automaticNeuronID - 1;
        }


        synchronized public int addNeuron(
            int electrode, IntegerList times, int nSpikes) throws IOException {

            addNeuron(electrode, automaticNeuronID, times, nSpikes);
            automaticNeuronID++;

            return automaticNeuronID - 1;
        }
     */

    /**
     * Adds a new neuron to the file. The id should be provided and be different that
     * any other ID allready present in the file.
     *
     * @param electrode
     * @param neuronID
     * @param times
     * @param nSpikes
     * @throws IOException
     */
    synchronized public void addNeuron(
        int electrode, int neuronID, int[] times, int nSpikes) throws IOException {

        if (version != INT_VERSION) {
            throw new Error("Wrong version");
        }

        int cellIndex = getIndexForID(neuronID);
        if (cellIndex == -1) {
            cellIndex = getSomeEmptySlot();
            if (cellIndex == -1) {
                throw new IOException("The header is full.");
            }
        } else {
//            System.out.println("Neuron " + neuronID + " allready exists. Override it.");
        }
        neuronIDList[cellIndex] = neuronID;
        electrodeList[cellIndex] = electrode;
        startLocationList[cellIndex] = file.length();

        byte[] writeBuffer = new byte[nSpikes * valueLength];
        for (int i = 0, byteIndex = 0; i < nSpikes; i++) {
            int time = times[i];
            writeBuffer[byteIndex++] = (byte) ( (time >>> 24) & 0xFF);
            writeBuffer[byteIndex++] = (byte) ( (time >>> 16) & 0xFF);
            writeBuffer[byteIndex++] = (byte) ( (time >>> 8) & 0xFF);
            writeBuffer[byteIndex++] = (byte) ( (time >>> 0) & 0xFF);
        }
        file.seek(startLocationList[cellIndex]);
        file.writeInt(nSpikes);
        file.write(writeBuffer, 0, nSpikes * valueLength);

        flush(cellIndex);
    }


    synchronized public int getNextID() {
        return automaticNeuronID;
    }


    synchronized public int addExactNeuron(
        int electrode, double[] times, int nSpikes) throws IOException {

        addExactNeuron(electrode, automaticNeuronID, times, nSpikes);
        automaticNeuronID++;

        return automaticNeuronID - 1;
    }


    synchronized public void addExactNeuron(
        int electrode, int neuronID, double[] times, int nSpikes) throws IOException {

        if (version != DOUBLE_VERSION) {
            throw new Error("Wrong version");
        }

        int cellIndex = getIndexForID(neuronID);
        if (cellIndex == -1) {
            cellIndex = getSomeEmptySlot();
            if (cellIndex == -1) {
                throw new IOException("The header is full.");
            }
        } else {
            System.out.println("Neuron " + neuronID + " allready exists. Override it.");
        }
        neuronIDList[cellIndex] = neuronID;
        electrodeList[cellIndex] = electrode;
        startLocationList[cellIndex] = file.length();

        byte[] writeBuffer = new byte[nSpikes * valueLength];
        for (int i = 0, byteIndex = 0; i < nSpikes; i++) {
            long v = Double.doubleToLongBits(times[i]);
            writeBuffer[byteIndex++] = (byte) ( (v >>> 56) & 0xFF);
            writeBuffer[byteIndex++] = (byte) ( (v >>> 48) & 0xFF);
            writeBuffer[byteIndex++] = (byte) ( (v >>> 40) & 0xFF);
            writeBuffer[byteIndex++] = (byte) ( (v >>> 32) & 0xFF);
            writeBuffer[byteIndex++] = (byte) ( (v >>> 24) & 0xFF);
            writeBuffer[byteIndex++] = (byte) ( (v >>> 16) & 0xFF);
            writeBuffer[byteIndex++] = (byte) ( (v >>> 8) & 0xFF);
            writeBuffer[byteIndex++] = (byte) ( (v >>> 0) & 0xFF);
        }
        file.seek(startLocationList[cellIndex]);
        file.writeInt(nSpikes);
        file.write(writeBuffer, 0, nSpikes * valueLength);

        flush(cellIndex);
    }


    synchronized public void addNeuron(
        int electrode, int neuronID, IntegerList times, int nSpikes) throws IOException {

        if (version != INT_VERSION) {
            throw new Error("Wrong version");
        }

        int cellIndex = getIndexForID(neuronID);
        if (cellIndex == -1) {
            cellIndex = getSomeEmptySlot();
            if (cellIndex == -1) {
                throw new IOException("The header is full.");
            }
        } else {
            System.out.println("Neuron " + neuronID + " allready exists. Override it.");
        }
        neuronIDList[cellIndex] = neuronID;
        electrodeList[cellIndex] = electrode;
        startLocationList[cellIndex] = file.length();

        byte[] writeBuffer = new byte[nSpikes * valueLength];
        for (int i = 0, byteIndex = 0; i < nSpikes; i++) {
            int time = times.get(i);
            writeBuffer[byteIndex++] = (byte) ( (time >>> 24) & 0xFF);
            writeBuffer[byteIndex++] = (byte) ( (time >>> 16) & 0xFF);
            writeBuffer[byteIndex++] = (byte) ( (time >>> 8) & 0xFF);
            writeBuffer[byteIndex++] = (byte) ( (time >>> 0) & 0xFF);
        }
        file.seek(startLocationList[cellIndex]);
        file.writeInt(nSpikes);
        file.write(writeBuffer, 0, nSpikes * valueLength);

        flush(cellIndex);
    }


    synchronized public void deleteNeuron(int id) throws IOException {
        int index = getIndexForID(id);
        if (index != -1 && electrodeList[index] != UNUSED_SLOT) {
            electrodeList[index] = -Math.abs(electrodeList[index]);
            flush(index);
        }
    }


    synchronized public void undeleteNeurons() throws IOException {
        for (int i = 0; i < headerCapacity; i++) {
            if (electrodeList[i] != UNUSED_SLOT && electrodeList[i] < 0) {
                electrodeList[i] = Math.abs(electrodeList[i]);
            }
        }
        flushAll();
    }


    synchronized private int getIDForIndex(int i) {
        return neuronIDList[i];
    }


    synchronized private int getIndexForID(int id) {
        for (int i = 0; i < headerCapacity; i++)
            if (neuronIDList[i] == id) return i;
        return -1;
    }


    synchronized public int getElectrode(int id) {
        int index = getIndexForID(id);

        if (index == -1) {
            return -1;
        } else {
            return electrodeList[index];
        }
    }


    synchronized public int getSpikeCount(int id) throws IOException {
        int index = getIndexForID(id);
        if (index == -1) {
            return -1;
        } else {
            file.seek(startLocationList[index]);
            return file.readInt();
        }
    }


    synchronized public double[] getExactSpikeTimes(int id) throws IOException {
        if (version != DOUBLE_VERSION) {
            throw new Error("Wrong version");
        }

        int index = getIndexForID(id);
//        System.out.println("index " + index);
        if (index == -1) {
            throw new IllegalArgumentException(
                "The neuron file " + fileName + " does not contain id " + id);
        } else {
            file.seek(startLocationList[index]);
            int nSpikes = file.readInt();
            // FIXME what is the nuffer is too small ?
            file.readFully(buffer, 0, nSpikes * valueLength);
            double[] times = new double[nSpikes];
            int ch1, ch2, ch3, ch4;

            int i = 0;
            for (int spike = 0; spike < nSpikes; spike++) {
                ch1 = buffer[i++] & 0xff;
                ch2 = buffer[i++] & 0xff;
                ch3 = buffer[i++] & 0xff;
                ch4 = buffer[i++] & 0xff;
                int int1 = (ch1 << 24) + (ch2 << 16) + (ch3 << 8) + (ch4 << 0);

                ch1 = buffer[i++] & 0xff;
                ch2 = buffer[i++] & 0xff;
                ch3 = buffer[i++] & 0xff;
                ch4 = buffer[i++] & 0xff;
                int int2 = (ch1 << 24) + (ch2 << 16) + (ch3 << 8) + (ch4 << 0);

                times[spike] = Double.longBitsToDouble(
                    ( (long) (int1) << 32) + (int2 & 0xFFFFFFFFL));
            }

            return times;
        }
    }


    synchronized public int[] getSpikeTimes(int id) throws IOException {
        if (version == INT_VERSION) {
            if (!containsID(id)) {
                throw new IOException(
                    "The neuron file " + fileName + " does not contain id " + id);
            }
            return getSpikeTimesForIndex(getIndexForID(id));
        } else if (version == DOUBLE_VERSION) {
            double[] exactTimes = getExactSpikeTimes(id);
            int[] times = new int[exactTimes.length];
            for (int i = 0; i < times.length; i++) {
                times[i] = (int) Math.round(exactTimes[i]);
            }
            return times;
        } else {
            throw new Error("wrong version");
        }
    }


    synchronized private int[] getSpikeTimesForIndex(final int index) throws IOException {
        if (index == -1) {
            return null;
        } else {
            file.seek(startLocationList[index]);
            int nSpikes = file.readInt();
            
            // FIXME what is the nuffer is too small ?
            file.readFully(buffer, 0, nSpikes * valueLength);
            int[] times = new int[nSpikes];
            int ch1, ch2, ch3, ch4, i = 0;

            for (int spike = 0; spike < nSpikes; spike++) {
                ch1 = buffer[i++] & 0xff;
                ch2 = buffer[i++] & 0xff;
                ch3 = buffer[i++] & 0xff;
                ch4 = buffer[i++] & 0xff;
                times[spike] = (ch1 << 24) + (ch2 << 16) + (ch3 << 8) + (ch4 << 0);
            }

            return times;
        }
    }


    /**
     * This should be designed so that it is safe to run it more than once.  file.close() implements
     * Closeable, which is specified to have no effect if run more than once.
     * 
     * @throws IOException
     */
    synchronized public void close() throws IOException {
        file.close();
    }


    synchronized public int getNumberOfSamples() {
        return nSamples;
    }


    synchronized public int getSamplingFrequency() {
        return samplingFrequency;
    }


    synchronized public int[] getTTLTimes() throws IOException {
//        return getSpikeTimesForIndex(0);
        return getSpikeTimes( -1);
    }


    public static class ExtendedSpike {
        public int neuronID;
        public int time;

        public ExtendedSpike() {}
        
        public ExtendedSpike(int neuronID, int time) {
            this.neuronID = neuronID;
            this.time = time;
        }
        
        public String toString() {
            return neuronID + ": " + time;
        }
    }
    
    
    public static class OrderlySkippingSpikeStream {
        private LinkedList<ExtendedSpike> spikes = new LinkedList<ExtendedSpike>();
        
        /**
         * Gets the spikes (skipping according to spikesToUse), presorts them.  This is not 
         * a particularly efficient implementation but it seems to be good enough.  A multithreaded
         * merging sort will do much better if needed.
         * 
         * @param nf
         * @param spikesToUse
         * @throws IOException
         */
        public OrderlySkippingSpikeStream(NeuronFile nf, int spikesToUse) throws IOException {    		
            Vision.getInstance().sendMessage("Loading spikes...");
            for (int id : nf.getIDList()) {
                int[] times = nf.getSpikeTimes(id);
                
                int step = 1;
                if (spikesToUse > 0) step = times.length / spikesToUse;
                if (step < 1) step = 1;
                
                for (int i = 0; i < times.length; i += step)
                    spikes.add(new ExtendedSpike(id, times[i]));
            }

            Vision.getInstance().sendMessage("Sorting spikes...");
            Collections.sort(spikes, new Comparator<ExtendedSpike>() {
                public int compare(ExtendedSpike s1, ExtendedSpike s2) {
                    return s1.time - s2.time;
                }
            });
        }
        
        public boolean hasNext() { return !spikes.isEmpty(); }
        
        /**
         * Pops from the spike list.
         * 
         * Popping is nice as it should free up some memory (?).  Has downside that
         * it is more expensive than just advancing a pointer, stream cannot be rewound.  
         * May be worth revisiting.
         */
        public ExtendedSpike next() { return spikes.remove(); }
        
        public int getNumberOfSpikes() { return spikes.size(); }
    }
    

    public static class SpikeStream1 {
        private int[][] times;
        private int[] index;
        private ExtendedSpike spike = new ExtendedSpike();


        public SpikeStream1(int[][] times) {
            this.times = times;
            this.index = new int[times.length];

            // set to null all time arrays that have zero length
            // this way the next() method will ignore them and not crash
            for (int id = 0; id < this.times.length; id++) {
                if (this.times[id] != null && this.times[id].length == 0) {
                    this.times[id] = null;
                }
            }
        }


        public boolean hasNext() {
            for (int id = 0; id < times.length; id++) {
                if (times[id] != null) {
                    return true;
                }
            }
            return false;
        }


        public ExtendedSpike next() {
            int minTime = Integer.MAX_VALUE;
            int neuronID = -1;
            for (int id = 0; id < times.length; id++) {
                if (times[id] != null) {
                    if (times[id][index[id]] < minTime) {
                        minTime = times[id][index[id]];
                        neuronID = id;
                    }
                }
            }
            if (neuronID == -1) return null;

            spike.neuronID = neuronID;
            spike.time = times[neuronID][index[neuronID]];
            index[neuronID]++;
            if (index[neuronID] == times[neuronID].length) {
                times[neuronID] = null;
            }

            return spike;
        }


        public int getNumberOfSpikes() {
            int n = 0;
            for (int id = 0; id < times.length; id++) {
                if (times[id] != null) {
                    n += times[id].length;
                }
            }
            return n;
        }
    }


    public VisionHeader getHeader() {
        return header;//(VisionHeader) header.clone();
    }


    public SpikeStream1 getSpikeStream() throws IOException {
        int[] idList = getIDList();
        int maxID = MathUtil.max(idList);
        int[][] times = new int[maxID + 1][];
        for (int i = 0; i < idList.length; i++) {
            times[idList[i]] = this.getSpikeTimes(idList[i]);
        }
        return new SpikeStream1(times);
    }
    
    public OrderlySkippingSpikeStream getOSSS(int spikesToUse) throws IOException {
        return new OrderlySkippingSpikeStream(this, spikesToUse);
    }
    
    public String toString() {
        String newline = System.getProperty("line.separator");
        String toReturn = 
            "headerVersion: "     + version           + newline +
            "headerCapacity: "    + headerCapacity    + newline +
            "nSamples: "          + nSamples          + newline + 
            "samplingFrequency: " + samplingFrequency + newline +
            "time: "              + time              + newline + 
            "maxContamination: "  + maxContamination  + newline +
            "minNeuronSpikes:"    + minNeuronSpikes   + newline + 
            "removeDuplicates: "  + removeDuplicates  + newline + 
            "covarianceType: "    + covarianceType    + newline + 
            "visionVersion: "     + visionVersion     + newline + newline +

            "# neurons: " + getIDList().length + newline;
        return toReturn;
    }


    public static void main(String[] args) throws IOException {
        NeuronFile nf = new NeuronFile(
            "F:\\Good Data\\2000-12-14-0\\data051\\data051.neurons-old");

        PrintStream f = new PrintStream("old.txt");

        f.println(nf.getNumberOfNeurons() + " neurons");
        int[] id = nf.getIDList();
        for (int i = 0; i < id.length; i++) {
            f.print("\n");
            f.print(id[i] + ": " + nf.getElectrode(id[i]) + ": ");

            int[] t = nf.getSpikeTimes(id[i]);
            for (int j = 0; j < t.length; j++) {
                f.print(t[j] + ", ");
            }
        }
        f.close();
    }


    @Override
    /**
     * This is a backup safeguard, mostly here because people seemed to forget to close STAFiles and 
     * EI files appropriately from Matlab.  Do not rely on this to close files; you should always run close
     * explicitly.
     */
    protected void finalize() throws Throwable {
        try {
            close();
        } catch (Exception e) {} finally {
            super.finalize();
        }
    }

}
