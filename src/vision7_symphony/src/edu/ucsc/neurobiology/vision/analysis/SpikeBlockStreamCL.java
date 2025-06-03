package edu.ucsc.neurobiology.vision.analysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;

import edu.ucsc.neurobiology.vision.io.NeuronFile;

public class SpikeBlockStreamCL {
    private final NeuronFile neuronFile;
    private final int nlpoints, nrpoints;
    private final int blockSize, nBlocks, spikesPerBlock;
    private final int[] ids;
    private final ArrayList<LinkedList<Integer>> spikeTimes;
    private int sampleIndex = 0;
    
    public SpikeBlockStreamCL(NeuronFile nf, int spikesToUse, int nlpoints, int nrpoints, int blockSize) throws IOException {
        this.neuronFile = nf;
        this.nlpoints = nlpoints;
        this.nrpoints = nrpoints;
        
        this.blockSize = blockSize;
        this.nBlocks = 1 + nf.getNumberOfSamples() / blockSize;
        this.spikesPerBlock = spikesToUse / nBlocks;
        
        this.ids = neuronFile.getIDList();
        this.spikeTimes = new ArrayList<LinkedList<Integer>>(neuronFile.getNumberOfNeurons());
        for (int n = 0; n < ids.length; n++) {
            int id = ids[n];
            int[] times = neuronFile.getSpikeTimes(id);
            LinkedList<Integer> spikeList = new LinkedList<Integer>();
            for (int t : times) spikeList.add(t);
            spikeTimes.set(n, spikeList);
        }
    }
    
    public List<NeuronSpikeBlock> nextRaw() {
        List<NeuronSpikeBlock> block = new LinkedList<NeuronSpikeBlock>();
        
        for (int n = 0; n < ids.length; n++) {
            LinkedList<Integer> spikeList = spikeTimes.get(n);
            int time;

            // Skip past spikes that are too close to beginning of block
            while((time = spikeList.pop()) < sampleIndex + nlpoints);

            // Collect spikes within block
            LinkedList<Integer> inBlock = new LinkedList<Integer>();
            while((time = spikeList.pop()) < sampleIndex + blockSize - nrpoints)
                inBlock.add(time);
            
            // Skip some spikes to get the right number
            int nInBlock = inBlock.size();
            int step = nInBlock / spikesPerBlock;
            if (step < 1) step = 1;
            
            // Collect selected spikes
            LinkedList<Integer> selected = new LinkedList<Integer>();
            for (int s = 0; s < nInBlock; s += step) selected.add(inBlock.get(s));
            
            block.add(new NeuronSpikeBlock(ids[n], sampleIndex, selected));
        }
        
        Collections.sort(block, new Comparator<NeuronSpikeBlock>() {
            public int compare(NeuronSpikeBlock nsb0, NeuronSpikeBlock nsb1) {
                return nsb1.length - nsb0.length;
            }
        });
        
        sampleIndex += blockSize;
        return block;
    }
    
    
    class NeuronSpikeBlock {
        final int id;
        final int sampleIndex;
        final LinkedList<Integer> spikeList;
        final int length;
        
        public NeuronSpikeBlock(int id, int sampleIndex, LinkedList<Integer> spikeList) {
            this.id = id;
            this.sampleIndex = sampleIndex;
            this.spikeList = spikeList;
            this.length = spikeList.size();
        }
    }
    
    
}