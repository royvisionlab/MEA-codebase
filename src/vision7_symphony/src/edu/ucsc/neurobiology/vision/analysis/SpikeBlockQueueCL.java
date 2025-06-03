package edu.ucsc.neurobiology.vision.analysis;

import java.io.IOException;
import java.util.Collections;
import java.util.LinkedList;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

import org.bridj.Pointer;

import cern.colt.Arrays;

import com.nativelibs4java.opencl.CLBuffer;
import com.nativelibs4java.opencl.CLContext;
import com.nativelibs4java.opencl.CLEvent;
import com.nativelibs4java.opencl.CLMem.Usage;
import com.nativelibs4java.opencl.CLQueue;

import edu.ucsc.neurobiology.vision.analysis.SpikeBlockBuilder.RawSpikeBlock;
import edu.ucsc.neurobiology.vision.analysis.SpikeBlockQueueCL.SpikeBlockCL;
import edu.ucsc.neurobiology.vision.io.AsyncBufferCL;
import edu.ucsc.neurobiology.vision.io.NeuronFile;
import edu.ucsc.neurobiology.vision.io.PooledBufferCL;
import edu.ucsc.neurobiology.vision.util.VisionParams;

public class SpikeBlockQueueCL extends Thread {
    private static final float MAX_SPIKES_PER_CELL_PER_SECOND = 200f;
    private static final float MAX_SPIKES_PER_SAMPLE = MAX_SPIKES_PER_CELL_PER_SECOND / (1000 * (float) VisionParams.SAMPLES_PER_MILLISECOND);
    
    private final CLContext context;
    private final CLQueue clQueue;
    private final SpikeBlockBuilder builder;
    private final BlockingQueue<SpikeBlockCL> emptyBlocks;
    private final BlockingQueue<SpikeBlockCL> fullBlocks;
    private final int maxSpikesPerCellPerBlock;

    public static SpikeBlockQueueCL create(CLContext context, CLQueue clQueue, int nBlockBuffers, 
            NeuronFile nf, int spikesToUse, int nlpoints, int nrpoints, int blockSize) throws IOException, InterruptedException {
        SpikeBlockBuilder sbb = new SpikeBlockBuilder(nf, spikesToUse, nlpoints, nrpoints, blockSize);
        return new SpikeBlockQueueCL(context, clQueue, sbb, nBlockBuffers);
    }
    
    public SpikeBlockQueueCL(CLContext context, CLQueue clQueue, SpikeBlockBuilder builder, int nBlockBuffers) throws InterruptedException {
        super("SpikeBlockQueueCL");
        this.context = context;
        this.clQueue = clQueue;
        this.builder = builder;
        this.emptyBlocks = new ArrayBlockingQueue<SpikeBlockCL>(nBlockBuffers);
        this.fullBlocks  = new ArrayBlockingQueue<SpikeBlockCL>(nBlockBuffers);
        
        this.maxSpikesPerCellPerBlock = Math.round(MAX_SPIKES_PER_SAMPLE*builder.getBlockSize());
        init(nBlockBuffers);
    }
    
    private void init(int nBlockBuffers) throws InterruptedException {
        int nNeurons = builder.getNumNeurons();
        for (int i = 0; i < nBlockBuffers; i++) {
            AsyncBufferCL<Short> cellNums     = new AsyncBufferCL<Short>(context.createBuffer(Usage.Input, Short.class, nNeurons));
            AsyncBufferCL<Short> spikeNums    = new AsyncBufferCL<Short>(context.createBuffer(Usage.Input, Short.class, nNeurons));
            AsyncBufferCL<Short> spikeIndices = new AsyncBufferCL<Short>(context.createBuffer(Usage.Input, Short.class, nNeurons));
            
            // Guess max number of spikes; if this is too small, the whole thing will crash
            CLBuffer<Short> spikeBuffer = context.createBuffer(Usage.Input, Short.class, nNeurons*maxSpikesPerCellPerBlock);
            AsyncBufferCL<Short> spikeList = new AsyncBufferCL<Short>(spikeBuffer);
            
            emptyBlocks.put(new SpikeBlockCL(cellNums, spikeNums, spikeIndices, spikeList, emptyBlocks));
        }
    }
    
    public int getMaxSpikesPerCellPerBlock() { return maxSpikesPerCellPerBlock; }
    
    @Override
    public void run() {
        while (true) {
            RawSpikeBlock rsb = builder.nextRaw();
            try {
                SpikeBlockCL block = emptyBlocks.take();
                block.ensureFill(rsb, clQueue, maxSpikesPerCellPerBlock);
                fullBlocks.put(block);
            } catch (InterruptedException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }
    }

    public SpikeBlockCL take() throws InterruptedException {
        return fullBlocks.take();
    }

    
    /**
     * Translates RawSpikeBlock objects into simple GPU memory buffers
     */
    class SpikeBlockCL {
        public final AsyncBufferCL<Short> cellNums;
        public final AsyncBufferCL<Short> spikeNums;
        public final AsyncBufferCL<Short> spikeIndices;
        public final AsyncBufferCL<Short> spikeList;
        private final BlockingQueue<SpikeBlockCL> owner;
        
        public SpikeBlockCL(AsyncBufferCL<Short> cellNums, AsyncBufferCL<Short> spikeNums, AsyncBufferCL<Short> spikeIndices, AsyncBufferCL<Short> spikeList) {
            this(cellNums, spikeNums, spikeIndices, spikeList, null);
        }

        public SpikeBlockCL(AsyncBufferCL<Short> cellNums, AsyncBufferCL<Short> spikeNums, AsyncBufferCL<Short> spikeIndices, AsyncBufferCL<Short> spikeList,
                BlockingQueue<SpikeBlockCL> owner) {
            this.cellNums     = cellNums;
            this.spikeNums    = spikeNums;
            this.spikeIndices = spikeIndices;
            this.spikeList    = spikeList;
            this.owner = owner;
        }
        
        /**
         * Fill, but throw error if there are more than maxSpikes spikes for any cell
         * @param rsb
         * @param queue
         * @param maxSpikes
         */
        public void ensureFill(RawSpikeBlock rsb, CLQueue queue, int maxSpikes) {
            short[] spikeNums = rsb.getSpikeNums();
            for (short n : spikeNums)
                if (n > maxSpikes) throw new Error("Too many spikes!");
            fill(rsb.getCellNums(), spikeNums, rsb.getSpikeIndices(), rsb.getSpikeList(), queue);
        }
        
        public void fill(RawSpikeBlock rsb, CLQueue queue) {
            fill(rsb.getCellNums(), rsb.getSpikeNums(), rsb.getSpikeIndices(), rsb.getSpikeList(), queue);
        }
        
        /**
         * Ideally this should use non-blocking writes but I'm having problems with memory consistency and
         * I believe blocking writes fixes it.
         * @param cellNums
         * @param spikeNums
         * @param spikeIndices
         * @param spikeList
         * @param queue
         */
        public void fill(short[] cellNums, short[] spikeNums, short[] spikeIndices, short[] spikeList, CLQueue queue) {
            this.cellNums.write(queue, Pointer.pointerToShorts(cellNums), true);
            this.spikeNums.write(queue, Pointer.pointerToShorts(spikeNums), true);
            this.spikeIndices.write(queue, Pointer.pointerToShorts(spikeIndices), true);

            // Make sure buffer is big enough for spikeList, don't try to write empty spike list as pointerToShorts doesn't like this.
            this.spikeList.errEnsure(spikeList.length);
            if (spikeList.length > 0) this.spikeList.write(queue, Pointer.pointerToShorts(spikeList), true);
        }
        
        public void giveBack() {
            if (owner.contains(this)) {
                System.err.println("SpikeBlockCL has already been returned to the queue");
                return;
            }
            
            try {
                owner.put(this);
            } catch (InterruptedException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }
        
        public CLEvent[] getEvts() {
            LinkedList<CLEvent> evts = new LinkedList<CLEvent>();
            for (CLEvent evt : cellNums.getEvts())     evts.add(evt);
            for (CLEvent evt : spikeNums.getEvts())    evts.add(evt);
            for (CLEvent evt : spikeIndices.getEvts()) evts.add(evt);
            for (CLEvent evt : spikeList.getEvts())    evts.add(evt);
            return evts.toArray(new CLEvent[0]);
        }
        
        
        public void withEvt(CLEvent evt) { withEvts(new CLEvent[]{ evt }); }
        public void withEvts(CLEvent[] evts) {
            cellNums.withEvts(evts);
            spikeNums.withEvts(evts);
            spikeIndices.withEvts(evts);
            spikeList.withEvts(evts);
        }

        public void waitFor() {
            cellNums.waitFor();
            spikeNums.waitFor();
            spikeIndices.waitFor();
            spikeList.waitFor();
        }
        
    }
    
}
