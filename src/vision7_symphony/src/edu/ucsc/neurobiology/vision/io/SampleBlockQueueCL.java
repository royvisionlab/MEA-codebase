package edu.ucsc.neurobiology.vision.io;

import java.io.IOException;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

import com.nativelibs4java.opencl.CLContext;
import com.nativelibs4java.opencl.CLEvent;
import com.nativelibs4java.opencl.CLQueue;

/**
 * Wraps BlockSampleDecompressorCL to provide an asynchronous blocking queue of sample blocks with 
 * efficient reuse of CLBuffers.
 * 
 * Note, currently this is implemented so that if you run out of samples in the file and you're halfway through
 * a buffer then the rest of the buffer will just be left zeros.  This is probably not usually a problem.  Good
 * practice to use buffer lengths that evenly divide the total number of samples though; 1000 is usually a good 
 * choice.
 * 
 * @author Peter H. Li, The Salk Institute
 * @owner Peter H. Li, The Salk Institute
 */
public class SampleBlockQueueCL extends Thread {
    private final BlockSampleDecompressorCL decompressor;
    private final BlockingQueue<AsyncBufferCL<Byte>> inBuffers;
    private final BlockingQueue<AsyncBufferCL<Short>> outBuffers;
    private final BlockingQueue<AsyncBufferCL<Short>> results;
    
    public static SampleBlockQueueCL create(CLContext context, CLQueue queue, String rawDataFileName, 
            int nBuffers, int bufferSize) 
            throws IOException, InterruptedException {
        
        BlockSampleDecompressorCL decompressor = BlockSampleDecompressorCL.create(context, queue, rawDataFileName);
        return new SampleBlockQueueCL(decompressor, nBuffers, nBuffers, bufferSize);
    }
    
    public static SampleBlockQueueCL create(CLContext context, CLQueue queue, String rawDataFileName, 
            int nInBuffers, int nOutBuffers, int bufferSize) 
            throws IOException, InterruptedException {
        
        BlockSampleDecompressorCL decompressor = BlockSampleDecompressorCL.create(context, queue, rawDataFileName);
        return new SampleBlockQueueCL(decompressor, nInBuffers, nOutBuffers, bufferSize);
    }
    
    public SampleBlockQueueCL(BlockSampleDecompressorCL decompressor, 
            int nInBuffers, int nOutBuffers, int bufferSize) 
            throws IOException, InterruptedException {
        super("SampleBlockQueueCL");
        
        this.decompressor = decompressor;
        
        // 1 extra space to accept poisons at the end!
        inBuffers = new ArrayBlockingQueue<AsyncBufferCL<Byte>>(nInBuffers + 1);
        outBuffers = new ArrayBlockingQueue<AsyncBufferCL<Short>>(nOutBuffers + 1);
        results = new ArrayBlockingQueue<AsyncBufferCL<Short>>(nOutBuffers + 1);
        init(nInBuffers, nOutBuffers, bufferSize);
    }
        
    void init(int nInBuffers, int nOutBuffers, int bufferSize) throws InterruptedException {
        for (int i = 0; i < nInBuffers; i++)
            inBuffers.put(new AsyncBufferCL<Byte>(decompressor.createInBuffer(bufferSize)));
        for (int i = 0; i < nOutBuffers; i++)
            outBuffers.put(new AsyncBufferCL<Short>(decompressor.createOutBuffer(bufferSize)));
    }
    
    
    public AsyncBufferCL<Short> take() throws InterruptedException {
        return results.take();
    }

    
    @Override
    public void run() {
        while (true) {
            try {
                AsyncBufferCL<Byte> in = inBuffers.take();
                AsyncBufferCL<Short> out = outBuffers.take();
                if (in.getBuffer() == null || out.getBuffer() == null) break;
                
//				long startTime = System.currentTimeMillis();
                CLEvent decompressEvt = decompressor.decompress(in, out);
//				System.out.println("Finished GPU decompression : " + 
//						out.getBuffer().getElementCount() / decompressor.getHeader().getNumberOfElectrodes() + 
//						" samples in " + (System.currentTimeMillis()-startTime) + " ms");				
                
                in.withEvt(decompressEvt);
                inBuffers.put(in);
                
                if (decompressEvt == null) {
                    close();
                    break;
                }

                results.put(new PooledBufferCL<Short>(out.getBuffer(), outBuffers).withEvt(decompressEvt));
            } catch (InterruptedException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            } catch (IOException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }
    }
    
    
    public void close() throws InterruptedException {
        // Poison the queues
        inBuffers.put(new AsyncBufferCL<Byte>(null));
        outBuffers.put(new AsyncBufferCL<Short>(null));
        results.put(new AsyncBufferCL<Short>(null));

        // TODO: other cleanup
    }

    
    public RawDataHeader getHeader() {
        return decompressor.getHeader();
    }
    
}