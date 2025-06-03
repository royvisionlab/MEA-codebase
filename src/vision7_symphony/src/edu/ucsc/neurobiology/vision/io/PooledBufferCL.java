package edu.ucsc.neurobiology.vision.io;

import java.util.concurrent.BlockingQueue;
import com.nativelibs4java.opencl.CLBuffer;

/**
 * This CLBuffer is used in pooled fashion, so rather than being released when it is closed, it simply 
 * returns itself to the owning BlockingQueue(s).
 * 
 * @author Peter H. Li, The Salk Institute
 * @owner Peter H. Li, The Salk Institute
 */
public class PooledBufferCL<T> extends AsyncBufferCL<T> {
    private final BlockingQueue<AsyncBufferCL<T>> owner;

    public PooledBufferCL(CLBuffer<T> buffer, BlockingQueue<AsyncBufferCL<T>> owner) {
        super(buffer);
        this.owner = owner;
    }
    
    @Override
    public void close() {		
        if (owner.contains(this)) {
            System.err.println("PooledBufferCL has already been returned to the queue");
            return;
        }
        
        try {
            owner.put(this);
        } catch (InterruptedException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }
}
