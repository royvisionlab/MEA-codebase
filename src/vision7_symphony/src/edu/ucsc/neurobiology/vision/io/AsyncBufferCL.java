package edu.ucsc.neurobiology.vision.io;

import org.bridj.Pointer;

import com.nativelibs4java.opencl.CLBuffer;
import com.nativelibs4java.opencl.CLEvent;
import com.nativelibs4java.opencl.CLMem.Usage;
import com.nativelibs4java.opencl.CLQueue;

/**
 * Wraps CLBuffer along with the CLEvents that indicate when the buffer is ready to be used.
 * 
 * @author Peter H. Li, The Salk Institute
 * @owner Peter H. Li, The Salk Institute
 */
public class AsyncBufferCL<T> {
    private static CLEvent[] EMPTY_EVTS = new CLEvent[0];
    
    private CLBuffer<T> buffer;
    private CLEvent[] evts = EMPTY_EVTS;
    private boolean released = false;
    
    public AsyncBufferCL(CLBuffer<T> buffer) { 
        this.buffer = buffer;
    }
    
    public CLBuffer<T> getBuffer() { return buffer; }
    public CLEvent[] getEvts() { return evts; }

    public AsyncBufferCL<T> withEvts(CLEvent[] evts) { 
        this.evts = evts; 
        return this;
    }
    public AsyncBufferCL<T> withEvt(CLEvent evt) { return withEvts(new CLEvent[]{ evt }); }
    
    public void waitFor() {
        for (CLEvent e : evts) e.waitFor();
        evts = EMPTY_EVTS;
    }
    
    public void close() {
        if (released) return;
        buffer.release();
        released = true;
    }
        
    public void write(CLQueue queue, Pointer<T> pnt, boolean blocking) {
        if (blocking) {
            buffer.write(queue, pnt, blocking, evts);
            evts = EMPTY_EVTS;
        } else {
            CLEvent writeEvt = buffer.write(queue, pnt, blocking, evts);
            withEvt(writeEvt);
        }
    }
    
    public void clearEnsure(long size) {
        if (buffer.getElementCount() >= size) return;
        buffer = buffer.getContext().createBuffer(Usage.InputOutput, buffer.getElementClass(), size);
    }

    public void errEnsure(long size) {
        if (buffer.getElementCount() < size) throw new Error("Buffer wasn't big enough!");
    }
    
}