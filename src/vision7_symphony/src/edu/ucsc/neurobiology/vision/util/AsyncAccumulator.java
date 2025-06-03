package edu.ucsc.neurobiology.vision.util;

import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

public class AsyncAccumulator<A> implements Runnable {
    private final Accumulator<A> accumulator;
    private final Handler handler;
    private final BlockingQueue<A> queue;
    
    public AsyncAccumulator(Accumulator<A> accumulator, Handler handler, BlockingQueue<A> queue) {
        this.accumulator = accumulator;
        this.handler = handler;
        this.queue = queue;
    }
    
    public void run() {
        try {
            while (!Thread.currentThread().isInterrupted()) queue.put(accumulator.accumulate());			
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        } catch (Throwable t) {
            handler.handle(t);
        }
    }

    public A take() throws InterruptedException { return queue.take(); }
    
    
    static public <A> AsyncAccumulator<A> build(Accumulator<A> accumulator, Handler handler, int queueLength) {
        BlockingQueue<A> queue = new ArrayBlockingQueue<A>(queueLength);
        return new AsyncAccumulator<A>(accumulator, handler, queue);
    }
    
    static public <A> AsyncAccumulator<A> build(Accumulator<A> accumulator, int queueLength) {
        return build(accumulator, new Handler(){
            public void handle(Throwable t) { 
                t.printStackTrace();
                throw new Error(t.getMessage());
            }
        }, queueLength);
    }
    
    
    public interface Accumulator<A> {
        A accumulate() throws Throwable;
    }
    
    public interface Handler {
        void handle(Throwable t);
    }
}