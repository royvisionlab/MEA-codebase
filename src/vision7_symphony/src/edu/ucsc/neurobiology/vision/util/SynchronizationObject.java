package edu.ucsc.neurobiology.vision.util;

/**
 * An object used for high level synchronization.
 *
 * 1. Create the object. <br>
 * 2. Call the setWorking() method to put the object in the working state. <br>
 * 3. Call done() to to put the object in the finished state <br>
 *
 * The waitUntilDone() method will return only when the object is in the "finished" state.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class SynchronizationObject {
    public static final int WORKING = 0;
    public static final int DONE = 1;
    private int state = DONE;


    public SynchronizationObject() {
    }


    synchronized public void setWorking() {
        state = WORKING;
    }


    synchronized public void waitUntilDone() {
        if (state == DONE) {
            return;
        } else {
            while (state == WORKING) {
                try {
                    this.wait();
                } catch (InterruptedException e) {}
            }
        }
    }


    synchronized public void done() {
        state = DONE;
        this.notifyAll();
    }


    synchronized public int getState() {
        return state;
    }

}
