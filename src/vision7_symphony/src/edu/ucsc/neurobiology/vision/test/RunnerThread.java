package edu.ucsc.neurobiology.vision.test;

import java.awt.event.*;


/**
 * This class represents a Thread that can be run to perform periodically
 * a certain action (assynchronous animation). The working principle of this
 * thread is the following:
 * <ol>
 *   <li> Invoke the assynchronous action (by calling the actionPerformad() method)
 *   <li> Sleep for a specified period of time
 *   <li> When waking up check if the state of the SynchronizationObject reflects
 *     the fact that the action was perfored.
 *     <ol>
 *       <li> if "yes" go to point 1.
 *       <li> if "no" wait() indefinitelly on the SynchronizationObject to be notified
 *         when the assynchronous action gets performed. When the action is really
 *         done go to point 1.
 *     </ol>
 * </ol>
 * <b>Note:</b> The thread that performs the assynchronous action MUST update the state
 *       of the RunnerThreadSync and must call notifyAll() on the RunnerThreadSync
 *       object of this thread after every evaluation of the action.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class RunnerThread
    extends Thread {
    private boolean alive = false;
    private boolean running;
    private int dt;
    private ActionListener l;
    private RunnerThreadSync syncObject;


    /**
     * Creates an RunnerThread with a certain time interval and a specific action
     * to be performed. The thread is NOT started and it needs to be started using the
     * startRunning() method.
     *
     * @param dt The time interval
     * @param l An implementation of the AWT ActionListener inrface.
     */
    public RunnerThread(int dt, ActionListener l) {
        this.dt = dt;
        this.l = l;
        syncObject = new RunnerThreadSync();
        alive = true;
        running = false;
        start();
    }


    /**
     * Used to stop this thread when it is not any more needed.
     * The stopping is clean.
     */
    public void destroyThread() {
        if (alive) {
            alive = false;
            running = false;
        }
    }


    /**
     * Starts the runner to ececute the action.
     */
    public void startRunning() {
        if (!running) {
            running = true;
            synchronized (this) {
                this.notifyAll();
            }
        }
    }


    /**
     * Stops the process without destroying the thread.
     * The process can be continued with the startRunning() method.
     */
    public void stopRunning() {
        if (running) {
            running = false;
            synchronized (this) {
                this.notifyAll();
            }
        }
    }


    /**
     * Returns a boolean wich is true if the thread is running (animating)
     * and false otherwise.
     */
    public boolean isRunning() {
        return running;
    }


    /**
     * Returns the SynchronizationObject of this thread.
     */
    public RunnerThreadSync getSync() {
        return syncObject;
    }


    /**
     * This method is public as an implementation side-effect, <b>DO NOT CALL</b>.
     * It is the run() method of Thread that is overriden here to provide
     * the necessary functionality.
     */
    public void run() {
        while (alive) {
            while (!running) {
                try {
                    synchronized (this) {
                        this.wait();
                    }
                } catch (InterruptedException e) {}
            }

            synchronized (syncObject) {
                syncObject.setOperationDone(false);
            }

            l.actionPerformed(null);

            try {
                sleep(dt);
            } catch (InterruptedException e) {}

            while (!syncObject.isOperationDone()) {
                try {
                    System.out.println("not yet done");
                    synchronized (syncObject) {
                        syncObject.wait();
                    }
                } catch (InterruptedException e) {}
            }
        }
    }

}
