package edu.ucsc.neurobiology.vision.test;

/**
 * This is a helper class used in colaboration with the RunnerThread to
 * provide syncronization and communication between the RunnerThread and
 * the thread which executes the action.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class RunnerThreadSync {
    private boolean operationDone = false;


    /**
     * Sets the state of the opperation as being done.
     */
    public synchronized void setOperationDone(boolean operationDone) {
        this.operationDone = operationDone;
    }


    /**
     * Returns the state of the opperation.
     */
    public synchronized boolean isOperationDone() {
        return operationDone;
    }
}
