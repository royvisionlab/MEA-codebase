package edu.ucsc.neurobiology.vision.test;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ParallelTaskRunner {
    public final int nTasks;
    private TaskRunner[] runners;


    public ParallelTaskRunner(int nTasks) {
        this.nTasks = nTasks;
        runners = new TaskRunner[nTasks];
        for (int i = 0; i < nTasks; i++) {
            runners[i] = new TaskRunner();
        }
    }


    public synchronized void runTasks(Runnable[] tasks) {
//        System.out.println("Running tasks");
        for (int i = 0; i < nTasks; i++) {
            runners[i].runTask(tasks[i]);
        }
    }


    public synchronized void waitUtilDone() {
        for (int i = 0; i < nTasks; i++) {
//            System.out.println("wait " + i);
            runners[i].waitUtilDone();
        }
    }


    public synchronized void die() {
        for (int i = 0; i < nTasks; i++) {
            runners[i].die();
        }
    }


    public static void main(String[] args) {
        final double[] r = new double[1000000];

        Runnable[] tasks = new Runnable[2];
        for (int j = 0; j < tasks.length; j++) {
            tasks[j] = new Runnable() {
                public void run() {
                    for (int i = 0; i < 1000000; i++) {
                        r[i] = Math.sqrt(i);
                    }
                }
            };
        }

        ParallelTaskRunner runner = new ParallelTaskRunner(tasks.length);

        for (int i = 0; i < 1000; i++) {
            System.out.println("" + i);
            runner.runTasks(tasks);
            runner.waitUtilDone();
        }

        runner.die();
    }
}
