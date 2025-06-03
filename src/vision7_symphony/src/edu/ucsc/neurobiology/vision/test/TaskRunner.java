package edu.ucsc.neurobiology.vision.test;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class TaskRunner
    extends Thread {

    private Runnable task;
    private boolean die;


    public TaskRunner() {
        start();
    }


    public synchronized void runTask(Runnable task) {
//        System.out.println("runTask: task scheduled");

        if (this.task != null) {
            throw new IllegalThreadStateException(
                "TaskRunner is allready executing a task.");
        } else {
            this.task = task;
            notifyAll();
        }
    }


    public synchronized void die() {
        die = true;
        notifyAll();
    }


    public synchronized void waitUtilDone() {
//        System.out.println("waitUtilDone: enter");

        if (task == null) {
            // the task is allready finished or not yet started
            return;
        }

        while (task != null) {
            try {
                wait();
            } catch (InterruptedException e) {}
        }

//        System.out.println("waitUtilDone: exit");
    }


    public synchronized void run() {
        while (!die) {
//            System.out.println("run: waiting for task");
            while (task == null) {
                try {
                    wait();
                } catch (InterruptedException e) {}
                if (die) {
                    return;
                }
            }

            task.run();
            task = null;
//            System.out.println("run: task executed");
            notifyAll();
        }
    }


    public static void main(String[] args) {
        TaskRunner runner = new TaskRunner();
        Runnable task = new Runnable() {
            public void run() {
                for (int i = 0; i < 1000000; i++) {
                    Math.sqrt(i);
                }
            }
        };

        try {
            Thread.sleep(1000);
        } catch (InterruptedException e) {}

        for (int i = 0; i < 100; i++) {
            System.out.println(i);
            runner.runTask(task);
            runner.waitUtilDone();
        }

        runner.die();
    }
}
