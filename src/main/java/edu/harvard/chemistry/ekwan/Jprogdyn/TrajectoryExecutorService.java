package edu.harvard.chemistry.ekwan.Jprogdyn;

import java.util.*;
import java.io.*;
import java.util.concurrent.*;

/**
 * This class provides static methods to run groups of Trajectories.
 */
public class TrajectoryExecutorService implements Singleton {
    
    /** Thread pool for doing the work. */
    private static final ThreadPoolExecutor EXECUTOR;

    /** Static initializer. */
    static {
        ArrayBlockingQueue<Runnable> queue = new ArrayBlockingQueue<>(100000, true); // fair queue
        int NUMBER_OF_THREADS = Loader.getInteger("number_of_simultaneous_trajectories");
        EXECUTOR = new ThreadPoolExecutor(NUMBER_OF_THREADS, NUMBER_OF_THREADS, Long.MAX_VALUE, TimeUnit.SECONDS, queue);    
    }

    /** Not instantiable. */
    private TrajectoryExecutorService() {
        throw new IllegalArgumentException("not instantiable");
    }

    /**
     * Runs a group of trajectories in parallel.
     * The calling thread blocks until the trajectories are finished.
     * If there is a problem running any of the trajectories, the original trajectory is returned.
     * @param trajectories the trajectories to run
     * @return the completed trajectories
     */
    public List<Trajectory> runTrajectories(List<Trajectory> trajectories) {
        // submit work
        List<Future<Trajectory>> futures = new LinkedList<>();
        for (Trajectory t : trajectories) {
            Future<Trajectory> f = EXECUTOR.submit(t);
            futures.add(f);
        }

        // wait until jobs are done
        while (true) {
            int done = 0;
            for (Future<Trajectory> f : futures) {
                if ( f.isDone() || f.isCancelled() )
                    done++;
            }
            if ( done == futures.size() )
                break;
            try { Thread.sleep(500); }
            catch ( InterruptedException e ) {}
        }

        // return result
        List<Trajectory> completedTrajectories = new LinkedList<>();
        for (int i=0; i < trajectories.size(); i++) {
            Trajectory originalTrajectory = trajectories.get(i);
            Trajectory finishedTrajectory = null;
            Future<Trajectory> f = futures.get(i);
            try { finishedTrajectory = f.get(); }
            catch (Exception e) { }
            if ( finishedTrajectory != null )
                completedTrajectories.add(finishedTrajectory);
            else
                completedTrajectories.add(originalTrajectory);
        }
        return completedTrajectories;
    }
}
