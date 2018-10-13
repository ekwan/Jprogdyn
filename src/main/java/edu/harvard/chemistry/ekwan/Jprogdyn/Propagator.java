package edu.harvard.chemistry.ekwan.Jprogdyn;

/**
 * This interface represents an object that calculates the next or previous trajectory point.
 */
public interface Propagator {
    /**
     * Calculates the next trajectory point.
     * @param trajectory the trajectory being propagated
     * @param point the last point
     * @param forwards whether to go forwards (true) or backwards (false)
     * @param isNMRpoint whether NMR shieldings should be evaluated for the next point
     * @return the next point
     */
    public TrajectoryPoint propagate(Trajectory trajectory, TrajectoryPoint point, boolean forwards, boolean isNMRpoint);
}
