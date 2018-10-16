package edu.harvard.chemistry.ekwan.Jprogdyn;

import org.apache.commons.math3.geometry.euclidean.threed.*;
import java.util.*;
import java.io.*;
import java.util.concurrent.*;
import com.google.common.collect.*;
import org.apache.commons.io.FileUtils;

/**
 * This class represents a molecular dynamics trajectory.  This object is mutable.
 */
public class Trajectory implements Callable<Trajectory>, Serializable {
    
    // Constants

    /** For serialization. */
    public static final long serialVersionUID = 1L;

    /** How often to save progress.  Note that the first 20 points are always saved point by point. */
    public static final int CHECKPOINT_INTERVAL = Loader.getInteger("checkpoint_interval");

    // Fields

    /** The molecule to run the dynamics on. */
    public final Molecule molecule;

    /** The points that have been calculated so far in the forward direction in chronological order (t=1,2,3...). */
    public final List<TrajectoryPoint> forwardPoints;

    /** The points that have been calculated so far in the reverse direction in reverse chronological order (t=-1,-2,-3). */
    public final List<TrajectoryPoint> backwardPoints;

    /** The first point. */
    public TrajectoryPoint initialPoint;

    /** The time between points in fs. */
    public final double timestep;

    /** How many forward trajectory points to do (the initial point counts as the zeroth point). */
    public final int desiredForwardPoints;

    /** How many backward trajectory points to do. */
    public final int desiredBackwardPoints;

    /**
     * The trajectory will be stopped if any of these geometric conditions are met.
     * This is meant for reaction trajectories, so it is not compatible with doing backward points.
     */
    public final List<InternalCoordinate.Condition> stoppingConditions;

    /** How many times to try initializing. */
    public final int maxInitializationAttempts;

    /** The object used to create the initial positions and velocities. */
    public final Initializer initializer;

    /** What method to use for the dynamics, like b3lyp/6-31g(d). */
    public final CalculationMethod dynamicsMethod;

    /** What method to use for the NMR calculations, like b3lyp/6-31g(d). */
    public final CalculationMethod nmrMethod;

    /** How often to take NMR points.  Use 0 for never, 1 for every point, 2 for every other point, and so forth. */
    public final int NMRinterval;

    /** The filename to checkpoint this job to. */
    public final String checkpointFilename;

    /**
     * Constructs a trajectory.
     * @param molecule the molecule to run the dynamics on
     * @param timestep the time between points in fs
     * @param desiredForwardPoints how many forward trajectory points to do (the initial point counts as the zeroth point)
     * @param desiredBackwardPoints how many backward trajectory points to do
     * @param stoppingConditions when to stop running points (for reaction trajectories only, not compatible with backward points) 
     * @param maxInitializationAttempts how many times to try initializing
     * @param initializer the object used to create the initial positions and velocities
     * @param dynamicsMethod the level of theory to use for propagating the trajectory
     * @param nmrMethod the level of theory to use for computing NMR shieldings (if applicable)
     * @param NMRinterval how often to compute NMR shieldings in number of points (if applicable)
     * @param checkpointFilename the filename to checkpoint this job to 
     */
    public Trajectory(Molecule molecule, double timestep, int desiredForwardPoints, int desiredBackwardPoints,
                      List<InternalCoordinate.Condition> stoppingConditions, int maxInitializationAttempts,
                      Initializer initializer, CalculationMethod dynamicsMethod, CalculationMethod nmrMethod,
                      int NMRinterval, String checkpointFilename) {
        
        // check invariants
        if ( molecule == null )
            throw new NullPointerException("null molecule");
        this.molecule = molecule;
        forwardPoints = new ArrayList<>();
        backwardPoints = new ArrayList<>();
        
        if ( timestep <= 0.0 )
            throw new IllegalArgumentException("timestep must be positive");
        this.timestep = timestep;
        
        // setup the trajectory points
        initialPoint = null; // will be initialized during run()

        if ( desiredForwardPoints < 0 || desiredBackwardPoints < 0 )
            throw new IllegalArgumentException("negative number of requested points");
        else if ( desiredForwardPoints == 0 && desiredBackwardPoints == 0 )
            throw new IllegalArgumentException("must request some traj points");
        this.desiredForwardPoints = desiredForwardPoints;
        this.desiredBackwardPoints = desiredBackwardPoints;

        if ( maxInitializationAttempts < 1 )
            throw new IllegalArgumentException("must give the initializer at least once chance");
        this.maxInitializationAttempts = maxInitializationAttempts;

        // check that the stopping conditions are reasonable
        if ( stoppingConditions == null )
            throw new NullPointerException("use a blank list instead of a null for no stopping conditions");
        if ( desiredBackwardPoints > 0 && stoppingConditions.size() > 0 )
            throw new IllegalArgumentException("cannot use reverse points with stopping conditions");
        for ( InternalCoordinate.Condition condition : stoppingConditions )
            {
                InternalCoordinate coordinate = condition.internalCoordinate;
                for (Integer i : coordinate.indices)
                    if ( i >= molecule.contents.size() )
                        throw new IllegalArgumentException("stopping condition index out of bounds: " + coordinate.toString());
            }
        this.stoppingConditions = stoppingConditions;

        // check that the initializer is appropriate for this trajectory
        if ( initializer == null )
            throw new NullPointerException("must provide an initializer");
        if ( initializer.timestep != timestep )
            throw new IllegalArgumentException("mismatch between initializer and trajectory timestep");
        if ( ! initializer.molecule.equals(molecule) )
            throw new IllegalArgumentException("mismatch between initializer and trajectory molecules");
        this.initializer = initializer;

        // check the level of theory
        if ( dynamicsMethod == null )
            throw new NullPointerException("null dynamics method");
        if ( ! dynamicsMethod.equals(initializer.calculationMethod) )
            throw new IllegalArgumentException("mismatch between initializer and dynamics method");
        if ( dynamicsMethod.calculationType != CalculationMethod.CalculationType.ENERGY_AND_FORCE )
            throw new IllegalArgumentException("unexpected calculation type for dynamics");
        if ( nmrMethod != null && nmrMethod.calculationType != CalculationMethod.CalculationType.NMR )
            throw new IllegalArgumentException("unexpected calculation type for NMR");
        this.dynamicsMethod = dynamicsMethod;
        this.nmrMethod = nmrMethod;

        if ( NMRinterval < 0 )
            throw new IllegalArgumentException("nmr interval can't be negative");
        if ( NMRinterval > 0 && nmrMethod == null )
            throw new NullPointerException("must provide NMR method when requesting NMR points");
        if ( NMRinterval == 0 && nmrMethod != null )
            throw new IllegalArgumentException("provided NMR method but did not request any NMR points");
        this.NMRinterval = NMRinterval;

        if ( checkpointFilename == null || checkpointFilename.trim().length() == 0 )
            throw new NullPointerException("checkpoint filename blank");
        this.checkpointFilename = checkpointFilename;
    }

    /**
     * Runs the trajectory.  Checkpoints are saved frequently in the background.
     */
    @Override
    public Trajectory call() {
        // check if a more recent copy of this trajectory is available
        if ( new File(checkpointFilename).exists() ) {
            Trajectory alternative = loadCheckpoint(checkpointFilename);
            if ( alternative == null && new File(checkpointFilename + ".bak").exists() )
                alternative = loadCheckpoint(checkpointFilename + ".bak");
            if ( alternative != null &&
                 ( alternative.forwardPoints.size() > forwardPoints.size() ||
                   alternative.backwardPoints.size() > backwardPoints.size()  ))
                {
                    System.out.printf("Trajectory auto-loaded from checkpoint %s.\n", checkpointFilename);
                    return alternative.call();
                }
        }

        // check if we are done
        if ( this.isDone() ) {
            System.out.println("This trajectory is already finished.");
            return this;
        }

        // initialize if necessary
        if ( initialPoint == null ) {
            initialPoint = initializer.initialize(maxInitializationAttempts);
            checkpoint();
        }

        // run forward points
        int pointsThisCall = 0;
        int extraPoints = desiredForwardPoints-forwardPoints.size();
        for (int i=0; i < extraPoints; i++) {
            // get the most recent point
            TrajectoryPoint oldPoint = null;
            if ( forwardPoints.size() == 0 )
                oldPoint = initialPoint;
            else
                oldPoint = forwardPoints.get(forwardPoints.size()-1);
        
            // ask for NMR point if requested
            boolean doNMR = false;
            if ( NMRinterval > 0 && forwardPoints.size() % NMRinterval == 0 )
                doNMR = true;

            // get the next point
            TrajectoryPoint newPoint = Propagators.VELOCITY_VERLET.propagate(this, oldPoint, true, doNMR); // true indicates forwards
            pointsThisCall++;

            // update
            forwardPoints.add(newPoint);

            // print update
            double totalEnergy = newPoint.totalEnergy;
            double evalTime = newPoint.evaluationTime;
            System.out.printf("Forward point %d of %d (TE = %.8f hartree, eval time = %.1f)\n",
                              forwardPoints.size(), desiredForwardPoints, totalEnergy, evalTime);

            // checkpoint if requested
            // fixed so that checkpointing is done if the trajectory is finished
            if ( forwardPoints.size() % CHECKPOINT_INTERVAL == 0 || forwardPoints.size() == desiredForwardPoints || 
                 pointsThisCall < 20 || checkTerminationConditions(newPoint, false) )
                checkpoint();
        
            // check geometric termination conditions
            if ( checkTerminationConditions(newPoint, true) )
                return this;
        }

        // run backward points
        extraPoints = desiredBackwardPoints-backwardPoints.size();
        for (int i=0; i < extraPoints; i++) {
            // get the most recent point
            TrajectoryPoint oldPoint = null;
            if ( backwardPoints.size() == 0 )
                oldPoint = initialPoint;
            else
                oldPoint = backwardPoints.get(backwardPoints.size()-1);
        
            // ask for NMR point if requested
            boolean doNMR = false;
            if ( NMRinterval > 0 && backwardPoints.size() % NMRinterval == 0 )
                doNMR = true;

            // get the next point
            TrajectoryPoint newPoint = Propagators.VELOCITY_VERLET.propagate(this, oldPoint, false, doNMR);
            pointsThisCall++;

            // update
            backwardPoints.add(newPoint);

            // print update
            double totalEnergy = newPoint.totalEnergy;
            double evalTime = newPoint.evaluationTime;
            System.out.printf("Backward point %d of %d (TE = %.8f hartree, eval time = %.1f)\n",
                              backwardPoints.size(), desiredBackwardPoints, totalEnergy, evalTime);

            // checkpoint if requested
            if ( backwardPoints.size() % CHECKPOINT_INTERVAL == 0 || backwardPoints.size() == desiredForwardPoints || pointsThisCall < 20 )
                checkpoint();
        }

        System.out.println("All trajectory points are complete.");
        return this;
    }
    
    /**
     * Checks to see if the geometric stopping conditions have been met on the specified point.
     * @param point the point to check
     * @param verbose print more information about what condition was reached
     * @return true if the stopping conditions have been met
     */
    public boolean checkTerminationConditions(TrajectoryPoint point, boolean verbose)
    {
        if ( point == null )
            throw new NullPointerException("cannot check termination conditions on null point");
        List<Vector3D> positions = point.positions;
        for (InternalCoordinate.Condition condition : stoppingConditions) {
            if ( condition.reached(positions) ) {
                if (verbose) {
                    double currentValue = condition.internalCoordinate.getValue(positions);
                    System.out.printf("Termination condition reached: %s (actual value = %.2f)\n", condition.toString(), currentValue);
                }
                return true;
            }
        }
        return false;
    }

    /**
     * Tells us if this trajectory is finished.
     * @return true if we are finished
     */
    public boolean isDone() {
        //System.out.println(initialPoint == null);
        //System.out.printf("forward %d reverse %d\n", forwardPoints.size(), backwardPoints.size());
        if ( initialPoint == null )
            return false;
        if ( forwardPoints.size() < desiredForwardPoints ) {
            if ( forwardPoints.size() > 0 ) {
                    TrajectoryPoint lastPoint = forwardPoints.get(forwardPoints.size()-1);
                    if ( checkTerminationConditions(lastPoint, false) )
                        return true;
                }
            return false;
        }
        if ( backwardPoints.size() < desiredBackwardPoints )
            return false;
        return true;
    }

    /**
     * Write this Trajectory to its checkpoint file.
     */
    private void checkpoint() {
        // if there is already a checkpoint, back it up
        System.out.printf("Serializing...");
        File oldFile = new File(checkpointFilename);
        if ( oldFile.exists() ) {
            File backupFile = new File(checkpointFilename + ".bak");
            try { FileUtils.copyFile(oldFile,backupFile); }
            catch (IOException e) { System.out.printf("Error backing up file: %s\n", checkpointFilename); }
            //System.out.printf("backup made...");
        }

        // serialize the results
        try ( FileOutputStream fileOut = new FileOutputStream(checkpointFilename);
              ObjectOutputStream out = new ObjectOutputStream(fileOut); ) {
            out.writeObject(this);
            System.out.printf("done.  Saved to %s.\n", checkpointFilename);
        }
        catch (Exception e) {
            System.out.println("problem serializing:");
            e.printStackTrace();
        }
    }

    /**
     * Load trajectory from disk.
     * If an error occurs, a message is printed and null is returned.
     * @param filename the checkpoint to load from
     * @return the reconstituted object
     */
    public static Trajectory loadCheckpoint(String filename) {
        try ( FileInputStream fileIn = new FileInputStream(filename);
              ObjectInputStream in = new ObjectInputStream(fileIn); ) {
            Object object = in.readObject();
            if ( ! ( object instanceof Trajectory ) )
                throw new IllegalArgumentException("not a Trajectory");
            Trajectory trajectory = (Trajectory)object;
            return trajectory;
        }
        catch (Exception e) {
            System.out.println("Can't reconstitute Trajectory:");
            e.printStackTrace();
        }
        return null;
    }

    /** 
     * Returns an ordered list of trajectory points. 
     * @return the list of points 
     */
    public List<TrajectoryPoint> getPoints() {
        List<TrajectoryPoint> returnList = new ArrayList<>();

        for ( int i = backwardPoints.size()-1; i >= 0; i-- )
            returnList.add(backwardPoints.get(i));
        if ( initialPoint != null )
            returnList.add( initialPoint );
        for ( int i = 0; i < forwardPoints.size(); i++ )
            returnList.add(forwardPoints.get(i));

        return returnList;
    }

    /**
     * Returns an ordered list of shielding lists.
     * The outer index is the trajectory point index.
     * The inner index is the atom index.
     * @return the list of list of shieldings
     */
    public List<List<Double>> getShieldings() {
        List<List<Double>> returnShieldings = new ArrayList<>();
        List<TrajectoryPoint> orderedList = this.getPoints();
        for ( TrajectoryPoint t : orderedList )
            if ( !(t.shieldings == null) )
                if ( t.shieldings.size() != 0 )
                    returnShieldings.add(t.shieldings);
        return returnShieldings;
    }

    /**
     * Writes a string representation of the positions of each point.
     * This format will allow a movie to be played in MOLDEN.
     * @param filename the filename to write to
     */
    public void writeTrajString(String filename) {
        if ( initialPoint == null )
            throw new NullPointerException("can't write traj string because trajectory not initialized");

        StringBuilder sb = new StringBuilder();
        for (int i=backwardPoints.size()-1; i >= 0; i--) {
            TrajectoryPoint p = backwardPoints.get(i);
            sb.append(p.toTrajString(molecule));
        }
        
        sb.append(initialPoint.toTrajString(molecule));

        for (int i=0; i < forwardPoints.size(); i++) {
            TrajectoryPoint p = forwardPoints.get(i);
            sb.append(p.toTrajString(molecule));
        }

        InputFileFormat.writeStringToDisk(sb.toString(), filename);
    }

    @Override
    public String toString() {
        return String.format("Trajectory (%d forward points, %d backward points, %s)", forwardPoints.size(), backwardPoints.size(), checkpointFilename);
    }

    @Override
    public int hashCode() {
        return Objects.hash(molecule, forwardPoints, backwardPoints, checkpointFilename);
    }

    @Override
    public boolean equals(Object obj) {
        if ( obj == null ) 
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof Trajectory) )
            return false;

        Trajectory t = (Trajectory)obj;
        if ( molecule.equals(t.molecule) &&
             forwardPoints.equals(t.forwardPoints) &&
             backwardPoints.equals(t.backwardPoints) &&
             checkpointFilename.equals(t.checkpointFilename) )
            return true;
        return false;
    }
}
