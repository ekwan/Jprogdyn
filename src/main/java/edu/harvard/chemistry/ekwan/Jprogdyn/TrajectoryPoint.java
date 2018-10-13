package edu.harvard.chemistry.ekwan.Jprogdyn;

import org.apache.commons.math3.geometry.euclidean.threed.*;
import java.util.*;
import java.io.*;
import com.google.common.collect.*;

/**
 * This class represents a point along a molecular dynamics trajectory.
 */
public class TrajectoryPoint implements Immutable, Serializable {

    /** For serialization. */
    public static final long serialVersionUID = 1L;

    /** The time in fs.  Can be negative for backwards points. */
    public final double time;

    /** The kinetic energy in hartree. */
    public final double kineticEnergy;

    /** The potential energy in hartree. */
    public final double potentialEnergy;

    /** The total energy in kcal. */
    public final double totalEnergy;

    /** The current positions of all atoms in A. */
    public final List<Vector3D> positions;

    /** Cartesian velocities of all the atoms in angstroms/femtosecond. */
    public final List<Vector3D> velocities;

    /** Cartesian accelerations of all atoms in angstroms/fs^2. */
    public final List<Vector3D> accelerations;

    /** Cartesian forces on all atoms in hartrees/bohr. */
    public final List<Vector3D> forces;

    /** Cartesian forces on all atoms in (1/2 acceleration timeStep^2 term) in angstroms. */
    public final List<Vector3D> forces2;

    /** Absolute NMR shieldings in ppm. */
    public final List<Double> shieldings;

    /** Time it took to evaluate this point in seconds. */
    public final double evaluationTime;

    /**
     * Constructs a TrajectoryPoint.
     * @param time the time in fs (negative for backwards points)
     * @param kineticEnergy the kinetic energy in hartree
     * @param potentialEnergy the potential energy in hartree
     * @param totalEnergy the total energy in kcal
     * @param positions the positions of all atoms in Angstroms (cannot be null)
     * @param velocities Cartesian velocities of all the atoms in angstroms/femtosecond
     * @param accelerations accelerations of all atoms in angstroms/fs^2
     * @param forces Cartesian forces on all atoms in hartrees/bohr
     * @param forces2 Cartesian forces on all atoms in (1/2 acceleration timeStep^2 term) in angstroms
     * @param shieldings absolute NMR shieldings in ppm
     * @param evaluationTime time it took to evaluate this point in seconds
    */
    public TrajectoryPoint(double time, double kineticEnergy, double potentialEnergy, double totalEnergy,
                           List<Vector3D> positions, List<Vector3D> velocities, List<Vector3D> accelerations,
                           List<Vector3D> forces, List<Vector3D> forces2, List<Double> shieldings,
                           double evaluationTime) {
        
        // check invariants
        this.time = time;
        
        if ( kineticEnergy < 0.0 )
            throw new IllegalArgumentException("negative kinetic energy");
        this.kineticEnergy = kineticEnergy;
        
        this.potentialEnergy = potentialEnergy;
        this.totalEnergy = totalEnergy;

        if ( positions == null || positions.size() == 0 )
            throw new NullPointerException("null positions");
        this.positions = ImmutableList.copyOf(positions);

        if ( velocities == null )
            this.velocities = null;
        else {
                if ( velocities.size() == 0 )
                    throw new IllegalArgumentException("empty velocities");
                this.velocities = ImmutableList.copyOf(velocities);
            }

        if ( accelerations == null )
            this.accelerations = null;
        else {
                if ( accelerations.size() == 0 )
                    throw new IllegalArgumentException("empty accelerations");
                this.accelerations = ImmutableList.copyOf(accelerations);
            }

        if ( forces == null )
            this.forces = null;
        else {
                if ( forces.size() == 0 )
                    throw new IllegalArgumentException("empty forces");
                this.forces = ImmutableList.copyOf(forces);
            }

        if ( forces2 == null )
            this.forces2 = null;
        else {
                if ( forces2.size() == 0 )
                    throw new IllegalArgumentException("empty forces2");
                this.forces2 = ImmutableList.copyOf(forces2);
            }

        if ( shieldings == null )
            this.shieldings = null;
        else
            this.shieldings = ImmutableList.copyOf(shieldings);

        if (evaluationTime < 0.0)
            throw new IllegalArgumentException("negative evaluation time");
        this.evaluationTime = evaluationTime;
    }

    /**
     * Makes a trajectory point but does not evaluate it.  Created with time and positions only.
     * @param time the time in fs (negative for backwards points)
     * @param positions the current positions of all atoms in Angstroms
     * @return a blank TrajectoryPoint
     */
    public static TrajectoryPoint create(double time, List<Vector3D> positions) {
        return new TrajectoryPoint(time, 0.0, 0.0, 0.0, positions, null, null, null, null, null, 0.0);
    }

    /**
     * Makes a trajectory point but does not evaluate it.  Created with time, positions, and velocities only.
     * @param time the time in fs (negative for backwards points)
     * @param positions the current positions of all atoms in Angstroms
     * @param velocities Cartesian velocities of all the atoms in angstroms/femtosecond
     * @return a blank TrajectoryPoint
     */
    public static TrajectoryPoint create(double time, List<Vector3D> positions, List<Vector3D> velocities) {
        return new TrajectoryPoint(time, 0.0, 0.0, 0.0, positions, velocities, null, null, null, null, 0.0);
    }

    /**
     * Evaluates the point and returns a replacement that has all the fields all filled out.
     * @param molecule needed for header information
     * @param calculationMethod level of theory
     * @param timestep the timestep between the last point and this point in fs
     * @return the evaluated point
     */
    public TrajectoryPoint evaluatePoint(Molecule molecule, CalculationMethod calculationMethod, double timestep) {
        if ( calculationMethod == null )
            throw new NullPointerException("must specify a calculation method!");
        if ( calculationMethod instanceof GaussianCalculationMethod )
            return evaluateGaussianPoint(molecule, (GaussianCalculationMethod)calculationMethod, timestep);
        else
            throw new IllegalArgumentException("unrecognized electronic structure program");
    }

    /**
     * Evaluates a TrajectoryPoint using Gaussian.
     * @param molecule needed for header information
     * @param gaussianCalculationMethod level of theory
     * @param timestep the timestep between the last point and this point in fs
     * @return the evaluated point
    */
    private TrajectoryPoint evaluateGaussianPoint(Molecule molecule, GaussianCalculationMethod gaussianCalculationMethod, double timestep) {
        // create and run job
        GaussianInputFile inputFile = new GaussianInputFile(molecule, positions, gaussianCalculationMethod);
        GaussianJob job = new GaussianJob(inputFile);
        GaussianResult result = job.call();
        GaussianOutputFile outputFile = result.out;
        Molecule m = outputFile.molecule;
        
        // If available, change from forces, which is in hartree/bohr, to forces2, which is in Angstroms.  forces2 is defined as a distance in angstroms:
        //
        // 0.5 * acceleration * timestep^2 term
        //
        // This is useful because we need to calculate x(t) = x(0) + velocity(t) * t + 0.5*a*t^2.
        // force = mass * acceleration, so divide force by mass to get acceleration and multiply by timestep^2.
        //
        // Here is the dimensional analysis.  Capital letters represent constants in the Units class.  Note that a J is kg m^2 s^-2, which cancels
        // out the timestep s^2 and the mass kg.  
        //
        //    hartree     timestep^2 1E-30 s^2     1000 g     KCAL_PER_HARTREE kcal      4184    kg 1E20 angstrom^2           bohr      
        //   --------- * ---------------------- * -------- * ----------------------- * -------- --------------------  * ----------------
        //     bohr           mass * g/mol           kg               hartree            kcal           s^2              BOHR angstroms 
        //
        //  where I replaced J with kg (1E20 angstrom^2) s^-2 in the third last term.  Also note that hartree is implicitly energy/mol.  The numbers are:
        // 
        //    timestep^2 * 1E-30 * 1000 * KCAL_PER_HARTREE * 4184 * 1E20    timestep^2 * 1E-7 * KCAL_PER_HARTREE * 4184
        // = ----------------------------------------------------------- = ---------------------------------------------
        //                                    mass * BOHR                                   mass * BOHR
        //
        // The actual conversion factor needs an extra 0.5.
        List<Vector3D> newForces = m.forces;                             // hartree/bohr
        List<Vector3D> newForces2 = null;                                // in angstroms
        List<Vector3D> newAccelerations = null;                          // in A / fs^2
        if ( newForces != null ) {
            if ( newForces.size() != molecule.contents.size() )
                throw new IllegalArgumentException("list size mismatch: check that forces were calculated");
            
            // calculate 0.5 * a * t^2 term
            newForces2 = new ArrayList<>(newForces.size());
            double conversionFactor = 0.5 * timestep * timestep * 1E-7 * Units.KCAL_PER_HARTREE * 4184 / Units.BOHR;
            for (int i=0; i < molecule.contents.size(); i++) {
                Vector3D rawForce = newForces.get(i);
                double mass = molecule.contents.get(i).mass; // in amu
                Vector3D newForce = rawForce.scalarMultiply(conversionFactor/mass);
                newForces2.add(newForce);
            }
            
            // calculate accelerations in A / fs^2
            //
            // force = mass * acceleration, so divide force by mass to get acceleration
            //
            // dimensional analysis:
            //
            //    hartree     mol          bohr          1000 g     J_PER_HARTREE J     kg m^2   1     1E20 A^2     1E-30 s^2
            //   --------- * ----- * ---------------- * -------- * ----------------- * -------- --- * ---------- * -----------
            //     bohr        g      BOHR angstroms       kg           hartree           s^2    J        m^2         fs^2
            //
            // = 1000 * J_PER_HARTREE * 1E-10 / BOHR
            newAccelerations = new ArrayList<>(newForces.size());
            conversionFactor = 1E-7 * Units.J_PER_HARTREE / Units.BOHR;
            for (int i=0; i < molecule.contents.size(); i++) {
                Vector3D rawForce = newForces.get(i);
                double mass = molecule.contents.get(i).mass; // in amu
                Vector3D acceleration = rawForce.scalarMultiply(conversionFactor/mass);
                newAccelerations.add(acceleration);
            }
        }

        // create new point
        // kinetic and total energy set by Propagator
        return new TrajectoryPoint(time, 0.0, m.potentialEnergy, 0.0,
                                   positions, velocities, newAccelerations, newForces, newForces2,
                                   m.shieldings, result.elapsedTime); 
    }

    /**
     * Writes a block for a MOLDEN movie.
     * @param molecule we need the atom symbols
     * @return the traj string
     */
    public String toTrajString(Molecule molecule) {
        String trajString = positions.size() + "\n" + potentialEnergy + " Jprogdyn t = " + String.format("%.1f\n", time);
        for (int i=0; i < positions.size(); i++) {
                String symbol = molecule.contents.get(i).symbol;
                Vector3D position = positions.get(i);
                trajString = trajString + String.format( "%s %.7f %.7f %.7f\n", symbol, position.getX(), position.getY(), position.getZ() );
            }
        return trajString;
    }

    @Override
    public String toString() {
        return String.format("Trajectory Point: %.1f fs, PE = %.8f, KE = %.8f, TE = %.8f\n", time, potentialEnergy, kineticEnergy, totalEnergy); 
    }

    @Override
    public int hashCode() {
        return Objects.hash(time, kineticEnergy, potentialEnergy, totalEnergy, positions, velocities, accelerations, forces, forces2, shieldings, evaluationTime); 
    }

    @Override
    public boolean equals(Object obj) {
        if ( obj == null ) 
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof TrajectoryPoint) )
            return false;

        TrajectoryPoint p = (TrajectoryPoint)obj;
        if ( time == p.time &&
             kineticEnergy == p.kineticEnergy &&
             potentialEnergy == p.potentialEnergy &&
             totalEnergy == p.totalEnergy &&
             positions.equals(p.positions) &&
             Objects.equals(velocities, p.velocities) &&
             Objects.equals(accelerations, p.accelerations) &&
             Objects.equals(forces, p.forces) &&
             Objects.equals(forces2, p.forces2) &&
             Objects.equals(shieldings, p.shieldings) &&
             evaluationTime == p.evaluationTime )
            return true;
        return false;
    }
}
