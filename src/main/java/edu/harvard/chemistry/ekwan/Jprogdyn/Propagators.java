package edu.harvard.chemistry.ekwan.Jprogdyn;

import java.util.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;

/**
 * This enum contains implementations of Propagator.
 * Only the Velocity Verlet algorithm is currently implemented.
 */
public enum Propagators implements Propagator {
    /**
     * Uses the velocity Verlet algorithm to calculate the next point.  The algorithm is:
     *
     * x(t + Dt) = x(t)  +  v(t) * Dt  +  0.5 * a(t) * Dt^2
     * v(t + Dt) = v(t)  +  0.5 * [ a(t) + a(t+Dt) ] Dt
     *
     * where Dt is the timestep, x(t) is the position at time t, v(t) is the velocity at time t, and
     * a(t) is the acceleration at time t.  The acceleration can be derived from the force divided by the mass.
     * Unlike the regular Verlet, this algorithm is self-starting.
     *
     * Note that because we need the accelerations from both points to calculate the new velocity, we have to
     * make an intermediate TrajectoryPoint that only has the new positions.  This intermediate is not exposed.
     */
    VELOCITY_VERLET {
        @Override
        public TrajectoryPoint propagate(Trajectory trajectory, TrajectoryPoint oldPoint, boolean forwards, boolean isNMRpoint) {
            // check for nulls
            if (trajectory == null)
                throw new NullPointerException("null trajectory");
            if (oldPoint == null)
                throw new NullPointerException("null old point");

            // get required information
            Molecule molecule = trajectory.molecule;
            CalculationMethod dynamicsMethod = trajectory.dynamicsMethod;
            CalculationMethod nmrMethod = trajectory.nmrMethod;

            // set timestep
            double timestep = trajectory.timestep; // in fs
            if ( timestep < 0.0 )
                throw new IllegalArgumentException("timestep must be positive");
            if ( ! forwards )
                timestep = timestep * -1.0; // go backwards

            // initialize lists to hold new verlet quantities
            List<Vector3D> x_t = oldPoint.positions;   // current positions in angstroms
            List<Vector3D> v_t = oldPoint.velocities;  // current velocities in angstroms/femtosecond
            List<Vector3D> a_t = oldPoint.forces2;     // 1/2 * a * Dt^2 in angstroms

            // calculate new positions
            List<Vector3D> x_t_plus = new ArrayList<Vector3D>(x_t);
            for (int i=0; i < x_t_plus.size(); i++) {
                    Vector3D v = x_t_plus.get(i);
                    v = v.add( v_t.get(i).scalarMultiply(timestep) );
                    v = v.add( a_t.get(i) );
                    x_t_plus.set(i, v);
                    Vector3D delta = v.subtract(x_t.get(i));
                }

            // evaluate forces with new geometry
            TrajectoryPoint tempPoint = TrajectoryPoint.create(oldPoint.time + timestep, x_t_plus);
            if ( isNMRpoint )
                tempPoint = tempPoint.evaluatePoint(molecule, nmrMethod, timestep);
            else
                tempPoint = tempPoint.evaluatePoint(molecule, dynamicsMethod, timestep);  

            // get accelerations in A/s^2
            List<Vector3D> oldAccelerations = oldPoint.accelerations;
            List<Vector3D> newAccelerations = tempPoint.accelerations;

            // calculate new velocities
            List<Vector3D> v_t_plus = new ArrayList<Vector3D>(v_t);
            for (int i=0; i < v_t_plus.size(); i++)
                {
                    Vector3D v = v_t_plus.get(i);
                    Vector3D temp = oldAccelerations.get(i).add( newAccelerations.get(i) );
                    temp = temp.scalarMultiply(0.5 * timestep);
                    v = v.add(temp);
                    v_t_plus.set(i, v);
                }

            // calculate kinetic and total energies
            // dimensional analysis:
            //  
            //   g        A^2              fs^2             1E-20 m^2       kg        s^2 J         hartree 
            // ----- * -------- * ---------------------- * ----------- * --------- * -------- * ----------------- = hartree 
            //  mol      fs^2             1E-30 s^2            A^2        1000 g      kg m^2     J_PER_HARTREE J
            //
            // = 1E7 / J_PER_HARTREE * timestep^2
            double kineticEnergy = 0.0;  // in hartree / mol
            double conversionFactor = 1E7 / ( Units.J_PER_HARTREE * timestep * timestep );
            for (int i=0; i < molecule.contents.size(); i++)
                {
                    double mass = molecule.contents.get(i).mass;                   // amu
                    double velocityNorm = v_t_plus.get(i).getNormSq();             // (A/fs)^2
                    kineticEnergy += 0.5 * mass * velocityNorm * conversionFactor; // 0.5 * m * v^2
                }
            double totalEnergy = kineticEnergy + tempPoint.potentialEnergy;  // in hartree
            
            // construct and return new point
            return new TrajectoryPoint(tempPoint.time, kineticEnergy, tempPoint.potentialEnergy, totalEnergy,
                                       tempPoint.positions, v_t_plus, tempPoint.accelerations, tempPoint.forces,
                                       tempPoint.forces2, tempPoint.shieldings, tempPoint.evaluationTime);
        }
    };
}
