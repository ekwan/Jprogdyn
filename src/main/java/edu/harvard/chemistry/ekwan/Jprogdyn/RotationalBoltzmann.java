package edu.harvard.chemistry.ekwan.Jprogdyn;

import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.special.*;
import com.google.common.collect.*;
import java.util.*;
import java.io.*;
import java.util.concurrent.*;

/**
 * A Boltzmann-like distribution of angular momenta for a given molecule.  
 * This object should be initialized using a list of Atoms, and 
 * will calculate and store the principal axes, inertia tensor, 
 * and parameters for a  distribution of angular frequency vectors.  When the 
 * getRotation method is called, an angular frequency vector 
 * is obtained from random distribution.
 */
public class RotationalBoltzmann implements Immutable {
    /**
     * Multiplication by this constant converts amu * angstrom^2 
     * to electronVolt * femtosecond^2.  Believe it or not, this is 
     * the only unit conversion factor required, as the Boltzmann
     * distribution is calculated in terms of a unitless parameter.
     */
    private final static double MAGIC_CONSTANT = 103.642692;
    
    /** Principal axes. */
    public final Vector3D axis1, axis2, axis3;

    /** Principal moments in amu * angstrom^2 */
    public final double I1, I2, I3;

    /**
     * Constructs the distribution from a Molecule.
     * @param molecule the input molecule
     */
    public RotationalBoltzmann(Molecule molecule) {
        this(molecule.contents);
    }
    
    /**
    * Constructs the distribution from a list of atoms. 
    * @param contents the molecular geometry to build the distribution from
    */
    public RotationalBoltzmann(List<Atom> contents)
    {
        // we want to move the center of mass to the origin for easy calculation
        Vector3D centerOfMass = RotationalBoltzmann.centerOfMass(contents); 
        List<Atom> shiftedContents = new ArrayList<Atom>();
        for ( Atom atom : contents )
            shiftedContents.add(atom.shift(centerOfMass.negate()));

        // we will now calculate moments
        double Ixx = 0;
        double Iyy = 0;
        double Izz = 0;
        double Ixy = 0;
        double Iyz = 0;
        double Ixz = 0;
        for ( Atom atom : shiftedContents )
            {
                if ( atom.mass == 0 )
                    throw new IllegalArgumentException("Atom mass is zero!");

                Ixx += atom.mass * (atom.position.getY()*atom.position.getY() + atom.position.getZ()*atom.position.getZ());
                Iyy += atom.mass * (atom.position.getX()*atom.position.getX() + atom.position.getZ()*atom.position.getZ());
                Izz += atom.mass * (atom.position.getX()*atom.position.getX() + atom.position.getY()*atom.position.getY());
                Ixy -= atom.mass * atom.position.getX() * atom.position.getY();
                Iyz -= atom.mass * atom.position.getY() * atom.position.getZ();
                Ixz -= atom.mass * atom.position.getX() * atom.position.getZ();
            }

        RealMatrix I = MatrixUtils.createRealMatrix(new double[][]{{Ixx, Ixy, Ixz},{Ixy, Iyy, Iyz},{Ixz, Iyz, Izz}});
        EigenDecomposition eigenI = new EigenDecomposition(I);
       
        // System.out.println(eigenI.getD());
       
        // get the principal axes from the V matrix
        this.axis1 = new Vector3D(eigenI.getV().getColumn(0)).normalize();
        this.axis2 = new Vector3D(eigenI.getV().getColumn(1)).normalize();
        this.axis3 = new Vector3D(eigenI.getV().getColumn(2)).normalize();

        // read the eigenvalues from the EigenDecomposition.
        // certain pathological cases will not have three eigenvalues
        // in this cases, we will store a zero value for some eigenvalues, 
        // and remember to deal with it when obtaining an angular frequency.
        // essentially, we will not rotate on axes with zero eigenvalue.
        double[] evalues = eigenI.getRealEigenvalues();
        if ( evalues.length == 0 )
            throw new IllegalArgumentException ("Given atom group has no valid rotations!");
        else if ( evalues.length == 1 )
            {// This case is silly, and will never occur.
                I1 = evalues[0];
                I2 = 0;
                I3 = 0;
            }
        else if ( evalues.length == 2 )
            {// This is the linear case.
                I1 = evalues[0];
                I2 = evalues[1];
                I3 = 0;
            }
        else
            {
                I1 = evalues[0];
                I2 = evalues[1];
                I3 = evalues[2];
            }
    }

    /**
     * Returns a random angular frequency from the distribution.  The user
     * inputs a temperature, and three random numbers.  These numbers correspond
     * to the rotational energy around each of the 3 principal axes.  The sign
     * of the number specifies the direction, while the magnitude of the number
     * represents the fraction of molecules that should have less rotational
     * energy in this mode in our distribution.  The returned angular momentum
     * is in regular Cartesian coordinates -- anyone outside this class
     * should not have to think about the principal axes.
     * @param kT a double representing the temperature in ELECTRONVOLTS
     * @param x1 a random number between -1 and 1 for the first angular frequency component
     * @param x2 a random number between -1 and 1 for the second angular frequency component
     * @param x3 a random number between -1 and 1 for the third angular frequency component
     * @return omega, the angular frequency in radians per FEMTOSECOND
     */
    public Vector3D getOmega(double kT, double x1, double x2, double x3)
    {
        double[] omega = {0,0,0};
        // we remove the signs and put them back later
        int sign1 = (int)Math.signum(x1);
        int sign2 = (int)Math.signum(x2);
        int sign3 = (int)Math.signum(x3);
        x1 = Math.abs(x1);
        x2 = Math.abs(x2);
        x3 = Math.abs(x3);
        double energy1 = 0;
        double energy2 = 0;
        double energy3 = 0;

        if (I1 > 0)
            {
                energy1 = kT * inverseCumulativeBoltzmann(x1);
                omega[0] = sign1 * Math.sqrt((2 * energy1)/(I1 * MAGIC_CONSTANT));
            }
        if (I2 > 0)
            {
                energy2 = kT * inverseCumulativeBoltzmann(x2);
                omega[1] = sign2 * Math.sqrt((2 * energy2)/(I2 * MAGIC_CONSTANT));
            }
        if (I3 > 0)
            {
                energy3 = kT * inverseCumulativeBoltzmann(x3);
                omega[2] = sign3 * Math.sqrt((2 * energy3)/(I1 * MAGIC_CONSTANT));
            }

        double conversion = 3.8267327959301E-23 * Units.AVOGADROS_NUMBER;
        double energy1_kcal = energy1 * conversion;
        double energy2_kcal = energy2 * conversion;
        double energy3_kcal = energy3 * conversion;
        double totalEnergy_kcal = energy1_kcal + energy2_kcal + energy3_kcal;
        //System.out.printf("Total rotational energy: %.4f kcal/mol (%.4f axis1, %.4f axis2, %.4f axis3)\n", totalEnergy_kcal, energy1_kcal, energy2_kcal, energy3_kcal);
        
        // transform these back into regular Cartesian coordinates
        return new Vector3D(omega[0], axis1, omega[1], axis2, omega[2], axis3);
    }

    /**
     * Calculates the center of mass vector for a list of atoms.
     */
    private static Vector3D centerOfMass(List<Atom> contents)
    {
        double totalMass = 0.0;
        Vector3D weightedPosition = Vector3D.ZERO;
        for ( Atom atom : contents )
            {
                totalMass += atom.mass;
                weightedPosition = weightedPosition.add(atom.mass, atom.position);
            }
        return weightedPosition.scalarMultiply(1.0/totalMass); 
    }

    /**
     * Returns a random energy drawn from a Boltzmann distribution for the specified temperature.
     * Will return a mean energy of 0.5 kT.  0.5 kT is about 0.3 kcal/mol at 298 K.
     * @param temperature the temperature in K
     * @return the random energy in kcal/mol
     */
    public static double getRandomBoltzmannEnergy(double temperature)
    {
        double randomNumber = ThreadLocalRandom.current().nextDouble();
        
        //
        //   kcal
        // ------- * K = kcal/mol ---> must use R, not k
        //  mol K
        //
        double kT = Units.R_GAS_KCAL * temperature;
        double energy = inverseCumulativeBoltzmann(randomNumber) * kT;
        return energy;
    }

    /** How many kT to search up to. */
    public static final double CUTOFF_ENERGY = 10.0;

    /**
     * Inverse Cumulative Boltzmann helper method.  Calculates the 
     * rotational kinetic energy for which x fraction of molecules have less
     * kinetic energy.  Normalized so that kT = 1.
     * Unit conversions are to be made outside the method.
     * @param x the fraction of molecules with less rotational energy
     * @return energy divided by kT
     */
    private static double inverseCumulativeBoltzmann(double x)
    {
        if (x > 1 || x < 0)
            throw new IllegalArgumentException("Illegal rotational energy specification.");
        
        // Find the inverse of this monotonically increasing function
        // in two steps.  First scan up to 10 kT with large tolerance,
        // then scan in smaller steps around the solution.  The cutoff at 10 kT
        // is arbitrary, but we should not be using energies that high anyway.
        double trialEnergy = -1;
        for ( double i = 0; i < CUTOFF_ENERGY; i += 0.01 )
            {
                if ( cumulativeBoltzmann(i) > x )
                    {
                        trialEnergy = i - 0.01;
                        break;
                    }
              //  System.out.println( cumulativeBoltzmann(i) );
            }

        // if we failed to find the inverse on this first pass,
        // we are in the high energy tail of the distribution, and 
        // we simply return the max energy
        if ( trialEnergy == -1 )
            return CUTOFF_ENERGY;

        // now we scan in small steps.
        for ( double i = trialEnergy; i < trialEnergy + 0.01; i += 0.0001 )
            {
                if ( cumulativeBoltzmann(i) > x )
                    {
                        trialEnergy = i - 0.00001;
                        break;
                    }
            }
        return trialEnergy - 0.00001;
    }

    /**
     * Cumulative Boltzmann helper method.  Calculates the cumulative
     * Boltzmann distribution for a positive energy x = E/kT.
     * @param x the energy
     * @return the probability of a molecule having kinetic energy less than x
     */
    private static double cumulativeBoltzmann(double x)
    {
        if ( x < 0 )
            throw new IllegalArgumentException("Tried to use Boltzmann cdf with negative energy!");
        return Erf.erf(Math.sqrt(x)); // Kinetic energy is normally distributed.

        //return Math.sqrt(x)*(-2*Math.exp(-x)/Math.sqrt(Math.PI) + Erf.erf(Math.sqrt(x))/Math.sqrt(x));
    }

    /**
     * This is a static utility method that calculates velocities for rotation.
     * Given a group of atoms and an origin, the method should return a matching
     * List of velocities for these atoms that performs a rotation at angular
     * frequency omega.
     * @param contents a list of atoms to be rotated
     * @param omega the angular frequency vector in radians per FEMTOSECOND
     * @param origin the origin of rotation; it is recommended that 
     *        this be set to the origin of the coordinate system
     *        so that center of mass rotations will not be confused
     *        with other rotations.  See alias methods.
     * @return vector of velocities in the same order as the original
    */
    public static List<Vector3D> getVelocities(List<Atom> contents, Vector3D omega, Vector3D origin)
    {
        List<Atom> shiftedContents = new ArrayList<Atom>();
        for ( Atom atom : contents )
            shiftedContents.add(atom.shift(origin.negate()));
        
        List<Vector3D> returnList = new ArrayList<>(contents.size());
        for ( Atom atom : contents )
            returnList.add(Vector3D.crossProduct(omega, atom.position));

        return ImmutableList.copyOf(returnList);
    }

    /**
     * Alias method.  Automatically use center of mass as rotation origin.
     * @param contents the molecular geometry
     * @param omega the angular frequency vector
     * @return vector of velocities in same order as the original
     */
    public static List<Vector3D> getVelocities(List<Atom> contents, Vector3D omega) {
        return getVelocities(contents, omega, RotationalBoltzmann.centerOfMass(contents));
    }

    /**
     * Checks equality by checking the contained atoms only.
     */
    @Override
    public boolean equals(Object obj) {
        if ( obj == null )
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof RotationalBoltzmann) )
            return false;
        
        RotationalBoltzmann r = (RotationalBoltzmann)obj;
        if ( axis1.equals(r.axis1) &&
             axis2.equals(r.axis2) && 
             axis3.equals(r.axis3) &&
             I1 == r.I1 &&
             I2 == r.I2 &&
             I3 == r.I3 )
            return true;
        return false;
    }

    /**
    * String representation of the distribution is just the axes and moments.
    */
    @Override
    public String toString() {
        return(I1 + "\n" + axis1 + "\n"
                + I2 + "\n" + axis2 + "\n"
                + I3 + "\n" + axis3 + "\n");
    }

    /**
    * To agree with this.equals, the hash is just that of the atoms.
    */
    public int hashCode() {
        return Objects.hash(axis1, axis2, axis3, I1, I2, I3);
    }
}
