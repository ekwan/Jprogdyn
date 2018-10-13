package edu.harvard.chemistry.ekwan.Jprogdyn;

import org.apache.commons.math3.geometry.euclidean.threed.*;
import com.google.common.collect.*;
import java.util.*;
import java.io.*;

/**
 * This immutable class represents a molecular vibration.
 */
public class NormalMode implements Immutable, Serializable {

    /** For serialization. */
    public static final long serialVersionUID = 1L;

    /** Harmonic frequency of this mode in inverse cm. */
    public final double frequency;

    /** Reduced mass in amu. */
    public final double reducedMass;

    /** Force constant in mDyne/A. */
    public final double forceConstant;

    /** The normal displacements in A for each atom.  Indexed by atom.  See {@link Molecule#contents}. */
    public final List<Vector3D> coordinates;

    /**
     * Constructs a NormalMode.
     * @param frequency the frequency in inverse cm
     * @param reducedMass the reduced mass in amu
     * @param forceConstant the force constant in mDyne/A
     * @param coordinates the normal displacements in Angstrom for each atom
     */
    public NormalMode(double frequency, double reducedMass, double forceConstant, List<Vector3D> coordinates) {
        this.frequency = frequency;
        this.reducedMass = reducedMass;
        this.forceConstant = forceConstant;
        
        if ( coordinates == null || coordinates.size() == 0 )
            throw new NullPointerException("no normal coordinates");
        this.coordinates = ImmutableList.copyOf(coordinates);
    }

    @Override
    public String toString() {
        StringBuilder s = new StringBuilder(String.format("\nFrequency (cm-1): %.2f\n", frequency));
                                   s.append(String.format("Reduced Mass (amu): %.2f\n", reducedMass));
                                   s.append(String.format("Force Constant (mDyne/A): %.2f\n",forceConstant));
        s.append("Normal Coordinates (xyz; A):\nAtom Index: X displacement:     Y displacement    Z displacement:\n");
        for (int i=0; i < coordinates.size(); i++)
            {
                Vector3D v = coordinates.get(i);
                s.append(String.format("%5d        %10.5f          %10.5f         %10.5f\n", (i+1), v.getX(), v.getY(), v.getZ()));
            }
        s.append("\n");
        return s.toString();
    }

    @Override
    public int hashCode() {
        return Objects.hash(frequency, reducedMass, forceConstant, coordinates);
    }

    @Override
    public boolean equals(Object obj) {
        if ( obj == null ) 
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof NormalMode) )
            return false;

        NormalMode n = (NormalMode)obj;
        if ( frequency == n.frequency &&
             reducedMass == n.reducedMass &&
             forceConstant == n.forceConstant &&
             coordinates.equals(coordinates) )
            return true;
        return false;
    }
}
