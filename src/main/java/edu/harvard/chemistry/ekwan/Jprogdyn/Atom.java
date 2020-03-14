package edu.harvard.chemistry.ekwan.Jprogdyn;

import org.apache.commons.math3.geometry.euclidean.threed.*;
import java.util.*;
import java.io.*;

/**
 * This class represents an atom.  It is immutable.
 */
public class Atom implements Immutable, Serializable {
    /** For serialization. */
    public static final long serialVersionUID = 1L;

    /** The atomic symbol. */ 
    public final String symbol;

    /** The atomic mass in amu. */
    public final double mass;
 
    /** The position of the atom. */
    public final Vector3D position;
   
    /**
     * Constructs an atom.
     * @param symbol the atomic symbol
     * @param mass the atomic mass
     * @param x the x position of the atom in angstroms
     * @param y the y position of the atom in angstroms
     * @param z the z position of the atom in angstroms
     */
    public Atom(String symbol, double mass, double x, double y, double z) {
        this(symbol, mass, new Vector3D(x,y,z));
    }

    /**
     * Constructs an atom.
     * @param symbol the atomic symbol
     * @param mass the atomic mass
     * @param position the position of the atom
     */
    public Atom(String symbol, double mass, Vector3D position) {
        if ( symbol == null )
            throw new NullPointerException("null symbol");
        else if ( symbol.length() == 0 )
            throw new IllegalArgumentException("zero length symbol");
        if ( mass < 0.0 )
            throw new IllegalArgumentException("negative mass");
        this.symbol = symbol;
        this.mass = mass;
        if ( position == null )
            throw new NullPointerException("null position");
        this.position = position;
    }

    /**
     * Returns a copy of this atom 
     * @param newPosition the new position for the atom
     * @return a copy of the atom
     */
    public Atom shift(Vector3D newPosition) {
        return new Atom(symbol, mass, newPosition);
    }

    @Override
    public String toString() {
        return String.format("%-5s  %12.7f   %15.10f   %15.10f   %15.10f", symbol, mass, position.getX(), position.getY(), position.getZ());
    }

    /**
     * Prints a MOLDEN geometry string.
     * @return the description
     */
    public String toTrajString() {
        return String.format("%s %.7f %.7f %.7f\n", symbol, position.getX(), position.getY(), position.getZ());
    }

    @Override
    public int hashCode() {
        return Objects.hash(symbol, mass, position);
    }

    @Override
    public boolean equals(Object obj) {
        if ( obj == null ) 
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof Atom) )
            return false;

        Atom a = (Atom)obj;
        if ( symbol.equals(a.symbol) &&
             mass == a.mass &&
             position.equals(a.position) )
            return true;
        return false;
    }
}
