package edu.harvard.chemistry.ekwan.Jprogdyn;

import java.util.*;
import java.io.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import com.google.common.collect.*;

/**
 * This immutable class represents a molecule.
 */
public class Molecule implements Immutable, Serializable {
    
    /** For serialization. */
    public static final long serialVersionUID = 1L;

    /** Geometry. */
    public final List<Atom> contents;

    /** Normal modes. */
    public final List<NormalMode> modes;

    /** Forces in hartree/bohr. */
    public final List<Vector3D> forces;

    /** Description of this molecule. */
    public final String name;

    /** Overall electrostatic charge. */
    public final int charge;

    /** Overall spin multiplicity. */
    public final int multiplicity;

    /** Potential energy in hartree. */
    public final double potentialEnergy;

    /** NMR absolute shieldings in ppm. */
    public final List<Double> shieldings;

    /**
     * Constructs a molecule.
     * @param contents the atoms that this molecule will contain
     * @param modes the normal modes of the molecule
     * @param forces the Cartesian forces on each atom
     * @param name the description for the molecule
     * @param charge the overall electrostatic charge
     * @param multiplicity the overall spin multiplicity
     * @param potentialEnergy the potential energy in hartree
     * @param shieldings the absolute chemical shieldings in ppm
     */
    public Molecule(List<Atom> contents, List<NormalMode> modes, List<Vector3D> forces,
                    String name, int charge, int multiplicity,
                    double potentialEnergy, List<Double> shieldings) {
        // check invariants
        if ( contents == null || contents.size() == 0 )
            throw new NullPointerException("null or zero length contents");
        this.contents = ImmutableList.copyOf(contents);
        
        if ( modes == null )
            throw new NullPointerException("null modes");
        this.modes = ImmutableList.copyOf(modes);

        if ( forces == null )
            throw new NullPointerException("forces");
        this.forces = ImmutableList.copyOf(forces);

        if ( name == null )
            throw new NullPointerException("null name");
        this.name = name;
        this.charge = charge;
        this.multiplicity = multiplicity;
        this.potentialEnergy = potentialEnergy;
        
        if ( shieldings == null )
            throw new IllegalArgumentException("zero size shieldings");
        this.shieldings = ImmutableList.copyOf(shieldings);
    }

    /**
     * Counts the number of imaginary frequencies.
     * @return the number of imaginary frequencies
     */
    public int numberOfImaginaryFrequencies() {
        int imag = 0;
        for (NormalMode n : modes)
            if ( n.frequency < 0 )
                imag++;
        return imag;
    }

    /**
     * Get the positions of the atoms.
     * @return the atomic positions
     */
    public List<Vector3D> getPositions() {
        List<Vector3D> positions = new ArrayList<>(contents.size());
        for (Atom a : contents)
            positions.add(a.position);
        return positions;
    }

    /**
     * Returns a new molecule with the specified positions.
     * @param newPositions the new positions
     * @return copy of this molecule with the new positions 
     */
    public Molecule setPositions(List<Vector3D> newPositions) {
        int numberOfAtoms = contents.size();
        if ( newPositions.size() != numberOfAtoms )
            throw new IllegalArgumentException("size of new positions array does not match number of atoms");
        List<Atom> newContents = new ArrayList<>(numberOfAtoms);
        for (int i=0; i < numberOfAtoms; i++) {
            Atom oldAtom = contents.get(i);
            Vector3D newPosition = newPositions.get(i);
            Atom newAtom = oldAtom.shift(newPosition);
            newContents.add(newAtom);
        }
        return new Molecule(newContents, modes, forces, name, charge, multiplicity, 0.0, new ArrayList<Double>());
    }


    /**
     * Counts the number of heavy atoms.
     * @return the number of heavy atoms
     */
    public int heavyAtoms() {
        int count = 0;
        for (Atom a : contents)
            if (! a.symbol.equals("H"))
                count++;
        return count;
    }

    /**
     * Writes the geometry in Gaussian format.
     * @return the geometry
     */
    public String geometryString() {
        String returnString = charge + " " + multiplicity + "\n";
        for (int i=0; i < contents.size(); i++)
            returnString += contents.get(i).toString() + "\n";
        return returnString;
    }

    /**
     * Writes the geometry in MOLDEN format.
     * @return the geometry
     */
    public String trajGeometryString() {
        String returnString = "";
        for (int i=0; i < contents.size(); i++)
            returnString += contents.get(i).toTrajString();
        return returnString;
    }

    /**
     * Returns a human-readable description of the normal modes.
     * @return the normal modes
     */
    public String modesString() {
        String returnString = String.format("Normal Modes (%d total):\n\n", modes.size());
        for (int i=0; i < modes.size(); i++) {
                returnString += String.format("Mode %d:\n\n", i);
                returnString += modes.get(i).toString() + "\n-------------\n";
            }
        return returnString;
    }

    @Override
    public String toString() {
        String returnString = String.format("Charge = %d  Multiplicity = %d\n\n", charge, multiplicity);
        returnString       += "Symbol   Weight (amu)          X (A)             Y (A)             Z (A)\n";
        for (Atom a : contents)
            returnString += a.toString() + "\n";
        returnString += String.format("%d normal modes read.", modes.size());
        return returnString;
    }

    @Override
    public int hashCode() {
        return Objects.hash(contents, modes, name, charge, multiplicity, potentialEnergy, shieldings);
    }

    @Override
    public boolean equals(Object obj) {
        if ( obj == null ) 
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof Molecule) )
            return false;

        Molecule m = (Molecule)obj;
        if ( contents.equals(m.contents) &&
             modes.equals(m.modes) &&
             name.equals(m.name) &&
             charge == m.charge && multiplicity == m.multiplicity &&
             potentialEnergy == m.potentialEnergy &&
             shieldings.equals(m.shieldings) )
            return true;
        return false;
    }
}
