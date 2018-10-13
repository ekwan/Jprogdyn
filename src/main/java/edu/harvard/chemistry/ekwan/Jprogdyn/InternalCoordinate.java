package edu.harvard.chemistry.ekwan.Jprogdyn;

import com.google.common.collect.*;
import java.util.*;
import java.io.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;

/**
 * This abstract class represents some feature of molecular geometry.
 */
public abstract class InternalCoordinate implements Immutable, Serializable {
    
    /** For serialization. */
    public static final long serialVersionUID = 1L;

    /** Atom indices. */
    public final List<Integer> indices;

    /** Description of this coordinate. */
    public final String description;

    /**
     * Generic constructor.
     * @param indices the atom indices
     * @param description description of this coordinate
     */
    public InternalCoordinate(List<Integer> indices, String description)
    {
        if ( indices == null || indices.size() < 2 )
            throw new NullPointerException("must specify at least two atom indices");
        for (Integer i : indices)
            if ( i < 0 )
                throw new IllegalArgumentException("negative index " + i);
        if ( ImmutableSet.copyOf(indices).size() != indices.size() )
            throw new IllegalArgumentException("duplicate atom numbers: " + indices.toString());
        this.indices = indices;
        this.description = description;
    }

    /**
     * Returns the numerical value of this coordinate given some atom positions.
     * @param positions the positions to use
     * @return the value of this internal coordinate
     */
    public abstract double getValue(List<Vector3D> positions);

    @Override
    public String toString()
    {
        if ( description != null && description.trim().length() > 0 )
            return description;
        String returnString = "";
        for (Integer i : indices)
            returnString += String.format("%d-", i);
        return returnString.substring(0, returnString.length()-1);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(indices);
    }

    @Override
    public boolean equals(Object obj)
    {
        if ( obj == null ) 
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof InternalCoordinate) )
            return false;

        InternalCoordinate c = (InternalCoordinate)obj;
        if ( Objects.equals(indices, c.indices) )
            return true;
        return false;
    }

    /** 
     * The condition will return true when the internal coordinate meets this requirement.
     */
    public enum ConditionType {

        /** True if the coordinate is greater than the given condition. */
        GREATER_THAN(">"),

        /** True if the coordinate is equal to the given condition. */
        EQUALS("="),

        /** True if the coordinate is less than or equal to the given condition. */
        LESS_THAN("<=");

        /** Description of this condition type.*/
        public final String description;

        /**
         * Constructor.
         * @param description description of this condition type
         */
        ConditionType(String description) {
            this.description = description;
        }

        @Override
        public String toString() {
            return description;
        }
    }

    /**
     * Represents a stopping condition for a trajectory.
     */
    public static class Condition implements Serializable, Immutable {
        /** For serialization. */
        public static final long serialVersionUID = 1L;

        /** Values within this tolerance of the critical value will be considered equal. */
        public static final double EQUALS_TOLERANCE = 0.20;

        /** The internal coordinate to check. */
        public final InternalCoordinate internalCoordinate;

        /** The type of condition to check. */
        public final ConditionType conditionType;

        /** The critical value. */
        public final double criticalValue;

        /**
         * Create a stopping condition for a trajectory.
         * @param internalCoordinate the internal coordinate to check
         * @param conditionType the type of condition to check
         * @param criticalValue the critical value
         */
        public Condition(InternalCoordinate internalCoordinate, ConditionType conditionType, double criticalValue) {
            if ( internalCoordinate == null )
                throw new NullPointerException("must specify a coordinate");
            this.internalCoordinate = internalCoordinate;

            if ( conditionType == null )
                throw new NullPointerException("must specificy a condition type");
            this.conditionType = conditionType;

            this.criticalValue = criticalValue;
        }

        /**
         * Check if the stopping condition has been reached.
         * @param positions the given atomic positions
         * @return true if the stopping conditions has been reached
         */
        public boolean reached(List<Vector3D> positions) {
            double value = internalCoordinate.getValue(positions);
            if ( conditionType == ConditionType.GREATER_THAN )
                return value > criticalValue;
            else if ( conditionType == ConditionType.LESS_THAN )
                return value < criticalValue;
            else if ( conditionType == ConditionType.EQUALS )
                return Math.abs( value - criticalValue ) < EQUALS_TOLERANCE;
            else
                throw new IllegalArgumentException("unsupported condition type");
        }

        @Override
        public String toString() {
            return String.format("%s %s %.3f", internalCoordinate.toString(), conditionType.toString(), criticalValue);
        }

        @Override
        public int hashCode() {
            return Objects.hash(internalCoordinate, conditionType, criticalValue);
        }

        @Override
        public boolean equals(Object obj) {
            if ( obj == null ) 
                return false;
            if ( obj == this )
                return true;
            if ( !(obj instanceof Condition) )
                return false;

            Condition c = (Condition)obj;
            if ( Objects.equals(internalCoordinate, c.internalCoordinate) &&
                 Objects.equals(conditionType, c.conditionType) &&
                 criticalValue == c.criticalValue )
                return true;
            return false;
        }

    }

    /**
     * Represents a bond length.
     */
    public static class Length extends InternalCoordinate
    {
        /** For serialization. */
        public static final long serialVersionUID = 1L;

        /**
         * Creates a bond length.
         * @param atom1 the index of the first atom (0-indexed)
         * @param atom2 the index of the second atom (0-indexed)
         * @param description the description of this bond length
         */
        public Length(int atom1, int atom2, String description) {
            super(ImmutableList.of(atom1, atom2), description);
        }

        /**
         * Compute the bond length.
         * @param positions the atomic coordinates
         * @return the bond length in Angstroms
         */
        public double getValue(List<Vector3D> positions) {
            Vector3D v1 = positions.get(indices.get(0));
            Vector3D v2 = positions.get(indices.get(1));
            return Vector3D.distance(v1, v2);
        }
    }

    /**
     * Represents a bond angle.
     */
    public static class Angle extends InternalCoordinate {

        /** For serialization. */
        public static final long serialVersionUID = 1L;

        /**
         * Creates a bond angle.
         * @param atom1 the index of the first atom (0-indexed)
         * @param atom2 the index of the second atom (0-indexed)
         * @param atom3 the index of the third atom (0-indexed)
         * @param description the description of this bond angle
          */
        public Angle(int atom1, int atom2, int atom3, String description) {
            super(ImmutableList.of(atom1, atom2, atom3), description);
        }

        /**
         * Compute the bond angle.
         * @param positions the atomic coordinates
         * @return the bond angle in degrees        
         */
        public double getValue(List<Vector3D> positions) {
            Vector3D v1 = positions.get(indices.get(0));
            Vector3D v2 = positions.get(indices.get(1));
            Vector3D v3 = positions.get(indices.get(2));
            Vector3D v1prime = v1.subtract(v2);
            Vector3D v3prime = v3.subtract(v2);
            return Math.toDegrees(Vector3D.angle(v1prime, v3prime));
        }
    }

    /**
     * Represents a bond torsion.
     */
    public static class Torsion extends InternalCoordinate {
        
        /** For serialization. */
        public static final long serialVersionUID = 1L;

        /**
         * Creates a bond torsion.
         * @param atom1 the index of the first atom (0-indexed)
         * @param atom2 the index of the second atom (0-indexed)
         * @param atom3 the index of the third atom (0-indexed)
         * @param atom4 the index of the fourth atom (0-indexed)
         * @param description the description of this bond torsion
        */
        public Torsion(int atom1, int atom2, int atom3, int atom4, String description) {
            super(ImmutableList.of(atom1, atom2, atom3, atom4), description);
        }

        /**
         * Compute the bond torsion.
         * @param positions the atomic coordinates
         * @return the dihedral angle in degrees        
         */
        public double getValue(List<Vector3D> positions) {
            Vector3D v1 = positions.get(indices.get(0));
            Vector3D v2 = positions.get(indices.get(1));
            Vector3D v3 = positions.get(indices.get(2));
            Vector3D v4 = positions.get(indices.get(3));
        
            Vector3D b1 = v2.add(-1.0, v1);
            Vector3D b2 = v3.add(-1.0, v2);
            Vector3D b3 = v4.add(-1.0, v3);

            // make sure the vectors are not collinear
            if (Vector3D.angle(b1,b2) == 0 || Vector3D.angle(b2,b3) == 0) {
                System.out.println("Warning! Collinear dihedral angle!");
                return 0.0;
            }

            // compute dihedral angle
            double term1 = Vector3D.dotProduct(b1.scalarMultiply( b2.getNorm() ), Vector3D.crossProduct(b2, b3) );
            double term2 = Vector3D.dotProduct(Vector3D.crossProduct(b1, b2), Vector3D.crossProduct(b2, b3) );
            return Math.toDegrees( Math.atan2(term1, term2) );
        }
    }
}
