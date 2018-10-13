package edu.harvard.chemistry.ekwan.Jprogdyn;

import java.io.*;

/**
 * This is an abstract class that represents the information needed to perform an
 * electronic structure calculation.  It contains the program, the level of theory,
 * and how much memory/how many processors to use.
 */
public abstract class CalculationMethod implements Immutable, Serializable {
    
    /** For serialization. */
    public static final long serialVersionUID = 1L;

    /** Represents the kind of calculation that will be performed. */
    public enum CalculationType
    {
        /** For trajectory points. */
        ENERGY_AND_FORCE,

        /** For NMR calculations. */
        NMR;
    }

    /** Represents an electronic structure program. */
    public enum Program
    {
        /** Use Gaussian. */
        GAUSSIAN;
    }

    /** Represents the kind of calculation that will be performed. */
    public final CalculationType calculationType;

    /** Which program to use. */
    public final Program program;

    /** Memory to use in GB. */
    public final int memory;

    /** Number of rocessors to use. */
    public final int processors;

    /**
     * Constructs a CalculationMethod.
     * @param calculationType the kind of calculation to run
     * @param program the electornic structure program to use
     * @param memory the amount of memory to use in GB
     * @param processors the number of processors to use
     */
    public CalculationMethod(CalculationType calculationType, Program program,
                             int memory, int processors)
    {
        this.calculationType = calculationType;
        this.program = program;
        if ( memory < 1 )
            throw new IllegalArgumentException("must use at least 1 GB");
        this.memory = memory;
        if ( processors < 1 )
            throw new IllegalArgumentException("must use at least one processor");
        this.processors = processors;
    }
}
