package edu.harvard.chemistry.ekwan.Jprogdyn;

import java.util.*;
import java.io.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import com.google.common.collect.*;
import java.util.concurrent.*;

/**
 * This class collects together some static methods for creating the first trajectory point
 * in a simulation.
 */
public class Initializer implements Immutable, Serializable
{
    // Constants

    /** For serialization. */
    public static final long serialVersionUID = 1L;

    /** Below this frequency in inverse cm, no displacements will be made. */
    public static final double DISPLACEMENT_FREQUENCY_THRESHOLD = 50.0;

    /** Imaginary frequencies will be treated as positive frequencies with this value in cm^-1. */
    public static final double MINIMUM_FREQUENCY = 2.0;

    /** If the temperature is below this number in K, modes will only be placed in the ground vibrational state. */
    public static final double MINIMUM_TEMPERATURE = 10.0;

    /** The maximum ZPE ratio (probability of being in level i vs i+1). */
    public static final double MAX_ZPE_RATIO = 0.999999;

    /** The maximum QMHO vibrational level. */
    public static final int MAX_QMHO_LEVEL = 1000;

    // Fields

    /** The original molecule. */
    public final Molecule molecule;

    /** The temperature of the simulation in K. */
    public final double temperature;

    /** The timestep of the simulation in fs. */
    public final double timestep;

    /** What kind of vibrational initialization to do. */
    public final VibrationalInitializationType vibrationType;

    /** What kind of rotational initialization to do. */
    public final RotationalInitializationType rotationType;

    /**
     * What kind of vibrational initialization to do for transition state modes (referred to by index 0...n-1 in molecule.modes).
     * Use this to turn off modes or initialize them as transition state modes.  Overrides vibrationalType.
     * Modes are zero-indexed: 0, 1, ..., n-1.
     */
    public final Map<Integer,VibrationalInitializationType> specialModeInitializationMap;

    /** The difference between the desired and actual potential energies of the initialized geometry must be within this percentage. */
    public final double tolerance;

    /** The level of theory to do the initialization at. */
    public final CalculationMethod calculationMethod;

    /** The linear scaling factor of vibrational frequencies.  This must be positive, and is used to crudely adjust for anharmonicity. */
    public final double scaleFactor;

    /** Keeps a record of the maximum possible displacement in each mode in A.  Indexed by mode 0, 1, ..., n-1. */
    public final Map<Integer,Double> maximumDisplacementRecord;

    /** Keeps a record of the displacement made in each mode in A. */
    public final Map<Integer,Double> displacementRecord;

    /** Keeps a record of the velocity in each mode. */
    public final Map<Integer,Double> velocityRecord;

    /**
     * Loads the parameters we need for initialization but doesn't actually do any work.
     * @param molecule original molecule
     * @param temperature temperature in K
     * @param timestep timestep in fs
     * @param vibrationType what kind of vibrational initialization to do by default
     * @param rotationType what kind of rotational initialization to do
     * @param specialModeInitializationMap map from mode index (0-indexed) to what kind of vibrational initialization to do in that mode 
     * @param tolerance the difference between the desired and actual potential energies must be within this percentage
     * @param calculationMethod the level of theory to do the initializaiton at
     * @param scaleFactor linear scaling factor of vibrational frequencies
     */
    public Initializer(Molecule molecule, double temperature, double timestep,
                       VibrationalInitializationType vibrationType, RotationalInitializationType rotationType,
                       Map<Integer,VibrationalInitializationType> specialModeInitializationMap,
                       double tolerance, CalculationMethod calculationMethod, double scaleFactor)
    {
        if ( molecule == null )
            throw new NullPointerException("molecule cannot be null");
        if ( molecule.modes == null || molecule.modes.size() == 0 )
            throw new NullPointerException("no modes were found for this molecule, cannot initialize");
        this.molecule = molecule;

        if ( temperature < 0.0 )
            throw new IllegalArgumentException("temperature cannot be negative");
        else if ( temperature < 10.0 )
            System.out.printf("Warning: low temperature of %.1f K selected.  Proceeding.\n", temperature);
        this.temperature = temperature;

        if ( timestep < 0.0 )
            throw new IllegalArgumentException("timestep cannot be negative");
        else if ( timestep > 10.0 )
            System.out.printf("Warning: unusually large timestep of %.1f fs selected.\n", timestep);
        else if ( timestep < 0.1 )
            System.out.printf("Warning: unusually small timestep of %.1f fs selected.\n", timestep);
        this.timestep = timestep;

        if ( vibrationType == VibrationalInitializationType.NONE && rotationType == RotationalInitializationType.NONE )
            throw new IllegalArgumentException("at least one of vibration or rotation must be initialized");
        this.vibrationType = vibrationType;
        this.rotationType = rotationType;

        if ( specialModeInitializationMap == null )
            throw new NullPointerException("use an empty map to treat all modes normally");
        for (Integer i : specialModeInitializationMap.keySet())
            {
                if ( i > molecule.modes.size() - 1 || i < 0 )
                    throw new IllegalArgumentException(String.format("special mode out of range: found %d, but valid range is 0-%d", i, molecule.modes.size()-1));
            }
        this.specialModeInitializationMap = ImmutableMap.copyOf(specialModeInitializationMap);

        if ( tolerance < 0.0 )
            throw new IllegalArgumentException("tolerance cannot be negative as we only compare abs(desired-actual) potential energy");
        this.tolerance = tolerance;

        if ( calculationMethod == null )
            throw new NullPointerException("level of theory must be specified");
        this.calculationMethod = calculationMethod;

        if ( scaleFactor < 0 )
            throw new IllegalArgumentException("cannot use negative scale factor");
        if ( scaleFactor < 0.5 || scaleFactor > 2 )
            System.out.printf("Warning: extreme vibational scale factor used: %.2f.");
        this.scaleFactor = scaleFactor;

        this.maximumDisplacementRecord = new HashMap<>(molecule.modes.size());
        this.displacementRecord = new HashMap<>(molecule.modes.size());
        this.velocityRecord = new HashMap<>(molecule.modes.size());
    }

    /**
     * Represents different kinds of possible vibrational initialization strategies.
     */
    public enum VibrationalInitializationType
    {
        /** Draw the initial conditions from a quantum mechanical canonical ensemble. */
        QUASICLASSICAL,

        /** Draw the initial conditions from a classical canonical ensemble. */
        CLASSICAL,

        /** Draw the initial displacements from a uniform classical distribution. */
        UNIFORM,

        /** Treat this as a transition state vector.  This means zero displacements and positive velocity. */
        TS_POSITIVE,

        /** Treat this as a transition state vector.  This means zero displacements and negative velocity. */
        TS_NEGATIVE,

        /** Treat this as a transition state vector.  This means zero displacements and random velocity. */
        TS_RANDOM,

        /** Do not displace this mode. */
        NONE;
    }

    /**
     * Represents different kinds of possible rotational initialization strategies.
     */
    public enum RotationalInitializationType
    {
        /** Draw the initial conditions from a classical canonical ensemble. */
        CLASSICAL,

        /** Do not add rotations. */
        NONE;
    }

    /**
     * A helper class that keeps track of some quantities we need to pass around between methods.
     * This object is mutable, but is disposed of after each initialization attempt.
     */
    private static class ScratchPaper
    {
        /** The geometry that will be perturbed (angstroms). */
        public final List<Vector3D> startingGeometry;

        /** The displacements for each atom in angstroms. */
        public final List<Vector3D> displacements;

        /** The velocities for each atom in angstroms/fs. */
        public final List<Vector3D> velocities;

        /** The desired total potential energy in hartree.  This is the energy of the unperturbed molecule plus any extra potential energy from vibrations. */
        public double desiredPotentialEnergy;

        /**
         * Start a new piece of scratch paper.
         */
        public ScratchPaper(Molecule molecule)
        {
            // copy starting geometry
            List<Vector3D> tempList = new ArrayList<>(molecule.contents.size());
            for (Atom a : molecule.contents)
                tempList.add(a.position);
            startingGeometry = ImmutableList.copyOf(tempList);

            // initialize displacements to zero
            displacements = new ArrayList<>(molecule.contents.size());
            for (Vector3D v : startingGeometry)
                displacements.add(Vector3D.ZERO);
            
            // initialize velocities to zero
            velocities = new ArrayList<>(displacements);

            // potential energy of unperturbed molecule
            desiredPotentialEnergy = molecule.potentialEnergy;
        }

        /**
         * Creates a time zero trajectory point using the displaced positions and xyz velocities that have been calculated during initialization.
         */
        public TrajectoryPoint createTrajectoryPoint()
        {
            List<Vector3D> positions = new ArrayList<>(startingGeometry);
            for (int i=0; i < positions.size(); i++)
                positions.set(i, positions.get(i).add(displacements.get(i)));
            return TrajectoryPoint.create(0.0, positions, velocities);
        }
    }

    /**
     * Makes a t=0 trajectory point for a trajectory but does not evaluate it.  Can be called multiple
     * times on the same object; repeated calls are statistically independent.  If the desired and actual
     * potential energy are within the specified tolerance, then a TrajectoryPoint will be returned.
     * @return a set of proposed positions and velocities
     */
    private ScratchPaper createCandidate()
    {
        // make a blank piece of scratch paper
        ScratchPaper scratchPaper = new ScratchPaper(molecule);
        maximumDisplacementRecord.clear();
        displacementRecord.clear();
        velocityRecord.clear();

        // perform vibrational initialization if requested
        for (int i=0; i < molecule.modes.size(); i++)
            {
                NormalMode mode = molecule.modes.get(i);
                VibrationalInitializationType type = vibrationType;
                if ( specialModeInitializationMap.containsKey(i) )
                    type = specialModeInitializationMap.get(i);

                if ( type == VibrationalInitializationType.QUASICLASSICAL )
                    doQuasiclassicalVibration(scratchPaper, mode, i);
                else if ( type == VibrationalInitializationType.CLASSICAL )
                    doClassicalVibration(scratchPaper, mode, i);
                else if ( type == VibrationalInitializationType.UNIFORM )
                    doUniformVibration(scratchPaper, mode, i);
                else if ( type == VibrationalInitializationType.TS_POSITIVE )
                    doTSVibration(scratchPaper, mode, i, VelocitySign.POSITIVE);
                else if ( type == VibrationalInitializationType.TS_NEGATIVE )
                    doTSVibration(scratchPaper, mode, i, VelocitySign.NEGATIVE);
                else if ( type == VibrationalInitializationType.TS_RANDOM )
                    doTSVibration(scratchPaper, mode, i, VelocitySign.RANDOMIZE);
                else if ( type == VibrationalInitializationType.NONE )
                    {
                        System.out.printf("Not displacing mode %d (%.0f cm^-1).\n", i, mode.frequency);
                    }
                else
                    throw new IllegalArgumentException("unsupported vibrational initialization type, should not be possible");
            }

        // perform rotational initialization if requested
        if ( rotationType == RotationalInitializationType.CLASSICAL )
            doClassicalRotations(scratchPaper);

        // construct TrajectoryPoint from scratch paper
        TrajectoryPoint candidate = scratchPaper.createTrajectoryPoint();
        double extraEnergy = (scratchPaper.desiredPotentialEnergy - molecule.potentialEnergy) * Units.KCAL_PER_HARTREE;
        System.out.printf("Desired potential energy is %.8f hartree (%.4f kcal extra).\n", scratchPaper.desiredPotentialEnergy, extraEnergy);
        return scratchPaper;
    }

    /**
     * Tries to create initial conditions for molecular dynamics.  Each proposed geometry is checked for its potential energy.
     * The potential energy must be within the tolerance specified when creating the Initializer. 
     * @param maxAttempts maximum number of tries
     * @return the zeroth trajectory point, which is guaranteed to be harmonic
     */
    public TrajectoryPoint initialize(int maxAttempts)
    {
        for (int i=0; i < maxAttempts; i++)
            {
                System.out.printf("Initialization attempt %d of a possible %d\n", i+1, maxAttempts);

                // create a candidate set of positions and velocities
                ScratchPaper scratchPaper = createCandidate();
                TrajectoryPoint candidate = scratchPaper.createTrajectoryPoint();

                // evaluate the candidate
                candidate = candidate.evaluatePoint(molecule, calculationMethod, timestep);
                
                // check the candidate for validity
                //System.out.printf("%.2f kcal of kinetic energy were added.\n", candidate.kineticEnergy * Units.KCAL_PER_HARTREE);
                double desiredPotentialEnergy = scratchPaper.desiredPotentialEnergy;  // in hartree
                double actualPotentialEnergy  = candidate.potentialEnergy;            // in hartree
                double difference = 100.0 * Math.abs((desiredPotentialEnergy - actualPotentialEnergy) / desiredPotentialEnergy);
                double absoluteDifference = (actualPotentialEnergy - desiredPotentialEnergy)*Units.KCAL_PER_HARTREE;
                double absoluteTolerance = (tolerance/100.0) * Math.abs(desiredPotentialEnergy) * Units.KCAL_PER_HARTREE;
                if ( difference > tolerance )
                    {
                        System.out.printf("Initialization failed by %.4f%% = %.1f kcal (%.4f desired PE in hartree, %.4f actual, %.1f kcal tolerance).\n",
                                          difference, absoluteDifference, desiredPotentialEnergy, actualPotentialEnergy, absoluteTolerance);
                    }
                else
                    {
                        System.out.printf("Initialization passed! (Deviation was %.4f%% = %.1f kcal, %.4f desired PE in hartree, %.4f actual, %.1f kcal tolerance)\n",
                                          difference, absoluteDifference, desiredPotentialEnergy, actualPotentialEnergy, absoluteTolerance);
                        return candidate;
                    }
            }
        throw new IllegalArgumentException(String.format("Maximum number of initialization attempts has been exceeded."));
    }

    /** How to treat velocities in each mode. */
    public enum VelocitySign
    {
        /** Set the velocity to positive. */
        POSITIVE,

        /** Set the velocity to negative. */
        NEGATIVE,

        /** Choose a random velocity sign. */
        RANDOMIZE;
    }

    /**
     * For a specified energy and displacement, add displacements to the initial geometry and
     * add vibrational kinetic energy.  The amount of PE/KE to add is decided elsewhere.
     * @param mode the normal mode to use
     * @param modeIndex the 0-indexed number of this mode
     * @param thisTotalEnergy total energy to add to this mode in kcal/mol
     * @param shift the randomly selected displacement in this mode in angstroms
     * @param scratchPaper where to update the results to
     * @param velocitySign the sign of the velocities in this mode
     */
    private void doVibration(NormalMode mode, int modeIndex, double thisTotalEnergy, double shift, ScratchPaper scratchPaper, VelocitySign velocitySign)
    {
        // calculate energies
        double reducedMass = mode.reducedMass;           // amu
        double forceConstant = mode.forceConstant;       // mDyne/A
        double frequency = scaleFactor * mode.frequency; // cm^-1

        // compute relative displacement
        double maxShift = HarmonicOscillatorDistribution.getClassicalTurningPoint(thisTotalEnergy, forceConstant);
        if ( maximumDisplacementRecord.containsKey(modeIndex) || displacementRecord.containsKey(modeIndex) || velocityRecord.containsKey(modeIndex) )
            throw new IllegalArgumentException("unexpected duplicate entry in displacement record");
        maximumDisplacementRecord.put(modeIndex, maxShift);
        double relativeDisplacement = shift/maxShift;
        if ( mode.frequency < DISPLACEMENT_FREQUENCY_THRESHOLD )
            {
                System.out.printf("Relative displacement in mode %d set from %4.0f%% to ZERO because of low frequency!\n", modeIndex, relativeDisplacement * 100.0);
                relativeDisplacement = 0.0;
            }
        displacementRecord.put(modeIndex, maxShift*relativeDisplacement);

        // compute actual displacement
        List<Vector3D> normalModeCoordinates = mode.coordinates;
        
        for (int i=0; i < normalModeCoordinates.size(); i++)
            {
                Vector3D direction = normalModeCoordinates.get(i);
                Vector3D displacement = direction.scalarMultiply(relativeDisplacement * maxShift);
                scratchPaper.displacements.set(i, displacement.add(scratchPaper.displacements.get(i)));
            }

        // update desired potential energy
        //
        // V = 0.5 * k * x^2
        //
        // dimensional analysis (note that N * m = J):
        //
        //  mDyne           1E-8 N     1E-10 m       kcal       AVOGADROS_NUMER hartree
        // ------- * A^2 * -------- * --------- * ---------- * ------------------------- 
        //    A              mDyne        A         4184 J      KCAL_PER_HARTREE kcal
        //
        // this gives the potential energy in this mode in hartree:
        double thisPotentialEnergy = 0.5 * forceConstant * shift * shift * Units.AVOGADROS_NUMBER / (1E18 * 4184 * Units.KCAL_PER_HARTREE);
        thisPotentialEnergy = thisPotentialEnergy * scaleFactor * scaleFactor; // since the force constant is proportional to frequency^2
        
        //System.out.println("PE: " + thisPotentialEnergy*Units.KCAL_PER_HARTREE);
        //System.out.println(thisTotalEnergy);
        
        // since we don't make any displacements in very low frequency or imaginary modes, they have no potential energy
        if ( relativeDisplacement == 0.0 )
            thisPotentialEnergy = 0.0;
        
        scratchPaper.desiredPotentialEnergy += thisPotentialEnergy;
        
        // update kinetic energy (hartree)
        double thisKineticEnergy = ( thisTotalEnergy / Units.KCAL_PER_HARTREE ) - thisPotentialEnergy;
        
        // compute mode velocity as sqrt(2*kineticEnergy / reduced mass), want result in A/fs
        //
        // hartree is implicitly per mol
        //
        //  hartree     J_PER_HARTREE J     kg m^2     1000 g        s^2        1E20 A^2     J_PER_HARTREE * 1E-7 A^2
        // --------- * ----------------- * -------- * -------- * ----------- * ---------- = --------------------------
        //  g / mol         hartree         s^2 J       kg        1E30 fs^2       m^2                  fs^2
        double modeVelocity = Math.sqrt(Units.J_PER_HARTREE * 1E-7 * 2.0 * thisKineticEnergy / reducedMass);

        // set sign of mode velocity
        if ( velocitySign == VelocitySign.POSITIVE )
            {
                modeVelocity = Math.abs(modeVelocity);
            }
        else if ( velocitySign == VelocitySign.NEGATIVE )
            {
                modeVelocity = -1.0 * Math.abs(modeVelocity);
            }
        else if ( velocitySign == VelocitySign.RANDOMIZE )
            {
                double randomPhase = ThreadLocalRandom.current().nextDouble();
                if ( randomPhase < 0.5 )
                    modeVelocity = Math.abs(modeVelocity) * -1.0;
                else
                    modeVelocity = Math.abs(modeVelocity);
            }
        else
            throw new IllegalArgumentException("impossible");
        velocityRecord.put(modeIndex, modeVelocity);

        // compute atom velocities in A / fs
        // velocity is a magnitude times a direction
        // the direction is the normal mode direction, the magnitude is the magnitude of the mode velocity
        for (int i=0; i < scratchPaper.velocities.size(); i++)
            {
                Vector3D currentVelocity = scratchPaper.velocities.get(i);
                Vector3D direction = normalModeCoordinates.get(i);
                Vector3D extraVelocity = direction.scalarMultiply(modeVelocity);
                Vector3D newVelocity = currentVelocity.add(extraVelocity);
                scratchPaper.velocities.set(i, newVelocity);
            }
    }

    /**
     * Performs quasi-classical initialization.  Changes will be made to the scratch paper.
     * The algorithm is as follows:
     * 1. Choose a vibrational level for each vibrational mode.  This determines the total energy in each mode.
     *    Selection is done using a Boltzmann process appropriate to the selected temperature.
     * 2. Select a random displacement in each mode by using the probability density appropriate to the vibrational
     *    eigenstate.  Displacements are limited to the classically allowed region.  Compute the relative displacement
     *    as the displacement divided by the maximum possible displacement.
     * 3. Use the relative displacements and the normal vectors for each mode to determine the actual displacements.
     * 4. Determine the desired potential energy and the kinetic energy.
     * @param scratchPaper the initialization to change
     * @param mode the mode to deal with
     * @param modeIndex the 0-indexed number of this mode
     */
    private void doQuasiclassicalVibration(ScratchPaper scratchPaper, NormalMode mode, int modeIndex)
    {
        // choose vibrational level
        double frequency = scaleFactor * mode.frequency; // cm-1
        if ( frequency < MINIMUM_FREQUENCY )
            frequency = MINIMUM_FREQUENCY;
        int level = chooseQHOVibrationalLevel(frequency, temperature);
        
        // calculate energies
        double reducedMass = mode.reducedMass;          // amu
        double forceConstant = mode.forceConstant;      // mDyne/A
        double ZPE = getZeroPointEnergy(frequency);     // kcal/mol
        double thisTotalEnergy = ZPE * ( 2.0 * level + 1.0 ); // kcal/mol (total energy for this mode)

        // get random displacement in angstroms
        double shift = HarmonicOscillatorDistribution.drawRandomQuantumDisplacement(level, reducedMass, frequency, forceConstant);
        double maxShift = HarmonicOscillatorDistribution.getClassicalTurningPoint(thisTotalEnergy, forceConstant);

        System.out.printf("Selected qc level %1d (%4.2f kcal total) for mode %4d (%4.0f cm^-1; unscaled %4.0f cm^-1).  Shift is %5.0f%% (%5.2f of a possible %5.2f A).\n",
                          level, thisTotalEnergy, modeIndex, frequency, frequency/scaleFactor, 100.0*(shift/maxShift), shift, maxShift);

        // displace and add velocities
        doVibration(mode, modeIndex, thisTotalEnergy, shift, scratchPaper, VelocitySign.RANDOMIZE);
    }

    /**
     * Performs classical initialization.  Changes will be made to the scratch paper.
     * Same algorithm as the quasi-classical algorithm (see {@link #doQuasiclassicalVibration(ScratchPaper, NormalMode)}),
     * except we choose the total energy from a continuous Boltzmann distribution and displacements from
     * the classical probability distribution.
     * @param modeIndex the 0-indexed number of this mode
     */
    private void doClassicalVibration(ScratchPaper scratchPaper, NormalMode mode, int modeIndex)
    {
        // draw the total energy for this mode from a Boltzmann distribution
        double thisTotalEnergy = RotationalBoltzmann.getRandomBoltzmannEnergy(temperature);  // kcal/mol
        double frequency = scaleFactor * mode.frequency; // cm-1
        if ( frequency < MINIMUM_FREQUENCY )
            frequency = MINIMUM_FREQUENCY;
        
        // calculate energies
        double reducedMass = mode.reducedMass;          // amu
        double forceConstant = mode.forceConstant;      // mDyne/A

        // get random displacement in angstroms
        double shift = HarmonicOscillatorDistribution.drawRandomClassicalDisplacement(thisTotalEnergy, forceConstant);
        System.out.printf("Selected classical initialization for mode %4d (%4.0f cm^-1; unscaled %4.0f cm^-1).\n", modeIndex, frequency, frequency/scaleFactor);

        // displace and add velocities
        doVibration(mode, modeIndex, thisTotalEnergy, shift, scratchPaper, VelocitySign.RANDOMIZE);
    }

    /**
     * Performs classical initialization with uniform displacement distribution.  
     * Changes will be made to the scratch paper.
     * Same algorithm as the quasi-classical algorithm (see {@link #doQuasiclassicalVibration(ScratchPaper, NormalMode)}),
     * except we choose the total energy from a continuous Boltzmann distribution and displacements from
     * a uniform distribution between turning points.
     */
    private void doUniformVibration(ScratchPaper scratchPaper, NormalMode mode, int modeIndex)
    {
        // draw the total energy for this mode from a Boltzmann distribution
        double thisTotalEnergy = RotationalBoltzmann.getRandomBoltzmannEnergy(temperature);  // kcal/mol
        double frequency = scaleFactor * mode.frequency; // cm-1
        
        // calculate energies
        double reducedMass = mode.reducedMass;          // amu
        double forceConstant = mode.forceConstant;      // mDyne/A

        // get random displacement in angstroms
        double shift = HarmonicOscillatorDistribution.drawRandomUniformDisplacement(thisTotalEnergy, forceConstant);
        System.out.printf("Selected uniform initialization for mode %4d (%4.0f cm^-1; unscaled %4.0f cm^-1).\n", modeIndex, frequency, frequency/scaleFactor);

        // displace and add velocities
        doVibration(mode, modeIndex, thisTotalEnergy, shift, scratchPaper, VelocitySign.RANDOMIZE);
    }

    /**
     * Initializes the specified mode as a transition state vector.  This means no displacement,
     * with a classical amount of kinetic energy in the positive or negative direction.
     * @param scratchPaper the initialization to perturb
     * @param mode the mode to deal with
     * @param modeIndex the 0-indexed number of this mode
     * @param velocitySign what sign to give the velocity along the imaginary mode
     */
    private void doTSVibration(ScratchPaper scratchPaper, NormalMode mode, int modeIndex, VelocitySign velocitySign)
    {
        // draw the total energy for this mode from a Boltzmann distribution
        double thisTotalEnergy = RotationalBoltzmann.getRandomBoltzmannEnergy(temperature);  // kcal/mol
        double frequency = scaleFactor * mode.frequency;

        // do not displace this mode
        double shift = 0.0;
        
        // displace and add velocities
        doVibration(mode, modeIndex, thisTotalEnergy, shift, scratchPaper, velocitySign);
        
        String velocityString = velocityRecord.get(modeIndex) > 0.0 ? "positive" : "negative";
        if ( velocityRecord.get(modeIndex) == 0.0 )
            velocityString = "zero.";
        System.out.printf("Did not displace transition state mode %4d (%4.0f cm^-1; unscaled %4.0f cm^-1).  Velocity is %s.\n",
                          modeIndex, frequency, frequency/scaleFactor, velocityString);

    }

    /**
     * Adds velocities for rotation, treating the rotations classically.
     */
    private void doClassicalRotations(ScratchPaper scratchPaper)
    {
        RotationalBoltzmann rotationalBoltzmann = new RotationalBoltzmann(molecule.contents);
        double kT = Units.BOLTZMANN_EV * temperature; // temperature in electron volts
        double x1 = ThreadLocalRandom.current().nextDouble(-1.0,1.0);
        double x2 = ThreadLocalRandom.current().nextDouble(-1.0,1.0);
        double x3 = ThreadLocalRandom.current().nextDouble(-1.0,1.0);
        Vector3D omega = rotationalBoltzmann.getOmega(kT, x1, x3, x3); // this is the axis of rotation
        List<Vector3D> velocities = RotationalBoltzmann.getVelocities(molecule.contents, omega);
        Vector3D momentum = Vector3D.ZERO;
        for (int i=0; i < scratchPaper.velocities.size(); i++)
            {
                Vector3D currentVelocity = scratchPaper.velocities.get(i);
                Vector3D extraVelocity = velocities.get(i);
                Vector3D newVelocity = currentVelocity.add(extraVelocity);
                scratchPaper.velocities.set(i, newVelocity);
                //System.out.printf("Rotational velocity norm in atom %3d is %.3f mA/fs.\n", i, extraVelocity.getNorm()*1000);
                double mass = molecule.contents.get(i).mass;
                momentum = momentum.add(extraVelocity.scalarMultiply(mass)); 
            }
        //System.out.println("Sum of rotational momenta: " + momentum);
    }

    /**
     * Helper method that converts a frequency to its corresponding quantum harmonic oscillator zero point energy.
     * @param frequency the frequency of the mode in cm^-1
     * @return the zero-point energy in kcal per mol
     */
    public static double getZeroPointEnergy(double frequency)
    {
        // prevent negative energies
        if ( frequency < 0.0 )
            throw new IllegalArgumentException("negative frequencies not allowed");

        // calculate the ZPE
        //   ZPE = 0.5 * hbar * angular frequency
        //       = 0.5 * planck's constant * frequency
        //       = 0.5 * planck's constant * speed of light * wavenumber
        // units =            J * s              cm / s         1 / s    = J, so convert answer to kcal/mol
        double ZPE = 0.5 * Units.H * Units.C * frequency * Units.J_TO_KCAL_PER_MOL;
        return ZPE;
    }

    /**
     * Chooses a random vibrational level for a quantum harmonic oscillator.  The level is chosen
     * from a Boltzmann distribution appropriate to the specified temperature.  This method is
     * thread safe.
     * @param frequency the frequency of this mode in cm^-1
     * @param temperature the temperature in K
     * @return the vibrational level (0, 1, ...)
     */
    public static int chooseQHOVibrationalLevel(double frequency, double temperature)
    {
        // if the temperature is low, put everything in the ground state
        if (temperature < MINIMUM_TEMPERATURE)
            return 0;

        double ZPE = getZeroPointEnergy(frequency); // ZPE in kcal/mol

        // the probability of being in level i vs level i+1 is exp(-spacing/RT),
        // where the level spacings are equal (spacing = 2 * ZPE)
        double ZPEratio = Math.exp( (-2.0 * ZPE) / (Units.R_GAS_KCAL * temperature) );

        if (ZPEratio > MAX_ZPE_RATIO)
            ZPEratio = MAX_ZPE_RATIO;

        // the probability of being in all states must sum to 1:
        // 1 = Pr(n=0) + Pr(n=1) + ... --> call P the probablity of being in the ground state
        // 1 = P + P*ZPEratio + P*ZPEratio^2 + ... --> factor out geometric series
        // 1 = P / (1-ZPEratio)
        // P = 1-ZPEratio
        double P = 1.0 - ZPEratio;

        // get a new random number to see which level to put this mode in
        double randomNumber = ThreadLocalRandom.current().nextDouble();
        int thisLevel = 0;
        while (true)
            {
                if (randomNumber < P)
                    return thisLevel;
                thisLevel++;
                
                // random number has to be...
                // Pr(n=0) = P --> [0,P)
                // Pr(n=1) = P*ZPEratio --> [P,P + P*ZPEratio)
                // Pr(n=2) = P*ZPEratio^2 --> [P + P*ZPEratio, P + P*ZPEratio + P*ZPEratio^2)
                // etc.
                P = P + P * ZPEratio;

                // only allow a certain number of excitations
                if (thisLevel > MAX_QMHO_LEVEL)
                    return MAX_QMHO_LEVEL;
            }
    }

    @Override
    public String toString()
    {
        String returnString = String.format("Initializer for: %s\n", molecule.name);
        returnString       += String.format("Temperature = %.1f K, timestep = %.1f fs, tolerance = %.1f percent\n", temperature, timestep, tolerance);
        returnString       += String.format("Level of theory: %s\n", calculationMethod.toString());
        returnString       += String.format("Scale factor: %.2f\n", scaleFactor);
        returnString       += String.format("Vibrational initialization: %s\n", vibrationType.toString());
        returnString       += String.format("Rotational initialization: %s\n", rotationType.toString());
        if ( specialModeInitializationMap.size() == 0 )
            returnString   +=               "No modes are treated specially.";
        else
            {
                for (Integer i : specialModeInitializationMap.keySet())
                    returnString   += String.format("   Special Mode %d: %s\n", i, specialModeInitializationMap.get(i).toString());
            }
        return returnString;
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(molecule, temperature, timestep, vibrationType, rotationType,
                            specialModeInitializationMap, tolerance, calculationMethod, scaleFactor);
    }

    @Override
    public boolean equals(Object obj)
    {
        if ( obj == null ) 
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof Initializer) )
            return false;

        Initializer i = (Initializer)obj;
        if ( Objects.equals(molecule, i.molecule) &&
             temperature == i.temperature &&
             timestep == i.timestep &&
             vibrationType == i.vibrationType &&
             rotationType == i.rotationType &&
             Objects.equals(specialModeInitializationMap, i.specialModeInitializationMap) &&
             tolerance == i.tolerance &&
             Objects.equals(calculationMethod, i.calculationMethod) &&
             scaleFactor == i.scaleFactor )
            return true;
        return false;
    }
}
