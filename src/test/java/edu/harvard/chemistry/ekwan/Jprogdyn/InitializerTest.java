package edu.harvard.chemistry.ekwan.Jprogdyn;

import java.util.*;
import java.lang.Math;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import com.google.common.collect.*;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * This test generates thermal initializations from a file and writes the
 * geometries to a set of files.  These files are intended as input for
 * other programs like Gaussian.
 *
 * The first structure will be unperturbed.
 *
 * To run: mvn -Dtest=InitializerTest test
 */
public class InitializerTest extends TestCase
{
    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public InitializerTest( String testName )
    {
        super( testName );
    }

    /**
     * 
     */
    public void testInitializer()
    {
        // displacement parameters
        String moleculeFilename = "test_files/simple-bridging-ts-b3lyp_d3bj-juntz-gas.out";      // filename to read modes from
        double temperature = 298.0;                                                              // in K
        int numberOfInitializations = 100;                                                       // how many files to make
        Map<Integer,Initializer.VibrationalInitializationType> specialModeInitializationMap = new HashMap<>();  // zero-indexed mode number --> initialization type
        specialModeInitializationMap.put(0,Initializer.VibrationalInitializationType.NONE);      // don't displace the TS mode

        // Gaussian parameters
        String outputPrefix = "analysis/simple-bridging-ts-b3lyp_d3bj-juntz-gas";
        String footer = "\n";

        // read molecule
        System.out.println("Loading data...\n");
        GaussianOutputFile frequenciesOutputFile = new GaussianOutputFile(moleculeFilename);
        Molecule molecule = frequenciesOutputFile.molecule;
        System.out.printf("Molecule read from %s:\n", moleculeFilename);
        System.out.println(molecule);
        System.out.println();
		System.out.println("Normal modes:");
        for (int i=0; i < molecule.modes.size(); i++) {
            NormalMode mode = molecule.modes.get(i);
            System.out.printf("Mode %4d : %.0f cm-1\n", i, mode.frequency);
        }

        // make dummy dynamics method
        CalculationMethod dynamicsMethod = new GaussianCalculationMethod(CalculationMethod.CalculationType.ENERGY_AND_FORCE,
                                                                         3, 4, "#p", footer);

        // make Initializer object
        Initializer initializer = new Initializer(molecule,                                                  // the molecule with frequencies to initialize with
												  temperature,                                               // in K
												  1.0, 			                                             // timestep in fs (irrelevant here)
												  Initializer.VibrationalInitializationType.QUASICLASSICAL,  // default vibrational initialization type
                                                  Initializer.RotationalInitializationType.NONE,             // no need for rotations here
                                                  specialModeInitializationMap,                              // treat some modes differently as specified here
                                                  0.1,                                                       // harmonic tolerance in percent (irrelevant)
                                                  dynamicsMethod,                                            // dynamics method (irrelevant)
                                                  1.0);                                                      // frequency scaling factor (irrelevant)

        // generate input files
        // first structure will be unperturbed
        for (int i=0; i < numberOfInitializations; i++) {
            System.out.printf("\n>>> Iteration %3d of %3d <<<\n", i+1, numberOfInitializations);

            // generate a new initialization
            Molecule newMolecule = molecule;
            if ( i>0 )
                newMolecule = initializer.generateStructure(molecule);
            else
                System.out.println("[ Not perturbing molecule for first iteration. ]");

            // write out molecule
            StringBuilder s = new StringBuilder();
            s.append(String.format("%d %d\n", molecule.charge, molecule.multiplicity));
            for (int j=0; j < molecule.contents.size(); j++) {
                String symbol = newMolecule.contents.get(j).symbol;
                double x = newMolecule.contents.get(j).position.getX();
                double y = newMolecule.contents.get(j).position.getY();
                double z = newMolecule.contents.get(j).position.getZ();
                s.append(String.format("   %-5s     %15.10f    %15.10f    %15.10f\n", symbol, x, y, z));
            }
            String filename = String.format("%s-init_%03d.template", outputPrefix, i);
            InputFileFormat.writeStringToDisk(s.toString(),filename);
            System.out.printf("> Wrote to %s.\n\n", filename);
        }

        assertTrue( true );
    }
}
