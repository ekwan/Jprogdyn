package edu.harvard.chemistry.ekwan.Jprogdyn;

import java.util.*;
import java.lang.Math;
import org.apache.commons.math3.geometry.euclidean.threed.*;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * This test takes a single normal mode from a molecule, creates positive and negative
 * displacements at regular intervals along that mode (up to the classical turning point
 * defined by the zero-point energy), and writes the results to a set of files.
 *
 * To run: mvn -Dtest=HarmonicTest test
 */
public class HarmonicTest extends TestCase
{
    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public HarmonicTest( String testName )
    {
        super( testName );
    }

    /**
     * 
     */
    public void testHarmonic()
    {
        // displacement parameters
        String moleculeFilename = "test_files/methane_b3lyp_midix.out";        // filename to read modes from
        int modeIndex = 1;                                                     // zero-indexed
        double relativeShiftMin = -1.0;                                        // most negative fraction of classical turning point to displace by
        double relativeShiftMax = 1.0;                                         // most positive fraction of classical turning point to displace by
        double relativeShiftInterval = 0.1;                                    // the displacement step interval in fractions of the classical turning point
        if ( relativeShiftMin >= relativeShiftMax )
            throw new IllegalArgumentException("check shift min/max");
        if ( relativeShiftInterval < 0 )
            throw new IllegalArgumentException("check shift interval");
        
        // Gaussian parameters
        String outputPrefix = "analysis/methane";
        String outputSuffix = "b3lyp-midix";
        int processors = 4;
        int memory = 3; // in GB
        String routeCard = "#p b3lyp midix";
        String footer = "\n";

        // read molecule
        System.out.println("Loading data...\n");
        GaussianOutputFile frequenciesOutputFile = new GaussianOutputFile(moleculeFilename);
        Molecule frequenciesMolecule = frequenciesOutputFile.molecule;
        System.out.printf("Molecule read from %s:\n", moleculeFilename);
        System.out.println(frequenciesMolecule);
        System.out.println();
		
        // read normal modes
        NormalMode mode = frequenciesMolecule.modes.get(modeIndex);
        List<Vector3D> normalModeCoordinates = mode.coordinates;
        System.out.printf("Mode %d:\n", modeIndex+1);
        System.out.println(mode);
        
        // get normal mode information
        double reducedMass = mode.reducedMass;           			// amu
        double forceConstant = mode.forceConstant;       			// mDyne/A
        double frequency = mode.frequency;			                // cm^-1
        double ZPE = Initializer.getZeroPointEnergy(frequency);     // kcal/mol
		double maxShift = HarmonicOscillatorDistribution.getClassicalTurningPoint(ZPE, forceConstant);  // Angstroms

        // generate input files
        for (double relativeShift = relativeShiftMin; relativeShift <= relativeShiftMax; relativeShift += relativeShiftInterval) {
            // calculate displacements
            List<Vector3D> finalPositions = new ArrayList<>(frequenciesMolecule.contents.size());
            for (int i=0; i < normalModeCoordinates.size(); i++) {
                Vector3D direction = normalModeCoordinates.get(i);
                Vector3D currentPosition = frequenciesMolecule.contents.get(i).position;
                Vector3D displacement = direction.scalarMultiply(relativeShift * maxShift);
                Vector3D finalPosition = currentPosition.add(displacement);
                finalPositions.add(finalPosition);
            }

            // write out molecule
            StringBuilder s = new StringBuilder();
            s.append(String.format("%%mem=%dGB\n", memory));
            s.append(String.format("%%nprocshared=%d\n", processors));
            s.append(routeCard + String.format("\n\ndisplacement: %.6f A\n\n", relativeShift * maxShift));
            s.append(String.format("%d %d\n", frequenciesMolecule.charge, frequenciesMolecule.multiplicity));
            for (int i=0; i < frequenciesMolecule.contents.size(); i++) {
                String symbol = frequenciesMolecule.contents.get(i).symbol;
                double x = finalPositions.get(i).getX();
                double y = finalPositions.get(i).getY();
                double z = finalPositions.get(i).getZ();
                s.append(String.format("   %-5s     %15.10f    %15.10f    %15.10f\n", symbol, x, y, z));
            }
            if ( footer.trim().length() > 0 )
            	s.append(footer + "\n");
        	s.append("\n\n");
            String shiftString = String.format("%03.0f", Math.abs(relativeShift) * 100.0);
            if ( relativeShift < 0 && Math.abs(relativeShift) > 0.0001 )
                shiftString = "m" + shiftString;
            String filename = String.format("%s-%s-%s.gjf", outputPrefix, shiftString, outputSuffix);
            InputFileFormat.writeStringToDisk(s.toString(),filename);
            System.out.printf("Wrote to %s (absolute displacement=%.6f A, relative displacement = %.3f).\n", filename, relativeShift * maxShift, relativeShift);
        }

        assertTrue( true );
    }
}
