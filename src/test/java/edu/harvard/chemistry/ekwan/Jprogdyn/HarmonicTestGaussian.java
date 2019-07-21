package edu.harvard.chemistry.ekwan.Jprogdyn;

import java.util.*;
import java.lang.Math;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import java.nio.file.*;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * This test takes a single normal mode from a molecule, creates positive and negative
 * displacements at regular intervals along that mode (up to the classical turning point
 * defined by the zero-point energy), and writes the results to a single Gaussian input file
 * as a set of --Link1-- commands.
 *
 * For imaginary frequencies, there is no classical turning point, so displacements are made
 * directly in angstroms.
 *
 * To run: mvn -Dtest=HarmonicTestGaussian test
 */
public class HarmonicTestGaussian extends TestCase
{
    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public HarmonicTestGaussian( String testName )
    {
        super( testName );
    }

    /**
     * Read the molecule, create displacements, then write the results to one Gaussian input file.
     */
    public void testHarmonicGaussian()
    {
        // displacement parameters
        int modeIndex = 0;                                                     // zero-indexed
        
        double relativeShiftMin = -1.0;                                        // most negative fraction of classical turning point to displace by
        double relativeShiftMax = 1.0;                                         // most positive fraction of classical turning point to displace by
        double relativeShiftInterval = 0.1;                                    // the displacement step interval in fractions of the classical turning point
        
        double absoluteShiftMin = -0.005;                                        // used for negative frequencies
        double absoluteShiftMax = 0.005;                                         // displacement in angstroms
        double absoluteShiftInterval = 0.0005;
        
        // check parameters
        if ( relativeShiftMin >= relativeShiftMax )
            throw new IllegalArgumentException("check relative shift min/max");
        if ( relativeShiftInterval < 0.0 )
            throw new IllegalArgumentException("check relative shift interval");
        if ( absoluteShiftMin >= absoluteShiftMax )
            throw new IllegalArgumentException("check absolute shift min/max");
        if ( absoluteShiftInterval < 0.0 )
            throw new IllegalArgumentException("check absolute shift interval");
         
        // Gaussian parameters
        String outputDirectory = "analysis";
        int processors = 36;
        int memory = 24;           // in GB
        String footer = "\n";

        // loop through all matching filenames
        //String path = "/Users/ekwan/research/grid/benchmark/simple-anionic";
        //String glob = "simple-anionic-ts-*-631+gd-gas.out";
        String path = "test_files";
        String glob = "simple*ts*m062x*631+gd*.out";
        try ( DirectoryStream<Path> dirStream = Files.newDirectoryStream( Paths.get(path), glob)) {
            for (Path p : dirStream) {
                // read molecule
                String moleculeFilename = p.toString();
				System.out.printf("\n\nReading data from %s...", moleculeFilename);                
        		GaussianOutputFile frequenciesOutputFile = new GaussianOutputFile(moleculeFilename);
        		Molecule frequenciesMolecule = frequenciesOutputFile.molecule;
				System.out.println("done.");	
    			String routeCard = frequenciesOutputFile.routeCard;        
	            System.out.printf("Original Route: %s\n", routeCard);
                routeCard = routeCard.replaceFirst("opt\\S*\\s+?","").replaceFirst("freq\\S*\\s+?","").replaceFirst("#p","#t");
                routeCard = routeCard.replaceFirst("opt\\S*$","").replaceFirst("freq\\S*$","");
                System.out.printf("Modified Route: %s\n", routeCard);

                // read normal modes
                NormalMode mode = frequenciesMolecule.modes.get(modeIndex);
                List<Vector3D> normalModeCoordinates = mode.coordinates;
                
                // get normal mode information
                double reducedMass = mode.reducedMass;           			// amu
                double forceConstant = mode.forceConstant;       			// mDyne/A
                double frequency = mode.frequency;			                // cm^-1
                double ZPE = -1.0;
                double maxShift = -1.0;
                System.out.printf("Will use mode %d, which has a frequency of %.1f cm-1.\n", modeIndex+1, mode.frequency);
                
                // determine shift bounds
                double thisShift = 0.0;
                double thisMaxShift = 0.0;
                double thisShiftStep = 0.0;
                double ZPEMaxShift = 0.0;
                if ( frequency < 0.0 ) {
                    thisShift = absoluteShiftMin;
                    thisMaxShift = absoluteShiftMax;
                    thisShiftStep = absoluteShiftInterval;
                }
                else {
                    thisShift = relativeShiftMin;
                    thisMaxShift = relativeShiftMax;
                    thisShiftStep = relativeShiftInterval;
                    ZPEMaxShift = HarmonicOscillatorDistribution.getClassicalTurningPoint(ZPE, forceConstant); // in Angstroms
                }
                
                // generate output file header
				StringBuilder outputStringBuilder = new StringBuilder();
                int numberOfGeometries = 0;
                for (; thisShift <= thisMaxShift+0.00001; thisShift += thisShiftStep) {
                    // calculate displacements
                    numberOfGeometries++;
                    List<Vector3D> finalPositions = new ArrayList<>(frequenciesMolecule.contents.size());
                    for (int i=0; i < normalModeCoordinates.size(); i++) {
                        Vector3D direction = normalModeCoordinates.get(i);
                        Vector3D currentPosition = frequenciesMolecule.contents.get(i).position;
                        double scalar = frequency < 0.0 ? thisShift : thisShift * ZPEMaxShift;
                        Vector3D displacement = direction.scalarMultiply(scalar);
                        Vector3D finalPosition = currentPosition.add(displacement);
                        finalPositions.add(finalPosition);
                    }

                    boolean firstGeometry = outputStringBuilder.length() == 0;
                    if ( ! firstGeometry )
                        outputStringBuilder.append("--Link1--\n");
                    outputStringBuilder.append("%chk=checkpoint.chk\n");
                    outputStringBuilder.append(String.format("%%mem=%dGB\n", memory));
                    outputStringBuilder.append(String.format("%%nprocshared=%d\n", processors));
                    if ( firstGeometry )
                        outputStringBuilder.append(routeCard + "\n");
                    else
                        outputStringBuilder.append(routeCard + " guess=read\n");
                    if ( frequency < 0.0 )
                        outputStringBuilder.append(String.format("\ndisplacement: %.6f A\n\n", thisShift));
                    else
                        outputStringBuilder.append(String.format("\ndisplacement: %.6f A\n\n", thisShift * ZPEMaxShift));
                    outputStringBuilder.append(String.format("%d %d\n", frequenciesMolecule.charge, frequenciesMolecule.multiplicity));
                    for (int i=0; i < frequenciesMolecule.contents.size(); i++) {
                        String symbol = frequenciesMolecule.contents.get(i).symbol;
                        double x = finalPositions.get(i).getX();
                        double y = finalPositions.get(i).getY();
                        double z = finalPositions.get(i).getZ();
                        outputStringBuilder.append(String.format("   %-5s     %15.10f    %15.10f    %15.10f\n", symbol, x, y, z));
                    }
                    if ( footer.trim().length() > 0 )
                        outputStringBuilder.append(footer + "\n");
                    outputStringBuilder.append("\n\n");
                }

				// write out file
                String outputFilename = p.getFileName().toString().replaceFirst(".out",".gjf").replaceFirst("-ts-","-displacements-");
                String outputFullFilename = String.format("%s/%s", outputDirectory, outputFilename);
                InputFileFormat.writeStringToDisk(outputStringBuilder.toString(),outputFullFilename);
                System.out.printf(">>> Wrote %d geometries to %s.\n", numberOfGeometries, outputFullFilename);
	
/*
                if ( frequency > 0.0 ) {
                    ZPE = Initializer.getZeroPointEnergy(frequency);     // kcal/mol
                    maxShift = HarmonicOscillatorDistribution.getClassicalTurningPoint(ZPE, forceConstant);  // Angstroms
                
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
                }
				else {
				}
                */
            }
        }
        catch (Exception e) { e.printStackTrace(); }
        assertTrue( true );

/*
		

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
*/
    }
}
