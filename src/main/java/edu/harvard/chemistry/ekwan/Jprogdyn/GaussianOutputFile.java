package edu.harvard.chemistry.ekwan.Jprogdyn;

import java.io.*;
import java.util.*;
import com.google.common.collect.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;

/**
 * This class reads Gaussian output files.  It reads the last energy and geometry that was
 * printed out.  If calculations have been performed at more than one level of theory, the
 * last energy for the first level of theory is used.  If freq=hpmodes is requested, the
 * normal modes are read.  If an NMR calculation is performed, the shieldings are read.
 * Atomic weights are also read (#p should be requested).  Multi-line title cards are
 * supported.  Charge and multiplicity are read.
 *
 * To prevent problems with energies from different levels of theory, Gaussian output files
 * that use Link1 should be used with caution.
 */
public class GaussianOutputFile extends OutputFileFormat {

    /** The molecule. */
    public final Molecule molecule;

    /**
     * Reads a Gaussian output file from disk.  Reads the last geometry and potential energy, as well
     * as forces and normal modes, if they are available.  Verbose output should be called with
     * #p so that the masses can be read.
     * @param filename the Gaussian output file to read from
     */
    public GaussianOutputFile(String filename)
    {
        super(filename);

        // check for Link1
        int link1count = 0;
        for (int i=0; i < fileContents.size(); i++) {
            List<String> fields = fileContents.get(i);
            if ( fields.size() > 4 && fields.get(0).equals("Entering") && fields.get(1).equals("Link") && fields.get(2).equals("1") )
                link1count++;
            if ( link1count > 1 )
                {
                    if ( Loader.getString("trajectory_type").equals("reaction") )
                        System.out.printf("Warning: %s has a Link1 directive, continuing anyways.\n", filename);
                    break;
                }
        }

        // read atom symbols
        List<String> atomSymbols = new ArrayList<>();
        for (int i=0; i < fileContents.size(); i++) {
            List<String> fields = fileContents.get(i);
            if ( fields.size() == 6 && fields.get(0).equals("Charge") && fields.get(3).equals("Multiplicity") ) {
                while (true) {
                    i++;
                    fields = fileContents.get(i);
                    if ( fields.size() < 4 || fields.size() > 5 )
                        break;
                    atomSymbols.add(fields.get(0));
                }
            break;
            }
        }

        // read atom weights
        List<Double> atomicWeights = new ArrayList<>(atomSymbols.size()); // in amu
        for (int i=0; i < fileContents.size(); i++) {
            List<String> fields = fileContents.get(i);
            if ( fields.size() > 1 && fields.get(0).equals("AtmWgt=") ) {
                if ( atomicWeights.size() == atomSymbols.size() )
                    break;
                else if ( atomicWeights.size() > atomSymbols.size() )
                    throw new IllegalArgumentException("too many weights");

                for (int j=1; j < fields.size(); j++)
                    atomicWeights.add(Double.valueOf(fields.get(j)));
            }
        }
        if ( atomicWeights.size() > 0 && atomicWeights.size() != atomSymbols.size() )
            throw new IllegalArgumentException("error reading atomic weights");

        // read last geometry
        List<Vector3D> positions = new ArrayList<>(atomSymbols.size());
        for (int i=0; i < fileContents.size(); i++) {
            List<String> fields = fileContents.get(i);
            if ( fields.size() == 2 && fields.get(0).equals("Standard") && fields.get(1).equals("orientation:") ) {
                i += 5;
                positions.clear();
                for (int j=i; j < i+atomSymbols.size(); j++) {
                    fields = fileContents.get(j);
                    int numberOfFields = fields.size();
                    if ( numberOfFields < 5 || numberOfFields > 6 )
                        throw new IllegalArgumentException("unexpected number of geometry fields:\n" + fields.toString());
                    double x = Double.valueOf(fields.get(numberOfFields-3));
                    double y = Double.valueOf(fields.get(numberOfFields-2));
                    double z = Double.valueOf(fields.get(numberOfFields-1));
                    positions.add(new Vector3D(x,y,z));
                }
                i += atomSymbols.size();
            }
        }
        if ( positions.size() == 0 || positions.size() != atomSymbols.size() )
            throw new IllegalArgumentException("error reading positions!");

        // construct atoms
        List<Atom> contents = new ArrayList<>();
        for (int i=0; i < atomSymbols.size(); i++) {
            String symbol = atomSymbols.get(i);
            double mass = 0.0;
            if ( atomicWeights.size() > 0 )
                mass = atomicWeights.get(i);
            Vector3D position = positions.get(i);
            Atom atom = new Atom(symbol, mass, position);
            contents.add(atom);
        }

        // read normal modes if available (assumes freq=hpmodes is set)
        List<Double> frequencies = new ArrayList<>();
        List<Double> reducedMasses = new ArrayList<>();
        List<Double> forceConstants = new ArrayList<>();
        List<List<Vector3D>> normalDisplacements = new ArrayList<>();
        for (int i=0; i < fileContents.size(); i++) {
            List<Integer> skipJ = new ArrayList<>();
            List<String> fields = fileContents.get(i);
            if ( fields.size() > 2 && fields.get(0).equals("Frequencies") && fields.get(1).equals("---") ) {
                // read frequencies in inverse cm
                for (int j=2; j < fields.size(); j++) {
                    if ( fields.get(j).indexOf("*") > -1 ) {
                            System.out.println("Warning: invalid frequency ignored.");
                            skipJ.add(j+1);  // have to account for the fact that frequencies is on word
                                             // but reduced masses and force constants are two words
                        }
                    else
                        frequencies.add(Double.valueOf(fields.get(j)));
                }

                // read reduced masses in amu
                i++;
                fields = fileContents.get(i);
                for (int j=3; j < fields.size(); j++)
                    if ( ! skipJ.contains(j) )
                        reducedMasses.add(Double.valueOf(fields.get(j)));

                // read force constants in mDyne/A
                i++;
                fields = fileContents.get(i);
                for (int j=3; j < fields.size(); j++)
                    if ( ! skipJ.contains(j) )
                        forceConstants.add(Double.valueOf(fields.get(j)));
                
                // read normal displacements in A
                i += 3;
                fields = fileContents.get(i);
                List<Integer> skipModeIndices = new ArrayList<>();
                double[][][] temp = new double[fields.size()-3][atomSymbols.size()][3]; // mode index, atom index, xyz index
                for (int j=i; j < i + atomSymbols.size()*3; j++) {
                    fields = fileContents.get(j);
                    int coordIndex = Integer.parseInt(fields.get(0))-1;
                    int atomIndex  = Integer.parseInt(fields.get(1))-1;
                    for (int k=3; k < fields.size(); k++) {
                        int modeIndex = k-3;
                        if ( skipJ.contains(k) ) {
                            skipModeIndices.add(modeIndex);
                            continue;
                        }
                        
                        temp[modeIndex][atomIndex][coordIndex] = Double.parseDouble(fields.get(k));
                    }
                }
                for (int modeIndex = 0; modeIndex < temp.length; modeIndex++) {
                    if ( skipModeIndices.contains(modeIndex) )
                        continue;

                    double[][] array = temp[modeIndex];
                    List<Vector3D> displacements = new ArrayList<>();
                    for (int atomIndex = 0; atomIndex < array.length; atomIndex++) {
                        Vector3D v = new Vector3D(array[atomIndex][0], array[atomIndex][1], array[atomIndex][2]);
                        displacements.add(v);
                    }
                    normalDisplacements.add(displacements);
                }
            }
        }

        if ( frequencies.size() != reducedMasses.size() || reducedMasses.size() != forceConstants.size() ||
             forceConstants.size() != normalDisplacements.size() ) {
            System.out.printf("frequencies: %d, red masses: %d, force constants: %d, displacements: %d\n",
                              frequencies.size(), reducedMasses.size(), forceConstants.size(), normalDisplacements.size());
            throw new IllegalArgumentException("size mismatch in normal mode parser");
        }
        List<NormalMode> modes = new ArrayList<>();
        for (int i=0; i < frequencies.size(); i++)
            modes.add( new NormalMode(frequencies.get(i), reducedMasses.get(i), forceConstants.get(i), normalDisplacements.get(i)) );

        // read forces if available
        List<Vector3D> forces = new ArrayList<>();
        for (int i=0; i < fileContents.size(); i++) {
            List<String> fields = fileContents.get(i);
            if ( fields.size() == 4 && fields.get(2).equals("Forces") && fields.get(3).equals("(Hartrees/Bohr)") ) {
                i += 3;
                for (int j=i; j < i+atomSymbols.size(); j++) {
                    fields = fileContents.get(j);
                    double x = Double.valueOf(fields.get(2));
                    double y = Double.valueOf(fields.get(3));
                    double z = Double.valueOf(fields.get(4));
                    Vector3D force = new Vector3D(x,y,z);
                    forces.add(force);
                    //System.out.println(force);
                }
                break;
            }
        }
        if ( forces.size() > 0 && forces.size() != atomSymbols.size() )
            throw new IllegalArgumentException("unexpected number of forces");

        // read name
        String name = "default";
        String[] lines = stringRepresentation.split("[\r\n]+");
        int nameStartIndex = 0;
        int nameEndIndex = 0;
        for (int i=0; i < lines.length; i++) {
            String line = lines[i];
            if ( line.trim().equals("Symbolic Z-matrix:") ) {
                if ( lines[i-1].startsWith("--") && lines[i-1].endsWith("--") )
                    nameEndIndex = i-1;
                else {
                    System.out.println(lines[i-1]);
                    throw new IllegalArgumentException("error parsing name");
                }

                int j = nameEndIndex - 1;
                while (j >= 0) {
                    if ( lines[j].startsWith("--") && lines[j].endsWith("--") ) {
                        nameStartIndex = j+1;
                        break;
                    }
                    j--;
                }
                if ( nameStartIndex == 0 )
                    throw new IllegalArgumentException("error parsing name");

                name = lines[nameStartIndex];
                for (j=nameStartIndex+1; j < nameEndIndex; j++)
                    name += "\n" + lines[j];

                break;
            }
        }

        // read charge and multiplicity
        int charge = 0;
        int multiplicity = 1;
        for (List<String> fields : fileContents) {
            if ( fields.size() == 6 && fields.get(0).equals("Charge") && fields.get(3).equals("Multiplicity") ) {
                charge = Integer.parseInt(fields.get(2));
                multiplicity = Integer.parseInt(fields.get(5));
                break;
            }
        }

        // read potential energy
        // read the last energy of the first entry into link1, so that the NMR calculation can be done at a different level of theory
        double potentialEnergy = 0.0;
        link1count = 0;
        for (List<String> fields : fileContents) {
            if ( fields.size() > 4 && fields.get(0).equals("Entering") && fields.get(1).equals("Link") && fields.get(2).equals("1") )
                link1count++;
            if ( link1count >= 2 )
                break;
            if ( fields.size() > 5 && fields.get(0).equals("SCF") && fields.get(1).equals("Done:") )
                potentialEnergy = Double.parseDouble(fields.get(4));
        }

        // read shieldings if available
        List<Double> shieldings = new ArrayList<>();
        for (List<String> fields : fileContents) {
            if ( fields.size() > 5 && fields.get(2).equals("Isotropic") ) {
                Double shielding = Double.valueOf(fields.get(4));
                shieldings.add(shielding);

                if ( shieldings.size() > atomSymbols.size() )
                    throw new IllegalArgumentException("too many shieldings parsed");
            }
        }

        // create object
        this.molecule = new Molecule(contents, modes, forces, name, charge, multiplicity, potentialEnergy, shieldings);
    }

    @Override 
    public String toString() {
        return molecule.toString();
    }
}
