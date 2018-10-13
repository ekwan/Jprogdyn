package edu.harvard.chemistry.ekwan.Jprogdyn;

import java.io.*;
import java.util.*;

/**
 * This abstract class represents files that are used as inputs to other programs.
 * This class is immutable.
 */
public abstract class InputFileFormat implements FileFormat, Immutable, Serializable {

    /** For serialization. */
    public static final long serialVersionUID = 1L;

    /** The string that will be written to a file */
    public final String stringRepresentation;

    /**
     * Constructor.
     * @param stringRperesentation the text of the file
     */
    public InputFileFormat(String stringRepresentation) {
        this.stringRepresentation = stringRepresentation;
    }

    /**
     * Writes the input file to disk.
     * @param filename the destination file
     */
    public void write(String filename) {
        writeStringToDisk(stringRepresentation, filename);
	}

    /**
     * Convenience method that writes a string to a file.
     * @param string the string to write
     * @param filename the filename to write to
     */
    public static void writeStringToDisk(String string, String filename) {
        try (PrintWriter outputFile = new PrintWriter(filename)) {
                outputFile.print(string);
            }
        catch (IOException e) {
                // abort program if there's a problem
                System.out.println("Error writing to " + filename + "!");
                e.printStackTrace();
            }
    }

    /**
     * Convenience method that appends a string to a file.
     * @param string the text to append
     * @param filename the file to append to
     */
    public static void appendStringToDisk(String string, String filename) {
        try {
                File file = new File(filename);
                if ( ! file.exists() )
                    file.createNewFile();
                FileWriter fileWriter = new FileWriter(filename,true);
                BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);
                bufferedWriter.write(string);
                bufferedWriter.close();
            }
        catch (IOException e) {
                System.out.println("Error appending to " + filename + "!");
                e.printStackTrace();
            }
    }

    @Override
    public String toString() {
        return stringRepresentation;
    }
}
