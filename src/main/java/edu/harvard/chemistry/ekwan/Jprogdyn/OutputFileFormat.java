package edu.harvard.chemistry.ekwan.Jprogdyn;

import java.io.*;
import java.util.*;
import com.google.common.collect.*;

/**
 * This abstract class represents the result of a program external to Jprogdyn.
 */
public abstract class OutputFileFormat implements FileFormat
{

    /**
     * The raw text in the file.
     */
    public final String stringRepresentation;

    /**
     * The parsed contents of the file.<p>
     * The outer list contains each line.<p>
     * Each inner list contains the space-separated tokens for each line.<p>
     */
    public final List<List<String>> fileContents;

    /**
     * Constructor.
     * @param stringRepresentation the raw text in the file
     * @param fileContents the parsed contents of the file (outer list contains each line, inner list contains space-separated tokens)
     */
    public OutputFileFormat(String stringRepresentation, List<List<String>> fileContents) {
        this.stringRepresentation = stringRepresentation;
        this.fileContents = fileContents;
    }

    /**
     * Constructs an instance by reading text from filename.
     * Fields are parsed by using spaces as delimeters.  Consecutive delimiters are ignored.
     * @param filename the file to read from
     */
    public OutputFileFormat(String filename) {
        // get file length
        File file = new File(filename);
        long length = file.length(); // in bytes
        int characters = (int)(length/2L);

        List<List<String>> tempList = new LinkedList<>();
        StringBuilder builder = new StringBuilder(characters);
        String line = null;
        try (BufferedReader reader = new BufferedReader(new FileReader(filename))) {
            while ( (line = reader.readLine()) != null )
                {
                    line               = line.trim();
                    String[] fields    = line.split("\\s+");
                    tempList.add(ImmutableList.copyOf(fields));
                    builder.append(line);
                    builder.append("\n");
                }
        }
        catch (Exception e) {
            throw new IllegalArgumentException(e.getMessage());
        }
        stringRepresentation = builder.toString();
        fileContents = ImmutableList.copyOf(tempList);
    }

    /**
     * Return the contents of the file.  Newlines will be present.
     * @return the text that was in the file
     */
    @Override
    public String toString()  {
        return stringRepresentation;
    }

    @Override
    public int hashCode() {
        return Objects.hash(stringRepresentation, fileContents);
    }

    @Override
    public boolean equals(Object obj) {
        if ( obj == null )
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof OutputFileFormat) )
            return false;

        OutputFileFormat o = (OutputFileFormat)obj;
        if ( stringRepresentation.equals(o.stringRepresentation) &&
             fileContents.equals(o.fileContents) )
            return true;
        return false;
    }
}
