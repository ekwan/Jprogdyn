package edu.harvard.chemistry.ekwan.Jprogdyn;

import java.util.*;
import java.io.*;

/**
 * This class represents the output of a Gaussian job.
 */
public class GaussianResult implements Result {

    /** The result of the calculation. */
    public final GaussianOutputFile out; 

    /** How long the calculation took in wallclock time (seconds). */
    public final double elapsedTime;

    /**
     * Reads the result of the job from the specified filename.
     * @param filename the filename of the Gaussian output file to read from
     * @param elapsedTime how long the calculation took in wallclock time (seconds)
     */
    public GaussianResult(String filename, double elapsedTime) {
        if ( ! new File(filename).exists() )
            throw new IllegalArgumentException(String.format("Filename %s not found!", filename));
        OutputFileFormat temp = new OutputFileFormat(filename) {};
        List<List<String>> fileContents = temp.fileContents;
        if ( fileContents.size() == 0 )
            throw new IllegalArgumentException("g09 output file is empty");

        String debugString = "";
        for (int i = Math.max(0, temp.fileContents.size()-5); i < temp.fileContents.size(); i++) {
            List<String> line = temp.fileContents.get(i);
            for (String s : line)
                debugString += s + " ";
            debugString += "\n";
        }
        debugString = debugString.substring(0, debugString.length()-1);

        List<String> lastLine = fileContents.get(fileContents.size()-1);
        if ( lastLine.size() < 2 )
            throw new IllegalArgumentException("Gaussian output last line truncated, end of file follows:\n" + debugString);
        if ( ! (lastLine.get(0).equals("Normal") && (lastLine.get(1).equals("termination")) ) )
            throw new IllegalArgumentException("Gaussian job did not terminate normally, end of file follows:\n" + debugString);
        out = new GaussianOutputFile(filename);
        this.elapsedTime = elapsedTime;
    }

    @Override
    public int hashCode() {
        return Objects.hash(out);
    }

    @Override
    public boolean equals(Object obj) {
        if ( obj == null )
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof GaussianResult) )
            return false;

        GaussianResult r = (GaussianResult)obj;
        if ( out.equals(r.out) )
            return true;
        return false;
    }

    @Override
    public String toString() {
        return "GaussianResult";
    }
}
