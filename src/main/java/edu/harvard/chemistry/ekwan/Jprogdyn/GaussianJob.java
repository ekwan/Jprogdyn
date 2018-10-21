package edu.harvard.chemistry.ekwan.Jprogdyn;

import java.util.*;
import java.io.*;
import java.util.concurrent.atomic.*;
import org.apache.commons.io.FileUtils;
import java.util.concurrent.*;
import com.google.common.collect.*;

/**
 * This class represents a Gaussian job.  The Gaussian executable must be in the path for this to work.
 */
public class GaussianJob implements Callable<GaussianResult>, Serializable {
    
    /** For serialization. */
    public static final long serialVersionUID = 1L;

    /** Creates temporary filename indices to run the jobs with. */
    public static final AtomicInteger INDEX_GENERATOR = new AtomicInteger();

    /** The input file to run. */
    public final GaussianInputFile gjf;

    /** The maximum number of filenames in one Gaussian directory. */
    public static final int MAX_FILENAMES = Loader.getInteger("gaussian_max_filenames");
    
    /** Where to run jobs. */
    public static final String GAUSSIAN_DIRECTORY = String.format("%s/%s", Loader.getString("working_directory"), Loader.getString("gaussian_directory"));

    /**
     * Constructor.
     * @param gjf the name of the Gaussian input file to run
     */
    public GaussianJob(GaussianInputFile gjf) {
        this.gjf = gjf;
    }

    /**
     * Auto-selects a filename and runs the analysis calculation.
     * @return the result of the calculation
     */
    public GaussianResult call() {
        // choose a base filename for this set of jobs
        String baseFilename = "";
        int index = 0;
        
        // expand environment variable if necessary
        
        counting:
        for (int i=0; i < MAX_FILENAMES; i++) {
            // get a new ID number for this job
            index = INDEX_GENERATOR.getAndIncrement();
            baseFilename = String.format("%s_gaussian_%010d", Loader.HOSTNAME, index);

            // don't allow this choice of filenames if any files with this prefix already exist
            for ( File f : new File(GAUSSIAN_DIRECTORY).listFiles() ) {
                if ( f.getName().startsWith(baseFilename) ) {
                    baseFilename = "";
                    continue counting;
                }
            }

            // reset counter if necessary
            if ( INDEX_GENERATOR.get() > MAX_FILENAMES )
                INDEX_GENERATOR.getAndSet(0);
            break;
        }
        if ( baseFilename.length() == 0 )
            throw new IllegalArgumentException("Unable to set filename!");

        // write input files to disk
        String jobDirectoryName = String.format("%s/%s", GAUSSIAN_DIRECTORY, baseFilename);
        File jobDirectory = new File(jobDirectoryName);
        boolean success = jobDirectory.mkdir();
        if ( !success )
            throw new IllegalArgumentException("failed to create directory " + jobDirectoryName);
        String gjfFilename = String.format("%s/gaussian.gjf", jobDirectoryName);
        gjf.write(gjfFilename);

        // call Gaussian
        double elapsedTime = 0.0;
        try {
                long startTime = System.currentTimeMillis();
                
                String runString = String.format("%s/run_gaussian.sh %s %s", GAUSSIAN_DIRECTORY, jobDirectoryName, baseFilename);
                //System.out.println(runString);
                Process process = Runtime.getRuntime().exec(runString);
                int exitValue = process.waitFor();
                long endTime = System.currentTimeMillis();
                elapsedTime = (endTime - startTime) / 1000.0;
        }
        catch (Exception e) {
            System.out.println("Error while running Gaussian job:");
            e.printStackTrace();
        }

        // retrieve output
        String outputFilename = String.format("%s/%s/gaussian.out", GAUSSIAN_DIRECTORY, baseFilename);
        GaussianResult result = new GaussianResult(outputFilename, elapsedTime);
	
        // remove files
        try {
            FileUtils.deleteDirectory(jobDirectory);
        }
        catch (Exception e) {
            System.out.println("Error while trying to delete directory: " + jobDirectory.getName());
            e.printStackTrace();
        }
        
        // return result
        return result;
    }

    @Override
    public int hashCode() {
        return Objects.hash(gjf);
    }

    @Override
    public boolean equals(Object obj) {
        if ( obj == null )
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof GaussianJob) )
            return false;

        GaussianJob j = (GaussianJob)obj;
        if ( gjf.equals(j.gjf) )
            return true;
        return false;
    }

    @Override
    public String toString() {
        return "GaussianJob";
    }
}
