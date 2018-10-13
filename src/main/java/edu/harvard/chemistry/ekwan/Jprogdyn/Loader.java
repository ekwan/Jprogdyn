package edu.harvard.chemistry.ekwan.Jprogdyn;

import java.io.*;
import java.util.*;
import com.google.common.collect.*;

/**
 * This class starts Jprogdyn.  It reads the specified configuration file,
 * loads any saved trajectories, and launches trajectories.
 */
public class Loader {
    
    /** Location of the configuration file. */
    public static final String CONFIG_FILENAME;

    /** The configuration data. */
    public static final Map<String,String> CONFIG_STRINGS_MAP;

    /** The name of the current host. */
    public static final String HOSTNAME;

	/**
     * This static initializer reads the specified configuration file.
	 */
    static {
        // set hostname
        String temp = "localhost";
        try {
                temp = java.net.InetAddress.getLocalHost().getHostName();
        }
        catch (Exception e) {
        }

        if ( temp.length() > 0 )
            temp = temp.split("\\.")[0];
        else
            temp = "localhost";
        HOSTNAME = temp;

        // check command line arguments
        System.out.println("\n=== Jprogdyn 1.0 ===");
        String configFilename = "Jprogdyn.config";
        String customConfigFilename = System.getProperty("config.filename");
        if ( customConfigFilename != null && customConfigFilename.trim().length() > 0 )
            configFilename = customConfigFilename;
        if ( ! new File(configFilename).isFile() ) {
            System.out.printf("Error: couldn't find configuration file: %s.\n", configFilename);
            System.exit(1);
        }
        CONFIG_FILENAME = configFilename;
        System.out.printf("Will read configuration data from %s.\n", CONFIG_FILENAME);
    
        // the configuration file must contain exactly these fields
        List<String> EXPECTED_KEYS = ImmutableList.of("bam_directory");

        // these fields contain strings that reference enums
        List<String> ENUM_KEYS = ImmutableList.of("ignore_duplicate_reads");
        for (String s : ENUM_KEYS) {
            if ( ! EXPECTED_KEYS.contains(s) )
                throw new IllegalArgumentException("check enum key " + s);
        }


        CONFIG_STRINGS_MAP = null;
    }

    /**
     * Convenience method for retrieving a string key.
     * Throws a NullPointerException for invalid keys.
     *
     * @param keyName the name of the key
     * @return the corresponding value
     */
    public static String getString(String keyName) {
        if ( CONFIG_STRINGS_MAP.containsKey(keyName) )
            return CONFIG_STRINGS_MAP.get(keyName);
        else
            throw new NullPointerException("no value found for key: " + keyName);
    }

    /**
     * Convenience method for retrieving an integer key.
     * Throws a NullPointerException for invalid keys.
     *
     * @param keyName the name of the key
     * @return the corresponding value
     */
    public static int getInteger(String keyName) {
        return Integer.parseInt(getString(keyName));
    }

    /**
     * Convenience method for retrieving a floating point key.
     * Throws a NullPointerException for invalid keys.
     *
     * @param keyName the name of the key
     * @return the corresponding value
     */
    public static double getDouble(String keyName) {
        return Double.parseDouble(getString(keyName));
    }

    /**
     * Executes the desired trajectories given the data from the configuration file.
     * @param args command-line arguments other than the configuration filename are ignored
     */
    public static void main(String[] args) {
        System.out.println("hello");
    }
}
