package edu.harvard.chemistry.ekwan.Jprogdyn;

import java.io.*;
import java.util.*;

/**
 * This class starts Jprogdyn.  It reads the specified configuration file,
 * loads any saved trajectories, and launches trajectories.
 */
public class Loader {
    /** Location of the configuration file. */
    public static final String CONFIG_FILENAME;

    static {
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
    }

    public static void main(String[] args) {
        System.out.println("hello");
    }
}
