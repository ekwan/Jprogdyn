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
        List<String> EXPECTED_KEYS = ImmutableList.of("working_directory","frequency_directory","frequency_file",
													  "gaussian_directory","gaussian_max_filenames","number_of_simultaneous_trajectories",
													  "number_of_processors_per_trajectory","memory_per_trajectory","gaussian_force_route_card",
													  "gaussian_force_footer","job_type","trajectory_type",
													  "number_of_total_trajectories","checkpoint_prefix","checkpoint_interval",
													  "temperature","maximum_number_of_initialization_attempts","harmonic_tolerance",
													  "scale_factor","vibrational_initialization_default","vibrational_initialization_override",
													  "rotational_initialization_type","termination_condition",
													  "nmr_point_interval","shieldings_file","gaussian_nmr_route_card",
													  "gaussian_nmr_footer","symmetry_groups");

        // read and parse configuration file
        File configFile = new File(CONFIG_FILENAME);
        Map<String,String> configStringsMap = new LinkedHashMap<>();
        try ( Scanner scanner = new Scanner(configFile); ) {
            while (scanner.hasNextLine()) {
                // read and tokenize
                String line = scanner.nextLine().trim();
                if ( line.startsWith("#") || line.length() == 0 )
                    continue;
                int colonIndex = line.indexOf(":");
                if ( colonIndex == -1 || line.length() == colonIndex-1 ) {
                    System.out.printf("Ignoring malformed configuration line:\n");
                    System.out.println(line);
                }

                // parse
                String key = line.substring(0, colonIndex).toLowerCase().trim();
                String value = line.substring(colonIndex+1, line.length()).trim();
                int hashtagIndex = value.indexOf("#");
                if ( hashtagIndex > -1 )
                    value = value.substring(0, hashtagIndex).trim();

                // store values
                if ( configStringsMap.containsKey(key) ) {
                    // this is a duplicate
                    if ( key.equals("termination_condition") ) {
                        String currentValue = configStringsMap.get(key);
                        String newValue = String.format("%s;%s", currentValue, value);
                        configStringsMap.put(key, newValue);
                    }
                    else {
                        System.out.printf("Aborting! Duplicate key found in configuration file: %s\n", key);
                        System.exit(1);
                    }
                }
                else {
                    // this is not a duplicate
                    if ( value.equals("@blank") )
                        configStringsMap.put(key,"");
                    else if ( key.equals("working_directory") && value.toLowerCase().equals("use_current") ) {
                        String workingDirectory = System.getProperty("user.dir");
                        configStringsMap.put(key, workingDirectory);
                    }
                    else
                        configStringsMap.put(key, value);
                }

			}		
        }
        catch (FileNotFoundException e) {
            System.out.printf("Error: configuration file (%s) does not exist!", CONFIG_FILENAME);
            System.exit(1);
        }
    
   		// store the configuration data 
        Set<String> currentKeys = configStringsMap.keySet();
        Set<String> expectedKeys = new HashSet<>(EXPECTED_KEYS);
        expectedKeys.removeAll(currentKeys);
        if ( expectedKeys.size() > 0 ) {
            System.out.println("String entries missing from configuration file:");
            for (String s : expectedKeys)
                System.out.println(s);
            System.exit(1);
        }
        CONFIG_STRINGS_MAP = ImmutableMap.copyOf(configStringsMap);
        for (String key : CONFIG_STRINGS_MAP.keySet()) {
            String value = CONFIG_STRINGS_MAP.get(key);
            //System.out.printf("%s : %s\n", key, value);
        }

        // check invariants
        String workingDirectory = getString("working_directory");
        if ( ! new File(workingDirectory).isDirectory() )
            quit(String.format("working directory %d not found", workingDirectory));

        String frequencyDirectory = String.format("%s/%s", workingDirectory, getString("frequency_directory"));
        if ( ! new File(frequencyDirectory).isDirectory() )
            quit(String.format("frequency directory %s not found", frequencyDirectory));
        
        String frequencyFilename = String.format("%s/%s", frequencyDirectory, getString("frequency_file"));
        if ( ! new File(frequencyFilename).isFile() )
            quit(String.format("frequency file %s not found", frequencyFilename));

        String gaussianDirectory = String.format("%s/%s", workingDirectory, getString("gaussian_directory"));
        if ( ! new File(gaussianDirectory).isDirectory() )
            quit(String.format("gaussian directory %s not found", gaussianDirectory));

        String gaussianScriptFilename = String.format("%s/%s", gaussianDirectory, "run_gaussian.sh");
        if ( ! new File(gaussianScriptFilename).isFile() )
            quit(String.format("could not find %s", gaussianScriptFilename));

        int gaussianMaxFilenames = getInteger("gaussian_max_filenames");
        if ( gaussianMaxFilenames < 10 )
            quit("gaussian_max_filenames must be at least 10");

        int numberOfSimultaneousTrajectories = getInteger("number_of_simultaneous_trajectories");
        if ( numberOfSimultaneousTrajectories < 1 )
            quit("number of simultaneous trajectories must be at least 1");

        int numberOfProcessorsPerTrajectory = getInteger("number_of_processors_per_trajectory");
        if ( numberOfProcessorsPerTrajectory < 1 )
            quit("number of processors per trajectory must be at least 1");
        int memoryPerTrajectory = getInteger("memory_per_trajectory");
        
        if ( memoryPerTrajectory < 1 )
            quit("memory per trajectory must be at least 1");

        String gaussianForceRouteCard = getString("gaussian_force_route_card");
        if ( gaussianForceRouteCard.length() == 0 )
            quit("gaussian force route card cannot be blank");

        String jobType = getString("job_type");
        if ( ! ( jobType.equals("trajectory") || jobType.equals("analysis") ) )
            quit("job_type must be 'trajectory' or 'analysis'");

        String trajectoryType = getString("trajectory_type");
        if ( ! ( trajectoryType.equals("reaction") || trajectoryType.equals("nmr") ) )
            quit("trajectory_type must be 'reaction' or 'nmr'");

        int numberOfTotalTrajectories = getInteger("number_of_total_trajectories");
        if ( numberOfTotalTrajectories < 0 )
            quit("number of total trajectories must be non-negative");

        String checkpointPrefix = getString("checkpoint_prefix");
        if ( checkpointPrefix.length() == 0 )
            quit("blank checkpoint_prefix not allowed");

        int checkpointInterval = getInteger("checkpoint_interval");
        if ( checkpointInterval < 1 || checkpointInterval > 50 )
            quit("checkpoint interval must be betweeen [1,50]");

        double temperature = getDouble("temperature");
        if ( temperature < 0.0 )
            quit("temperature must be non-negative");
        
        int maximumNumberOfInitializationAttempts = getInteger("maximum_number_of_initialization_attempts");
        if ( maximumNumberOfInitializationAttempts < 10 )
            quit("maximum number of initialization attempts must be at least 10");
        
        double harmonicTolerance = getDouble("harmonic_tolerance");
        if ( harmonicTolerance < 0.0001 || harmonicTolerance > 10.0 )
            quit("harmonic tolerance must be [0.0001,10.0]");
        
        double scaleFactor = getDouble("scale_factor");
        if ( scaleFactor < 0.5 || scaleFactor > 1.5 )
            quit("scale_factor must be [0.5, 1.5]");
        
        Set<String> allowedInitializationTypes = new HashSet<>();
        for ( Initializer.VibrationalInitializationType t : Initializer.VibrationalInitializationType.values() )
            allowedInitializationTypes.add(t.toString().toLowerCase());
        String vibrationalInitializationDefault = getString("vibrational_initialization_default").toLowerCase();
        if ( ! allowedInitializationTypes.contains(vibrationalInitializationDefault) )
            quit(String.format("unrecognized value for vibrational_initialization_default: %s"));
        
        String vibrationalInitializationOverride = getString("vibrational_initialization_override");
        if ( ! vibrationalInitializationOverride.equals("no_overrides") ) {
            try {
                String[] overrides = vibrationalInitializationOverride.split(";");
                for (String overrideString : overrides) {
                    String[] fields = overrideString.split(":");
                    if ( fields.length != 2 )
                        quit(String.format("invalid vibrational initialization override string: %s", overrideString));
                    int modeNumber = Integer.parseInt(fields[0]);
                    if ( modeNumber < 0 )
                        quit(String.format("negative mode number not allowed"));
                    if ( ! allowedInitializationTypes.contains( fields[1].trim().toLowerCase() ) )
                        quit(String.format("check initialization type for override: %s", fields[1]));
                }
            }
            catch (Exception e) {
                quit("check vibrational initialization string");
            }
        }

        allowedInitializationTypes = new HashSet<>();
        for ( Initializer.RotationalInitializationType t : Initializer.RotationalInitializationType.values() )
            allowedInitializationTypes.add(t.toString().toLowerCase());
        String rotationalInitializationType = getString("rotational_initialization_type").toLowerCase();
        if ( ! allowedInitializationTypes.contains(rotationalInitializationType) )
            quit(String.format("unrecognized value for rotational_initialization_type: %s", rotationalInitializationType));

        String terminationConditions = getString("termination_condition").toLowerCase();
        if ( ! terminationConditions.equals("no_termination_conditions") ) {
            String terminationConditionsSplit[] = terminationConditions.split(";");
            for (String terminationCondition : terminationConditionsSplit) {
                try {
                    String[] fields = terminationCondition.split(",");
                    String coordinateType = fields[0].toLowerCase();
                    
                    if ( ! ( ( coordinateType.equals("bond_length") && fields.length == 6 ) ||
                             ( coordinateType.equals("bond_angle") && fields.length == 7 ) ||
                             ( coordinateType.equals("torsion") && fields.length == 8 )        ) )
                        throw new Exception();

                    String comparisonOperatorString = fields[fields.length-2].trim().toLowerCase();
                    if ( ! ImmutableSet.of("greater_than","equals","less_than").contains( comparisonOperatorString ) )
                        throw new Exception();

                }
                catch (Exception e) {
                    quit(String.format("check termination conditions string: %s", terminationCondition));
                }
            }
        }

		if ( jobType.equals("nmr") ) {
			int nmrPointInterval = getInteger("nmr_point_interval");
			if ( nmrPointInterval < 1 )
				quit("check nmr_point_interval");
			
			String shieldingsFilename = String.format("%s/%s", frequencyDirectory, getString("shieldings_file"));
			if ( ! new File(shieldingsFilename).isFile() )
				quit(String.format("nmr shieldings file %s not found", shieldingsFilename));

            String gaussianNmrRouteCard = getString("gaussian_nmr_route_card");
            if ( gaussianNmrRouteCard.length() == 0 )
                quit("can't have blank nmr route card");

            try {
                String symmetryGroupsString = getString("symmetry_groups");
                String[] symmetryGroups = symmetryGroupsString.split(";");
                Set<Integer> alreadySeen = new HashSet<>();
                for (String symmetryGroup : symmetryGroups) {
                    String fields[] = symmetryGroup.split(",");
                    for (String s : fields) {
                        int atomNumber = Integer.parseInt(s.trim());
                        if ( atomNumber < 0 )
                            throw new Exception();
                        if ( alreadySeen.contains(atomNumber) )
                            throw new Exception();
                        alreadySeen.add(atomNumber);
                    }
                }
            }
            catch (Exception e) {
                quit("check symmetry groups string");
            }
		}

    }

    /**
     * Quits with an error message.
     * @param errorMessage the error message
     */
    public static void quit(String errorMessage) {
        System.out.println("Aborting due to error:");
        System.out.println(errorMessage);
        System.exit(1);
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
        System.out.println("Loaded configuration data successfully.");
        System.out.printf("Jprogdyn is running on host %s.\n", HOSTNAME);
    }
}
