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
													  "number_of_total_trajectories","checkpoint_directory", "checkpoint_prefix","checkpoint_interval",
													  "temperature","timestep","number_of_forward_points","number_of_backward_points", 
                                                      "maximum_number_of_initialization_attempts","harmonic_tolerance",
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

                // check for unexpected keys
                if ( ! EXPECTED_KEYS.contains(key) )
                    quit(String.format("unexpected line in configuration file:\n%s", line));

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

        String gaussianForceFooter = getString("gaussian_force_footer");

        String jobType = getString("job_type");
        if ( ! ( jobType.equals("trajectory") || jobType.equals("analysis") ) )
            quit("job_type must be 'trajectory' or 'analysis'");

        String trajectoryType = getString("trajectory_type");
        if ( ! ( trajectoryType.equals("reaction") || trajectoryType.equals("nmr") ) )
            quit("trajectory_type must be 'reaction' or 'nmr'");

        int numberOfTotalTrajectories = getInteger("number_of_total_trajectories");
        if ( numberOfTotalTrajectories < 0 )
            quit("number of total trajectories must be non-negative");

        String checkpointDirectory = String.format("%s/%s", workingDirectory, getString("checkpoint_directory"));
        if ( ! new File(checkpointDirectory).isDirectory() )
            quit(String.format("can't find checkpoint directory %s", checkpointDirectory));

        String checkpointPrefix = getString("checkpoint_prefix");
        if ( checkpointPrefix.length() == 0 )
            quit("blank checkpoint_prefix not allowed");
    
        int checkpointInterval = getInteger("checkpoint_interval");
        if ( checkpointInterval < 1 || checkpointInterval > 50 )
            quit("checkpoint interval must be betweeen [1,50]");

        double temperature = getDouble("temperature");
        if ( temperature < 0.0 )
            quit("temperature must be non-negative");
        
        double timestep = getDouble("timestep");
        if ( timestep < 0.0001 || timestep > 10 )
            quit("timestep out of range [0.0001,10]");

        int numberOfForwardPoints = getInteger("number_of_forward_points");
        if ( numberOfForwardPoints < 0 )
            quit("number of forward trajectory points must be non-negative");

        int numberOfBackwardPoints = getInteger("number_of_backward_points");
        if ( numberOfBackwardPoints < 0 )
            quit("number of backward trajectory points must be non-negative");

        if ( numberOfForwardPoints + numberOfBackwardPoints < 1 )
            quit("must calculate some trajectory points");

        int maximumNumberOfInitializationAttempts = getInteger("maximum_number_of_initialization_attempts");
        if ( maximumNumberOfInitializationAttempts < 10 )
            quit("maximum number of initialization attempts must be at least 10");
        
        double harmonicTolerance = getDouble("harmonic_tolerance");
        if ( harmonicTolerance < 0.0001 || harmonicTolerance > 10.0 )
            quit("harmonic tolerance must be [0.0001,10.0]");
        
        double scaleFactor = getDouble("scale_factor");
        if ( scaleFactor < 0.5 || scaleFactor > 1.5 )
            quit("scale_factor must be [0.5, 1.5]");
        
        String vibrationalInitializationDefaultString = getString("vibrational_initialization_default").toUpperCase();
        Initializer.VibrationalInitializationType vibrationalInitializationDefault = null;
        try {
            vibrationalInitializationDefault = Initializer.VibrationalInitializationType.valueOf(vibrationalInitializationDefaultString);
        }
        catch (Exception e) {
            quit("check default vibrational initialization type");
        }

        String vibrationalInitializationOverrideString = getString("vibrational_initialization_override").toUpperCase();
        Map<Integer,Initializer.VibrationalInitializationType> specialModeInitializationMap = new HashMap<>();
        if ( ! vibrationalInitializationOverrideString.equals("NO_OVERRIDES") ) {
            String[] overrides = vibrationalInitializationOverrideString.split(";");
            for (String overrideString : overrides) {
                String[] fields = overrideString.split(":");
                if ( fields.length != 2 )
                    quit(String.format("invalid vibrational initialization override string: %s", overrideString));
                int modeNumber = Integer.parseInt(fields[0]);
                if ( modeNumber < 0 )
                    quit(String.format("negative mode number not allowed"));
                if ( specialModeInitializationMap.containsKey(modeNumber) )
                    quit(String.format("duplicate mode number given:\n%s", overrideString));
                
                Initializer.VibrationalInitializationType type = null;
                try {    
                    type = Initializer.VibrationalInitializationType.valueOf(fields[1].trim());
                }
                catch (Exception e) {
                    quit(String.format("check vibrational initialization string\n%s", overrideString));
                }
                specialModeInitializationMap.put(modeNumber, type);
            }
        }

        String rotationalInitializationTypeString = getString("rotational_initialization_type").toUpperCase();
        Initializer.RotationalInitializationType rotationalInitializationType = null;
        try {
            rotationalInitializationType = Initializer.RotationalInitializationType.valueOf(rotationalInitializationTypeString);
        }
        catch (Exception e) {
            quit("check rotational initialization type");
        }
        
        String terminationConditionsString = getString("termination_condition");
        List<InternalCoordinate.Condition> terminationConditions = new LinkedList<>();
        
        if ( ! terminationConditionsString.equals("no_termination_conditions") && jobType.equals("trajectory") ) {
            String terminationConditionsSplit[] = terminationConditionsString.split(";");
            for (String terminationCondition : terminationConditionsSplit) {
                try {
                    String[] fields = terminationCondition.split(",");
                    for (int i=0; i < fields.length; i++)
                        fields[i] = fields[i].trim();
                    
                    String comparisonOperatorString = fields[fields.length-2].toUpperCase();
                    InternalCoordinate.ConditionType conditionType = InternalCoordinate.ConditionType.valueOf(comparisonOperatorString);
                    
                    double value = Double.parseDouble(fields[fields.length-1]);
                    String description = fields[fields.length-3];
                    
                    String coordinateType = fields[0].toLowerCase();
                    InternalCoordinate.Condition condition = null;
                    if ( coordinateType.equals("bond_length") && fields.length == 6 ) {
                        InternalCoordinate.Length bond = new InternalCoordinate.Length(Integer.parseInt(fields[1])-1, Integer.parseInt(fields[2])-1, description);
                        condition = new InternalCoordinate.Condition(bond, conditionType, value);
                    }
                    else if ( coordinateType.equals("bond_angle") && fields.length == 7 ) {
                        InternalCoordinate.Angle angle = new InternalCoordinate.Angle(Integer.parseInt(fields[1])-1, Integer.parseInt(fields[2])-1,
                                                                                      Integer.parseInt(fields[3])-1, description);
                        condition = new InternalCoordinate.Condition(angle, conditionType, value);
                    }
                    else if ( coordinateType.equals("torsion") && fields.length == 8 ) {
                        InternalCoordinate.Torsion torsion = new InternalCoordinate.Torsion(Integer.parseInt(fields[1])-1, Integer.parseInt(fields[2])-1,
                                                                                            Integer.parseInt(fields[3])-1, Integer.parseInt(fields[4])-1, description);
                        new InternalCoordinate.Condition(torsion, conditionType, value);
                    }
                    else
                        throw new Exception();
                    terminationConditions.add(condition);
                }
                catch (Exception e) {
                    quit(String.format("check termination conditions string: %s", terminationCondition));
                }
            }
        }

		int nmrPointInterval = 0;
		String shieldingsFilename = null;
	    String gaussianNmrRouteCard = null;
		if ( trajectoryType.equals("nmr") ) {
			nmrPointInterval = getInteger("nmr_point_interval");
			if ( nmrPointInterval < 1 )
				quit("check nmr_point_interval");
			
			shieldingsFilename = String.format("%s/%s", frequencyDirectory, getString("shieldings_file"));
			if ( ! new File(shieldingsFilename).isFile() )
				quit(String.format("nmr shieldings file %s not found", shieldingsFilename));

            gaussianNmrRouteCard = getString("gaussian_nmr_route_card");
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

        // ready to go
        System.out.println("Loaded configuration data successfully.");
        System.out.printf("Jprogdyn is running on host %s.\n", HOSTNAME);

        // determine what kind of job is being requested
        if ( jobType.equals("trajectory") && trajectoryType.equals("reaction") ) {
            // run reaction trajectories
            System.out.println("Will run reaction trajectories.");
        }
        else if ( jobType.equals("analysis") && trajectoryType.equals("reaction") ) {
            // run reaction trajectory analysis
            System.out.println("Will run reaction trajectory analysis.");
        }
        else if ( jobType.equals("trajectory") && trajectoryType.equals("nmr") ) {
            // run NMR trajectories
            System.out.println("Will run NMR trajectories.");
        
            // read frequency molecule
            GaussianOutputFile frequenciesOutputFile = new GaussianOutputFile(frequencyFilename);
            Molecule frequenciesMolecule = frequenciesOutputFile.molecule;
            System.out.printf("Read frequency data from %s (%d atoms, %d normal modes).\n", frequencyFilename,
                              frequenciesMolecule.contents.size(), frequenciesMolecule.modes.size());

            // read shieldings molecule
			GaussianOutputFile shieldingsOutputFile = new GaussianOutputFile(shieldingsFilename);
			Molecule shieldingsMolecule = shieldingsOutputFile.molecule;
            System.out.printf("Read chemical shift data from %s.\n", shieldingsFilename);
            if ( frequenciesMolecule.contents.size() != shieldingsMolecule.contents.size() )
                quit("mismatch between frequency and shieldings molecules");

			// make calculation methods
			String gaussianForceRouteCardFull = String.format("#p force %s", gaussianForceRouteCard);
            String gaussianForceFooterFull = String.format("%s\n", gaussianForceFooter);
            CalculationMethod dynamicsMethod = new GaussianCalculationMethod(CalculationMethod.CalculationType.ENERGY_AND_FORCE,
																			 memoryPerTrajectory, numberOfProcessorsPerTrajectory,
																			 gaussianForceRouteCardFull, gaussianForceFooterFull);
            String gaussianNMRFooterFull = String.format("--Link1--\n%%chk=Jprogdyn.chk\n%%nprocshared=%d\n%%mem=%dGB\n#p NMR geom=allcheck guess=read %s",
                                                         numberOfProcessorsPerTrajectory, memoryPerTrajectory, gaussianForceRouteCard);
			GaussianCalculationMethod nmrMethod = new GaussianCalculationMethod(CalculationMethod.CalculationType.NMR,
																		        memoryPerTrajectory, numberOfProcessorsPerTrajectory,
																		        gaussianForceRouteCardFull, gaussianNMRFooterFull);

            // create the requested trajectories
            System.out.println("Generating trajectories...");
			List<Trajectory> trajectories = new ArrayList<>(numberOfTotalTrajectories);
			for (int i=0; i < numberOfTotalTrajectories; i++)
				{
                    String checkpointFilename = String.format("%s/%s/%s_%04d.chk", workingDirectory, checkpointDirectory, checkpointPrefix, i);
                    Initializer initializer = new Initializer(frequenciesMolecule, temperature, timestep, vibrationalInitializationDefault,
                                                              rotationalInitializationType,
                                                              specialModeInitializationMap,
                                                              harmonicTolerance, dynamicsMethod, scaleFactor);
                    Trajectory trajectory = new Trajectory(frequenciesMolecule, timestep, numberOfForwardPoints, numberOfBackwardPoints,
                                                           new ArrayList<InternalCoordinate.Condition>(), maximumNumberOfInitializationAttempts,
                                                           initializer, dynamicsMethod, nmrMethod, nmrPointInterval, checkpointFilename); 
                    trajectories.add(trajectory);
                    trajectory.call();
                    System.out.printf("   Generated %s.\n", checkpointFilename);
				}

            // run the trajectories
            //TrajectoryExecutorService.runTrajectories(trajectories);
        }
        else if ( jobType.equals("analysis") && trajectoryType.equals("nmr") ) {
            // run NMR trajectory analysis
            System.out.println("Will run NMR trajectory analysis.");
        }
        else
            quit("unrecognized job type");
    
        // finished
        System.out.println("Jprogdyn has terminated successfully.");
    }
}
