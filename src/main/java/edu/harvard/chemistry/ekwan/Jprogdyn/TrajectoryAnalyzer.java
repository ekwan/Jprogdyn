package edu.harvard.chemistry.ekwan.Jprogdyn;

import java.util.*;
import java.io.*;
import com.google.common.collect.*;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.geometry.euclidean.threed.*;

/** 
 * This class combines some static methods for analyzing trajectories.
 */
public class TrajectoryAnalyzer implements Immutable, Singleton
{
    /** Not instantiable. */
    private TrajectoryAnalyzer()
    {
        throw new IllegalArgumentException("not instantiable");
    }

    /**
     * Takes a bunch of Trajectories in the form of serialized .chk files and makes MOLDEN movies.
     * Existing files will be overwritten.
     * @param trajectories the trajectories to analyze
     * @param directory where to output the movies to
     */
    public static void makeMovies(Map<String,Trajectory> trajectories, String directory)
    {
        int count = 0;
        for (String filename : trajectories.keySet())
            {
                Trajectory t = trajectories.get(filename);
                String trajFilename = directory + "/" + filename.replace(".chk", ".traj");
                t.writeTrajString(trajFilename);
                count++;
            }
        System.out.printf("%d movies written.\n", count);
    }

    /**
     * Convenience method for loading a bunch of related trajectories at once.
     * @param directory the directory in which to search for .chk files
     * @param prefix only take .chk files that start with this prefix (case-sensitive)
     * @return filenames mapped to trajectories
     */
    public static Map<String,Trajectory> loadTrajectories(String directory, String prefix)
    {
        int globalPoints = 0;
        Map<String,Trajectory> returnMap = new LinkedHashMap<>();
        File[] files = new File(directory).listFiles();
        Arrays.sort(files);
        for (File f : files)
            {
                String filename = f.getName();
                if ( filename.startsWith(prefix) && filename.endsWith("chk") )
                    {
                        Trajectory t = Trajectory.loadCheckpoint(directory + "/" + filename);
                        if ( t != null )
                            {
                                int totalPoints = t.backwardPoints.size() + t.forwardPoints.size();
                                if ( t.initialPoint != null )
                                    totalPoints++;
                                if ( totalPoints < 50 )
                                    {
                                        System.out.printf("Trajectory only has %d points, skipping (%s).\n", totalPoints, filename);
                                        continue;
                                    }
                                globalPoints += totalPoints;
                                System.out.printf("Loaded %d trajectory points from %s.\n", totalPoints, filename);
                                returnMap.put(filename, t);
                            }
                    }
            }
        if ( returnMap.size() == 0 )
            System.out.println("No trajectories loaded!");
        else
            System.out.printf("%d total points loaded from %d trajectories.\n", globalPoints, returnMap.size());
        return returnMap;
    }

    /** 
     * Prints out a report on the total energy stability of a bunch of trajectories.
     * Note that initial points are not assigned a total energy and are ignored.
     * @param trajectories the trajectories to analyze
     */
    public static void analyzeStability(Map<String,Trajectory> trajectories)
    {
        System.out.println("Trajectory                                         points      Std. Dev. (%)       Std. Dev. (kcal)");
        System.out.println("                                                                               0.0                  5.0"); 
        for (String filename : trajectories.keySet())
            {
                Trajectory t = trajectories.get(filename);
                DescriptiveStatistics stats = new DescriptiveStatistics();
                for (TrajectoryPoint p : t.backwardPoints)
                    stats.addValue(p.totalEnergy);
                for (TrajectoryPoint p : t.forwardPoints)
                    stats.addValue(p.totalEnergy);
                double mean = stats.getMean();
                double standardDeviation = stats.getStandardDeviation();
                double COV = Math.abs(standardDeviation*100.0/mean);
                standardDeviation = standardDeviation * 627.509469;
                long n = stats.getN();
                String asciiBar = AsciiBar.make(standardDeviation, 0.0, 5.0, "*", " ", 20);
                System.out.printf("%-40s           %4d        %9.6f%%      %s%-7.1f\n", filename, n, COV, asciiBar, standardDeviation);
            }
    }

    /**
     * Prints out a report on some geometric parameters for a bunch of trajectories.
     * @param trajectories the trajectories to analyze
     * @param coordinates the coordinates to print out
     * @param pointInterval how often to print out the coordinates
     * @param references known starting materials, products, etc. 
     */
    public static void analyzeGeometry(Map<String,Trajectory> trajectories, List<InternalCoordinate> coordinates,
                                       int pointInterval,
                                       Map<List<InternalCoordinate.Condition>,String> references)
    {
        if ( coordinates == null || coordinates.size() == 0 )
            throw new IllegalArgumentException("coordinates are null or empty");

        // determine min and max for each coordinate
        List<Double> minima = new ArrayList<>(coordinates.size());
        List<Double> maxima = new ArrayList<>(coordinates.size());
        for (InternalCoordinate coordinate : coordinates) {
            if ( coordinate == null )
                throw new IllegalArgumentException("null coordinates are not allowed");
            Double thisMinimum = null;
            Double thisMaximum = null;
            for (String filename : trajectories.keySet()) {
                Trajectory t = trajectories.get(filename);
                List<TrajectoryPoint> points = t.getPoints();
                for (TrajectoryPoint p : points) {
                    double value = coordinate.getValue(p.positions);
                    if ( thisMinimum == null || value < thisMinimum )
                        thisMinimum = value;
                    if ( thisMaximum == null || value > thisMaximum )
                        thisMaximum = value;
                }
            }
            if ( thisMinimum == null || thisMaximum == null )
                throw new NullPointerException("need some points");
            //System.out.printf("%s  min %.2f  max %.2f\n", coordinate.toString(), thisMinimum, thisMaximum);
            minima.add(thisMinimum);
            maxima.add(thisMaximum);
        }

        // print out the geometries for each trajectory in tabular format
        Map<String,Integer> outcomeMap = new LinkedHashMap<>();  // map from descriptions of outcomes to how many times the outcomes have been observed
        for (String description : references.values())
            outcomeMap.put(description, 0);
        outcomeMap.put("unknown (incomplete)", 0);
        outcomeMap.put("unknown (complete)", 0);
        System.out.println("\n=== Geometry Analysis ===\n");
        for (String filename : trajectories.keySet()) {
            Trajectory t = trajectories.get(filename);
            System.out.printf("%s (%s)\n\nTime(fs)", filename, t.isDone() ? "complete" : "incomplete");
            for (InternalCoordinate coordinate : coordinates)
                System.out.printf("        %-15s         ", coordinate.toString());
            System.out.println("\n");
            List<TrajectoryPoint> points = t.getPoints();
            for (int i=0; i < points.size(); i += pointInterval) {
                TrajectoryPoint p = points.get(i);
                List<Vector3D> positions = p.positions;
                double time = p.time;
                System.out.printf("%5.0f     ", time);
                for (int j=0; j < coordinates.size(); j++)
                    {
                        InternalCoordinate coordinate = coordinates.get(j);
                        double value = coordinate.getValue(positions);
                        double min = minima.get(j);
                        double max = maxima.get(j);
                        System.out.print(AsciiBar.make(value, min, max, "*", " ", 20));
                        System.out.printf("%-7.2f ", value);
                    }
                System.out.println();
            }

            // guarantees the endpoint will be printed
            if ( points.size() % pointInterval != 0 ) {
                int i = points.size()-1;
                TrajectoryPoint p = points.get(i);
                List<Vector3D> positions = p.positions;
                double time = p.time;
                System.out.printf("%5.0f     ", time);
                for (int j=0; j < coordinates.size(); j++) {
                    InternalCoordinate coordinate = coordinates.get(j);
                    double value = coordinate.getValue(positions);
                    double min = minima.get(j);
                    double max = maxima.get(j);
                    System.out.print(AsciiBar.make(value, min, max, "*", " ", 20));
                    System.out.printf("%-7.2f ", value);
                }
                System.out.println();
             }

            // find the earliest forward point, if any, that satisfies all conditions 
            String description = null;
            Double time = null;
            points:
            for (TrajectoryPoint p : t.forwardPoints) {
                // get position of each point
                List<Vector3D> positions = p.positions;
                
                // check each set of conditions
                for (List<InternalCoordinate.Condition> conditions : references.keySet()) {
                    boolean reached = true;
                    // each sub-condition must be true
                    for (InternalCoordinate.Condition condition : conditions) {
                        if ( !condition.reached(positions) ) {
                            reached = false;
                            break;
                        }
                    }
                   
                    if ( reached ) {
                        description = references.get(conditions);
                        time = p.time;
                        Integer currentNumber = outcomeMap.get(description);
                        currentNumber++;
                        outcomeMap.put(description, currentNumber);

                        // this is the earliest point where the condition has been satisified
                        // stop looking, even if future points no longer satisfy the conditions
                        //
                        // breaking out here prevents conditions from being double-counted,
                        // although what is defined as the "first" condition to be met is arbitrary
                        break points;
                    }
                }
            }

            // report results
            if ( time != null )
                // time is not null, so we were able to assign a forward point to a specified outcome
                System.out.printf("%s reached after %.0f fs.\n", description, time);
            else {
                // time is null, so we weren't able to assign any of the forward points to a specified outcome
                System.out.printf("No known species reached after %.0f fs.\n", time);
                
                // only note this if the trajectory is done
                if ( t.isDone() ) {
                    Integer currentNumber = outcomeMap.get("unknown (complete)");
                    currentNumber++;
                    outcomeMap.put("unknown (complete)", currentNumber++);
                }
                else {
                    Integer currentNumber = outcomeMap.get("unknown (incomplete)");
                    currentNumber++;
                    outcomeMap.put("unknown (incomplete)", currentNumber++);
                }
            }
            System.out.println();
        }

        // print out summary of trajectory outcomes
        System.out.println("\n=== Summary of Trajectories ===\n");
        System.out.println("Species                       Trajs          %");
        int total = 0;
        for (String description : outcomeMap.keySet()) {
            Integer currentNumber = outcomeMap.get(description);
            System.out.printf("%-30s  %3d     %6.0f\n", description, currentNumber, currentNumber*100.0/trajectories.size());
            total += currentNumber;
        }
        System.out.printf("\n%d trajectories total.\n", trajectories.size());
        if ( total != trajectories.size() )
            System.out.println("Warning, some double counting occurred because some trajectories matched more than one species.");
    }

    /**
     * Writes a text file of the internal coordinates as a function of time for all trajectories.
     * @param trajectoriesMap map from trajectory to CSV output filename
     * @param coordinates the coordinates to print out
     */
    
    public static void writeScatter(Map<Trajectory,String> trajectoriesMap, List<InternalCoordinate> coordinates)
    {
        if ( trajectoriesMap == null || trajectoriesMap.size() == 0 )
            throw new IllegalArgumentException("null or empty trajectories");
        if ( coordinates == null || coordinates.size() == 0 )
            throw new IllegalArgumentException("null or empty coordinates");

        for (Trajectory t : trajectoriesMap.keySet()) {
            String CSVfilename = trajectoriesMap.get(t);
            StringBuilder sb = new StringBuilder();
            List<TrajectoryPoint> points = t.getPoints();
            String header = "time";
            for (InternalCoordinate c : coordinates)
                header += "," + c.description;
            header += "\n";
            sb.append(header);
            for (TrajectoryPoint p : points) {
                sb.append(String.format("%.2f", p.time));
                for (InternalCoordinate c : coordinates)
                    sb.append(String.format(",%.4f", c.getValue(p.positions)));
                sb.append("\n");
            }
            InputFileFormat.writeStringToDisk(sb.toString(), CSVfilename);
            System.out.printf("Wrote internal coordinate data to %s.\n", CSVfilename);
        }
    }
}
