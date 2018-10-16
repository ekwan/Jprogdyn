package edu.harvard.chemistry.ekwan.Jprogdyn;

import java.util.*;
import java.io.*;
import com.google.common.collect.*;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

/**
 * This class collects together methods for analyzing NMR trajectories.
 */
public class NMRTrajectoryAnalyzer implements Immutable, Singleton {
    /** Not instantiable. */
    private NMRTrajectoryAnalyzer() {
        throw new IllegalArgumentException("not instantiable");
    }

    /** Represents the raw corrections obtained by averaging over a number of trajectories. */
    public static class NMRtrajectoryAnalysis implements Immutable {
        /** Parallel list containing the element of the nucleus. */
        public final List<Element> elements;
        
        /** Parallel list containing the atom numbers over which each correction has been averaged. 1, 2, ..., n. */
        public final List<List<Integer>> atomNumbers;
    
        /** Parallel list containing the mean raw corrections in ppm. */
        public final List<Double> means;

        /** Parallel list containing the standard deviations in ppm. */
        public final List<Double> standardDeviations;

        /** Parallel list containing the standard errors in ppm. */
        public final List<Double> standardErrors;

        /**
         * Make an analysis object.
         * @param elements the kinds of atoms being analyzed (parallel list)
         * @param atomNumbers list of 1-indexed atom numbers (parallel list)
         * @param means the mean raw corrections in ppm (parallel list)
         * @param standardDeviations the standard deviations in the mean raw corrections corrections (parallel list)
         * @param standardErrors the standard errors of the mean raw corrections (parallel list)
         */
        public NMRtrajectoryAnalysis(List<Element> elements, List<List<Integer>> atomNumbers,
                                     List<Double> means, List<Double> standardDeviations, List<Double> standardErrors) {
            if ( ImmutableSet.of(elements.size(), atomNumbers.size(), means.size(), standardDeviations.size()).size() != 1 )
                throw new IllegalArgumentException("size mismatch");
            if ( elements == null || atomNumbers == null || means == null || standardDeviations == null )
                throw new NullPointerException("no nulls allowed");
            for (List<Integer> list : atomNumbers)
                if ( list == null || list.size() == 0 )
                    throw new NullPointerException("no empties allowed");
            this.elements = ImmutableList.copyOf(elements);
            List<List<Integer>> tempList = new ArrayList<>(atomNumbers.size());
            for (List<Integer> list : atomNumbers)
                {
                    List<Integer> sorted = new ArrayList<Integer>(list);
                    Collections.sort(sorted);
                    tempList.add(ImmutableList.copyOf(sorted));
                }
            this.atomNumbers = ImmutableList.copyOf(tempList);
            this.means = ImmutableList.copyOf(means);
            this.standardDeviations = ImmutableList.copyOf(standardDeviations);
            this.standardErrors = ImmutableList.copyOf(standardErrors);
        }

        @Override
        public String toString() {
            String returnString = "Element:  Atom Numbers:                   Raw Correction:     Std. Err.:      Std. Dev.:\n";
            List<String> strings = new ArrayList<>(elements.size());
            for (int i=0; i < elements.size(); i++)
                {
                    Element element = elements.get(i);
                    List<Integer> thisAtomNumbers = atomNumbers.get(i);
                    Double mean = means.get(i);
                    Double standardDeviation = standardDeviations.get(i);
                    Double standardError = standardErrors.get(i);

                    String atomString = "";
                    for (Integer j : thisAtomNumbers)
                        atomString += String.format("%3d ", j);
                    returnString += String.format("  %2s      %-30s     %7.2f            %5.3f        %7.2f\n",
                                                  element.toString(), atomString, mean, standardError, standardDeviation);
                }
            Collections.sort(strings);
            for (String s : strings)
                returnString += s;
            return returnString;
        }
    }
    
    /**
     * Analyzes the given trajectories.
     * @param trajectories the trajectories to analyze
     * @param ignoreIncomplete whether to include incomplete trajectories in the analysis
     * @param symmetryList a list of lists containing the symmetry-equivalent atom _numbers_ (e.g. [ [1,2,3], [4,5] would mean average atoms 1,2,3 and 4,5)
     * @param molecule the molecule that contains the NMR shifts calculated at the stationary point geometry (needed to calculate raw corrections)
     * @return the analysis
     */
    public static NMRtrajectoryAnalysis analyze(List<Trajectory> trajectories, boolean ignoreIncomplete, List<Set<Integer>> symmetryList, Molecule molecule)
    {
        // initialize lists
        List<Element> elements = new ArrayList<>();
        int size = molecule.contents.size();
        List<List<Integer>> atomNumbers = new ArrayList<>();
        List<Integer> included = new ArrayList<>();
        for (Set<Integer> set : symmetryList)
            {
                if ( set.size() == 0 )
                    throw new IllegalArgumentException("unexpected empty");
                
                // put down the atom numbers in each list into a master list so we can make singly-occupied lists for all the other atom numbers
                List<Integer> list = new ArrayList<>(set);
                atomNumbers.add(list);
                included.addAll(set);
    
                // check we are averaging nuclei of the same element
                Set<Element> theseElements = new HashSet<>();
                for (Integer i : list)
                    {
                        Atom a = molecule.contents.get(i-1);
                        theseElements.add(Element.getElement(a.symbol));
                    }
                if ( theseElements.size() != 1 )
                    throw new IllegalArgumentException("check symmetry list: " + set.toString());

                // ok, so make a note of the element
                int atomIndex = list.get(0)-1;
                Atom a = molecule.contents.get(atomIndex);
                elements.add(Element.getElement(a.symbol));
            }
        for (int i=1; i <= size; i++)
            {
                if ( included.contains(i) )
                    continue;
                atomNumbers.add(ImmutableList.of(i));
                Atom a = molecule.contents.get(i-1);
                elements.add(Element.getElement(a.symbol));
            }
        List<List<Double>> meansByTrajectory = new ArrayList<>();       // means for each trajectory, index parallel to atomNumbers
        for (int i=0; i < atomNumbers.size(); i++)
            meansByTrajectory.add(new ArrayList<Double>());
        
        // collect data
        int count = 0;
        System.out.println(trajectories.size() + " total trajs");
        for (Trajectory trajectory : trajectories)
            {
                if ( ignoreIncomplete && !trajectory.isDone() ) 
                    {
                        System.out.printf("Skipping incomplete trajectory %s.\n", trajectory.checkpointFilename);
                        continue;
                    }
                
                // collect all the shieldings in this trajectory together
                List<TrajectoryPoint> allPoints = new ArrayList<>();
                allPoints.add(trajectory.initialPoint);
                allPoints.addAll(trajectory.forwardPoints);
                allPoints.addAll(trajectory.backwardPoints);
                List<List<Double>> trajectoryShieldings = new ArrayList<>();  // all shieldings in trajectory, outer index is trajectory point, inner index is atom
                for (TrajectoryPoint p : allPoints)
                    {
                        if (p.shieldings != null && p.shieldings.size() > 0)
                            trajectoryShieldings.add(p.shieldings);
                    }
                //System.out.println(trajectoryShieldings.size() + " shieldings in this traj");
                if (trajectoryShieldings.size() <= 20)  // if there aren't a minimum number of shieldings, keep going
                    continue;
                count++;

                // calculate the mean shielding accounting for symmetry in this trajectory

                for (int i=0; i < atomNumbers.size(); i++)
                    {
                        // collect all the shieldings for these atom numbers
                        List<Integer> atomNumbersToAverage = atomNumbers.get(i);
                        List<Double> shieldingsToAverage = new ArrayList<>();
                        for (Integer atomNumber : atomNumbersToAverage)
                            {
                                for (List<Double> pointShieldings : trajectoryShieldings)
                                    {
                                        Double thisShielding = pointShieldings.get(atomNumber-1);
                                        shieldingsToAverage.add(thisShielding);
                                    }
                            }

                        // add the mean
                        double thisMean = average(shieldingsToAverage);
                        meansByTrajectory.get(i).add(thisMean);
                    }
            }
        if ( count <= 1 )
            throw new IllegalArgumentException("not enough trajectories to proceed");

        // get the stationary point shieldings
        List<Double> stationaryPointShieldings = new ArrayList<>();
        for (List<Integer> atomNumbersToAverage : atomNumbers)
            {
                List<Double> shieldingsToAverage = new ArrayList<>();
                for (Integer atomNumber : atomNumbersToAverage)
                    shieldingsToAverage.add( molecule.shieldings.get(atomNumber-1) );
                double averageShielding = average(shieldingsToAverage);
                stationaryPointShieldings.add(averageShielding);
            }

        // do statistics
        List<Double> overallMeans = new ArrayList<>();
        List<Double> standardDeviations = new ArrayList<>();
        List<Double> standardErrors = new ArrayList<>();
        for (int i=0; i < meansByTrajectory.size(); i++)
            {
                DescriptiveStatistics stats = new DescriptiveStatistics();
                List<Double> theseMeans = meansByTrajectory.get(i); // shieldings for this set of atom numbers, indexed by trajectory
                double stationaryPointShielding = stationaryPointShieldings.get(i);

                // compute the raw correction
                for (Double d : theseMeans)
                    stats.addValue(d-stationaryPointShielding);
                overallMeans.add(stats.getMean());
                standardDeviations.add(stats.getStandardDeviation());
                double standardError = stats.getStandardDeviation() / Math.sqrt(theseMeans.size());
                standardErrors.add(standardError);
            }

        // return result
        return new NMRtrajectoryAnalysis(elements, atomNumbers, overallMeans, standardDeviations, standardErrors);
    }

    /**
     * Averages the numbers in a list.
     * @param list the numbers to average
     * @return the average
     */
    public static Double average(List<Double> list)
    {
        if ( list == null || list.size() == 0 )
            throw new NullPointerException("empty list for averaging");
        double sum = 0.0;
        for (Double d : list)
            sum += d;
        return sum / list.size();
    }
}
