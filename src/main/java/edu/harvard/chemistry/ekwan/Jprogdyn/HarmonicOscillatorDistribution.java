package edu.harvard.chemistry.ekwan.Jprogdyn;

import java.util.*;
import com.google.common.collect.*;
import org.apache.commons.math3.analysis.polynomials.*;
import java.util.concurrent.*;
import java.math.BigDecimal;

/**
 * This class draws random numbers from the eigenfunctions of the harmonic oscillator.
 * This class has no state.
 */
public class HarmonicOscillatorDistribution {

    /** Energy levels higher than this will be evaluated classically.  Attempting to go too high here will blow out the numerical precision of double. */
    public static final int MAX_LEVEL = 20;

    /** Maximum number of times to try drawing a quantum displacement. */
    public static final int MAX_QUANTUM_ATTEMPTS = 10000;

    /** Prevents two Hermite polynomials from being generated at the same time. */
    private static final Object HERMITE_LOCK = new Object();

    /**
     * Draws a random displacement from the quantum distribution for the specified oscillator.  The answer
     * must be within the classically allowed region.  Rejection sampling is used.  This class is thread-safe.
     * @param n the vibrational energy level (0, 1, ...)
     * @param mass the reduced mass in amu
     * @param wavenumber the frequency of this mode in cm^-1
     * @param forceConstant the force constant in mDyne/A
     * @return the random displacement in angstroms
     */
    public static double drawRandomQuantumDisplacement(int n, double mass, double wavenumber, double forceConstant)
    {
        //System.out.printf("n=%d  mass=%.1f  freq=%.1f  force_constant=%.1f\n", n, mass, wavenumber, forceConstant);
        if ( wavenumber < 0.0 )
            throw new IllegalArgumentException("negative frequencies not allowed");

        // compute the classical turning point
        double ZPE = Initializer.getZeroPointEnergy(wavenumber); // kcal/mol
        double energy = ZPE * (2.0 * n + 1.0);                   // kcal/mol
        double turningPoint = getClassicalTurningPoint(energy, forceConstant);

        // if this is a high energy level, draw a classical displacement
        if ( n > MAX_LEVEL )
            return drawRandomClassicalDisplacement(energy, forceConstant); 

        // the are the minimum and maximum y-values of the distribution
        double minValue = 0.0;
        double maxValue = getQuantumDistributionMaxValue(n, mass, wavenumber, turningPoint);
        
        int attempts = 0;
        while (attempts < MAX_QUANTUM_ATTEMPTS )
            {
                // draw a random x-value between -turningPoint and +turningPoint
                // use MAX_A so we don't go right to the edge
                double x = getDouble(-turningPoint, turningPoint);

                // evaluate the function value
                double functionValue = getQuantumDistributionValue(n, mass, wavenumber, x);

                // draw a random number between minValue and maxValue
                double y = getDouble(minValue, maxValue);

                // accept or reject
                if ( y <= functionValue )
                    return x;

                // keep track of how many times we've tried
                attempts++;
            }
        throw new IllegalArgumentException("classical rejection sampling failed");
    }

    /**
     * Finds the maximum value of the specified quantum harmonic oscillator distribution.
     * @param n the vibrational energy level (0, 1, ...)
     * @param mass the reduced mass in amu
     * @param wavenumber the frequency of this mode in cm^-1
     * @param maxShift the maximum positive displacement in angstroms
     * @return the maximum y-value
     */
    public static double getQuantumDistributionMaxValue(int n, double mass, double wavenumber, double maxShift)
    {
        if ( wavenumber < 0.0 )
            throw new IllegalArgumentException("negative frequencies not allowed");
        if ( n == 0 )
            return getQuantumDistributionValue(0, mass, wavenumber, 0);
        double stepSize = maxShift / (MAX_LEVEL * 100);
        double maxValue = 0.0;
        for (double x = -1.0 * maxShift; x <= 0; x += stepSize)
            {
                double value = getQuantumDistributionValue(n, mass, wavenumber, x);
                if ( value > maxValue )
                    maxValue = value;
            }
        return maxValue;
    }

    /**
     * Computes the value of |psi_n|^2, where psi is an eigenstate of the quantum harmonic oscillator.
     * @param n the vibrational energy level (0, 1, ...)  Should not exceed {@link #MAX_LEVEL}.
     *        If it does, use {@link #drawRandomClassicalDisplacement(double,double)}
     * @param mass the reduced mass in amu
     * @param wavenumber the frequency of this mode in cm^-1
     * @param x the displacement in angstroms
     * @return the value of the probability distribution in angstroms^-1
     */
    public static double getQuantumDistributionValue(int n, double mass, double wavenumber, double x)
    {
        if ( wavenumber < 0.0 )
            throw new IllegalArgumentException("frequency cannot be negative");
        if ( n < 0 )
            throw new IllegalArgumentException("n cannot be negative");
        if ( n > MAX_LEVEL )
            throw new IllegalArgumentException("exceeded maximum excitation level");

        // calculate 1 / (2^n n!)
        double term1 = 1.0 / ( Math.pow(2.0,n) * factorial(n) );

        // calculate sqrt(m*omega/pi*hbar)
        // where omega = 2 * pi * c * wavenumber
        //
        //  m * 2 * pi * c * wavenumber     4 * pi * m * c * wavenumber
        // ----------------------------- = -----------------------------
        //       pi * h / 2 * pi                       h
        //
        // dimensional analysis (this is the only term that has units):
        // we need the sqrt to work out to angstroms^-2
        // 
        //        mass g       kg             mol           C cm     wavenumber          s             m^2
        // 4pi * -------- * ------- * ------------------ * ------ * ------------ * ------------- * -----------
        //         mol       1000 g    AVOGADROS_NUMBER      s          cm          H m^2 * kg      1E20 A^2
        //
        //      4 * pi * mass * C * wavenumber
        // = ------------------------------------
        //    1000 * AVOGADROS_NUMBER * H * 1E20 
        //
        // note that we also need this quantity in the term3 argument, just times pi
        double temp = 4.0 * Math.PI * mass * Units.C * wavenumber / (Units.AVOGADROS_NUMBER * Units.H * 1E23);
        double term2 = Math.sqrt(temp);

        // exp[(-m omega / 2hbar) * x^2)]
        double term3 = Math.exp(-1.0 * temp * Math.PI * x * x);

        // Hermite polynomial of degree n, where x = sqrt(m * omega / hbar) x
        // note that the Hermite generator is not thread safe
        PolynomialFunction function = null;
        synchronized (HERMITE_LOCK)
            {
                function = PolynomialsUtils.createHermitePolynomial(n);
            }
        double argument = Math.sqrt(temp * Math.PI) * x;
        double term4 = function.value(argument);
        term4 = term4 * term4;
        return term1 * term2 * term3 * term4;
    }

    /**
     * Helper method that returns factorials.
     * @param n the number to return the factorial of
     * @return n!
     */
    public static double factorial(int n)
    {
        if ( n < 0 )
            throw new IllegalArgumentException("n must not be negative");

        double factorial = 1.0;
        for (int i=1; i <= n; i++)
            factorial = factorial * i;
        return factorial;
    }

    /** Maximum number of times to try drawing a classical displacement. */
    public static final int MAX_CLASSICAL_ATTEMPTS = 10000;

    /**
     * Draws a random displacement from a uniform distribution between classical turning points.
     * It is recommended that you do not use the program if you cannot understand this function.
     * 
     * @param energy the total energy in the mode in kcal/mol
     * @param forceConstant the force constant of the mode in mDyne/A
     * @return a random x-value in angstroms, which should be uniformly distributed.
     */
    public static double drawRandomUniformDisplacement(double energy, double forceConstant)
    {
        double A = getClassicalTurningPoint(energy, forceConstant);
        return getDouble(-A, A);
    }

    /**
     * Draws a random displacement from the classical distribution for the specified oscillator.  The answer
     * must be within the classically allowed region.  Rejection sampling is used.  This class is thread-safe.
     *
     * @param energy the total energy in the mode in kcal/mol
     * @param forceConstant the force constant of the mode in mDyne/A
     * @return a random x-value in angstroms, which should be concentrated near the classical turning points
     */
    public static double drawRandomClassicalDisplacement(double energy, double forceConstant)
    {
        // compute the classical turning point
        double A = getClassicalTurningPoint(energy, forceConstant);
        double minA = -1.0 * MAX_A * A;
        double maxA = MAX_A * A;

        // compute the minimum value of the distribution, which occurs at x=0
        double minValue = getClassicalDistributionValue(A, 0.0);

        // compute the maximum value of the distribution, using A*MAX_A so we don't go right to the edge
        double maxValue = getClassicalDistributionValue(A, maxA);

        int attempts = 0;
        while (attempts < MAX_CLASSICAL_ATTEMPTS )
            {
                // draw a random x-value between -A and +A
                // use MAX_A so we don't go right to the edge
                double x = getDouble(minA, maxA);

                // evaluate the function value
                double functionValue = getClassicalDistributionValue(A, x);

                // draw a random number between minValue and maxValue
                double y = getDouble(minValue, maxValue);

                // accept or reject
                if ( y <= functionValue )
                    return x;

                // keep track of how many times we've tried
                attempts++;
            }
        throw new IllegalArgumentException("classical rejection sampling failed");
    }

    /** Prevents getClassicalDistributionValue(double,double,double) from blowing up if x is set to A. */
    public static final double MAX_A = 0.999;

    /**
     * Computes the classical probability distribution function for the quantum harmonic oscillator.
     * The classical distribution is (Robinett, R.W.  Am. J. Phys.  1995, 63(9), 823):
     *
     *  1             1
     * ---- * ----------------- , where A is the classical turning point.  Note this diverges as x approaches A.
     *  pi     sqrt(A^2 - x^2)
     *
     * @param A the classical turning point in angstroms
     * @param x the x-coordinate in angstroms to evaluate the value of the probability distribution
     * @return the value of the probability distribution
     */
    public static double getClassicalDistributionValue(double A, double x)
    {
        // throw an exception if outside classically allowed region
        double x1 = x;
        if ( x1 > A || x1 < -1.0 * A )
            throw new IllegalArgumentException("outside classically allowed region");

        // if we are close to the edge, return something close to the edge to prevent the function value from blowing up
        if ( x1 > MAX_A * A )
            x1 = MAX_A*A;
        else if ( x1 < -1.0 * MAX_A * A )
            x1 = -1.0 * MAX_A * A;

        return (1.0 / Math.PI) * ( 1.0 / Math.sqrt( A*A - x1*x1 ) );
    }

    /**
     * Computes the classical turning point in a harmonic oscillator.
     *
     * total energy = 1 / ( 2.0 * forceConstant * turningPoint^2)
     * turningPoint = sqrt(2E/k)
     *
     * dimensional analysis (1 dyne = 1E-5 N)
     *
     *  energy     kcal       A        mdyne     N s^2        m       KCAL_PER_MOL_TO_J J mol     kg m^2     1E20 A^2
     * -------- = ------ * ------- * -------- * ------- * -------- * ------------------------- * -------- * ----------
     *    k        mol      mDyne     1E-8 N     kg m      1E10 A              kcal               s^2 J        m^2
     *
     * = KCAL_PER_MOL_TO_J * 1E18
     *
     * @param energy the total energy in the mode in kcal/mol
     * @param forceConstant the force constant of the mode in mDyne/A
     * @return the classical turning point (in the +x direction, angstroms)
     */
    public static double getClassicalTurningPoint(double energy, double forceConstant)
    {
        return Math.sqrt( (2.0 * energy / forceConstant) * (Units.KCAL_PER_MOL_TO_J * 1E18) );
    }

    public static double getDouble(double min, double max)
    {
        if ( min > max )
            throw new IllegalArgumentException("min greater than max" + min + "," + max);
        return min + (max-min) * ThreadLocalRandom.current().nextDouble();
    }
}
