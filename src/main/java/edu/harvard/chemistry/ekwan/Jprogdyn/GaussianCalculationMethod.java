package edu.harvard.chemistry.ekwan.Jprogdyn;

import java.util.*;

/**
 * This class represents a Gaussian calculation.
 * Note that NMR calculations often require forces.
 */
public class GaussianCalculationMethod extends CalculationMethod {
    
    /** For serialization. */
    public static final long serialVersionUID = 1L;

    /** The route card.  Can contain multiple lines.  A newline will automatically be appended. */
    public final String routeCard;

    /** Will be appended at the end of the input file.  A newline will automatically be appended. */
    public final String footer;

    /**
     * Constructor.
     * @param calculationType whether this is an NMR or force calculation
     * @param memory memory to use in GB
     * @param processors number of processors to use
     * @param routeCard the route card
     * @param footer any text to append after the geometry
     */
    public GaussianCalculationMethod(CalculationMethod.CalculationType calculationType, int memory,
                                     int processors, String routeCard, String footer) {
        super(calculationType, CalculationMethod.Program.GAUSSIAN, memory, processors);
        if ( routeCard == null || routeCard.trim().length() == 0 )
            throw new IllegalArgumentException("blank route card");
        this.routeCard = routeCard;
        this.footer = footer;
    }

    @Override
    public String toString() {
        String returnString = String.format("G09 calculation (%s, %d GB, %d processors)\n", calculationType.toString(), memory, processors);
        returnString += String.format("Route card: %s\nFooter:%s\n", routeCard, footer);
        return returnString;
    }

    @Override
    public int hashCode() {
        return Objects.hash(calculationType, program, memory, processors, routeCard, footer);
    }

    @Override
    public boolean equals(Object obj) {
        if ( obj == null )
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof GaussianCalculationMethod) )
            return false;

        GaussianCalculationMethod m = (GaussianCalculationMethod)obj;
        if ( m.calculationType == calculationType &&
             m.program == program &&
             m.memory == memory &&
             m.processors == processors &&
             Objects.equals(m.routeCard, routeCard) &&
             Objects.equals(m.footer, footer) )
            return true;
        return false;
    }
}
