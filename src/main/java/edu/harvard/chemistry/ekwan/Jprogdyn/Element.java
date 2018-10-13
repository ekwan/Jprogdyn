package edu.harvard.chemistry.ekwan.Jprogdyn;

import java.util.*;
import com.google.common.collect.*;

/**
 * This enum represents an element.
 */
public enum Element
{
    /** Hydrogen. */
    HYDROGEN("H", 1),

    /** Carbon. */
    CARBON  ("C", 6),

    /** Nitrogen. */
    NITROGEN("N", 7),

    /** Oxygen. */
    OXYGEN  ("O", 8),

    /** Sulfur. */
    SULFUR  ("S", 16);

    /** The atomic symbol. */
    public final String symbol;

    /** The atomic number. */
    public final int atomicNumber;

    /**
     * Constructs an element.
     * @param symbol the atomic symbol
     * @param atomicNumber the atomic number
     */
    Element(String symbol, int atomicNumber)
    {
        this.symbol = symbol;
        this.atomicNumber = atomicNumber;
    }

    /**
     * Identifies the element corresponding to a string.
     * @param symbol the symbol for the requested element (case-sensitive)
     * @return the corresponding Element enum element
     */
    public static Element getElement(String symbol)
    {
        for (Element e : Element.values()) {
            if ( e.symbol.equals(symbol) || Integer.toString(e.atomicNumber).equals(symbol) )
                return e;
        }
        throw new IllegalArgumentException(String.format("element not found for symbol %s", symbol));
    }

    /**
     * Returns a string representation of this Element.
     * @return the String
     */
    public String toString()
    {
        return symbol;
    }
}
