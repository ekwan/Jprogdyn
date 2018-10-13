package edu.harvard.chemistry.ekwan.Jprogdyn;

import java.util.*;

/**
 * This class lets takes a value and a range and draws a bar like this:
 *
 * [------*-----]
 *
 */
public class AsciiBar implements Immutable, Singleton {
    
    /** Not instantiable. */
    private AsciiBar() {
        throw new IllegalArgumentException("not instantiable");
    }

    /**
     * Prints a bar with a marker indicating the position of the value.
     * The returned string will have length+4 characters to accomodate the edges of the bar. 
     * @param value where the marker is
     * @param min the value of the left end of the bar
     * @param max the value of the right end of the bar
     * @param marker what the marker should be (one character)
     * @param notMarker what to put on the bar in places other than the marker
     * @param length the length of the bar not including the two characters that indicate the edges of the bar
     * @return the bar as a string
     */
    public static String make(double value, double min, double max, String marker, String notMarker, int length) {
        // check invariants
        if ( max <= min )
            throw new IllegalArgumentException("check bounds, max is less than min");
        if ( marker == null || marker.length() != 1 || notMarker == null || notMarker.length() != 1 )
            throw new NullPointerException("one character markers required");
        if ( length < 5 )
            throw new IllegalArgumentException("bar too short");

        // return special bars if the bounds are exceeded
        if ( value < min )
            {
                StringBuilder s = new StringBuilder(marker + "[");
                for (int i=0; i < length; i++)
                    s.append(notMarker);
                s.append("] ");
                return s.toString();
            }
        else if ( value > max )
            {
                StringBuilder s = new StringBuilder(" [");
                for (int i=0; i < length; i++)
                    s.append(notMarker);
                s.append("]" + marker);
                return s.toString();
            }

        // figure out where to put the marker
        double temp = (value-min)*(length-1) / (max-min);
        int markerPosition = (int)temp;

        // build the string
        StringBuilder s = new StringBuilder(" [");
        for (int i=0; i < markerPosition; i++)
            s.append(notMarker);
        s.append(marker);
        for (int i=markerPosition+1; i < length; i++)
            s.append(notMarker);
        s.append("] ");
        return s.toString();
    }

    /** Prints a bar with standard parameters. */
    public static String make(double value, double min, double max) {
        return make(value, min, max, "*", " ", 20);
    }
}
