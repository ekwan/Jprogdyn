package edu.harvard.chemistry.ekwan.Jprogdyn;

import java.util.concurrent.*;

/**
 * This interface represents a piece of work that can be done in parallel.
 */
public interface WorkUnit extends Callable<Result>
{
    /** Perform the work in this unit. */
    public Result call();
}
