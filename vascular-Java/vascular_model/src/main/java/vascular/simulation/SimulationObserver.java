package vascular.simulation;

import vascular.model.Lattice;

/**
 * Callback interface for monitoring simulation progress.
 * Observers should complete quickly and not block the simulation thread.
 */
@FunctionalInterface
public interface SimulationObserver {
    /**
     * Called periodically during simulation.
     *
     * @param eventCount total events processed so far
     * @param simTime    current simulation time
     * @param lattice    current lattice state (read-only snapshot)
     */
    void onStep(long eventCount, double simTime, Lattice lattice);

    /** Called when simulation is paused (e.g., for checkpointing). */
    default void onCheckpoint(long eventCount, double simTime) {}

    /** Called when simulation completes successfully. */
    default void onComplete(long eventCount, double simTime) {}
}
