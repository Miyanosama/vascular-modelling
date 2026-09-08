package vascular.simulation;

import vascular.model.Lattice;
import vascular.model.ModelParameters;

/**
 * Mutable snapshot of the current simulation state.
 */
public class SimulationState {

    private final Lattice lattice;
    private final ModelParameters params;
    private double simTime;
    private long eventCount;

    public SimulationState(Lattice lattice, ModelParameters params) {
        this.lattice = lattice;
        this.params = params;
        this.simTime = 0.0;
        this.eventCount = 0;
    }

    public Lattice lattice() { return lattice; }
    public ModelParameters params() { return params; }
    public double simTime() { return simTime; }
    public long eventCount() { return eventCount; }

    public void advanceTime(double tau) { this.simTime += tau; }
    public void incrementEvents() { this.eventCount++; }
}
