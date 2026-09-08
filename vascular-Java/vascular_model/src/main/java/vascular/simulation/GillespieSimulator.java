package vascular.simulation;

import vascular.diffusion.DiffusionEvent;
import vascular.diffusion.DiffusionManager;
import vascular.model.*;
import vascular.reaction.ReactionPropensity;
import vascular.reaction.ReactionType;

import java.util.random.RandomGenerator;

/**
 * Gillespie Stochastic Simulation Algorithm for the spatial reaction-diffusion model.
 *
 * <p>Core algorithm:
 * 1. Compute total propensity a0 across all voxels (reactions + diffusion)
 * 2. Sample time to next event: tau = -ln(U1) / a0
 * 3. Select voxel by rejection sampling, then select event within voxel
 * 4. Execute event (reaction or diffusion)
 * 5. Update propensities only in affected voxels (O(1) update)
 * 6. Notify observers periodically
 *
 * <p>This implements the full spatial SSA with constant-time event selection.
 */
public class GillespieSimulator {

    private final DomainMask domain;
    private final NeighborTable neighborTable;
    private final AnnularRegionMask annularMask;
    private final DiffusionManager diffusionManager;
    private final ModelParameters params;
    private final Lattice lattice;

    private final RandomGenerator rng;

    private double totalPropensity; // a0: sum of all reaction + diffusion propensities
    private double simTime;
    private long eventCount;

    private boolean running;
    private volatile boolean paused;
    private volatile boolean stopped;

    // Observation settings
    private SimulationObserver observer;
    private long observerInterval = 1_000_000; // notify every N events
    private long nextObserverEvent;

    // Maximums
    private double maxTime = Double.POSITIVE_INFINITY;
    private long maxEvents = Long.MAX_VALUE;

    // For forward-lookup: map voxel -> region bit (outer shell or not)
    private final boolean[] isOuterShell;

    /**
     * Creates a new Gillespie simulator.
     *
     * @param domain     spatial domain mask
     * @param params     model parameters (rate constants, diffusion coefficients)
     * @param seed       random seed for reproducibility
     */
    public GillespieSimulator(DomainMask domain, ModelParameters params, long seed) {
        this.domain = domain;
        this.params = params;
        this.lattice = new Lattice(domain);
        this.neighborTable = new NeighborTable(domain);
        this.annularMask = new AnnularRegionMask(domain, 0.96); // outer 4%
        this.diffusionManager = new DiffusionManager(domain, neighborTable, params);
        this.rng = new java.util.Random(seed);

        int totalVoxels = domain.totalVoxels();
        this.isOuterShell = new boolean[totalVoxels];
        for (int v = 0; v < totalVoxels; v++) {
            VoxelIndex vi = domain.fromLinear(v);
            this.isOuterShell[v] = annularMask.isInOuterShell(vi.x(), vi.y(), vi.z());
        }

        this.totalPropensity = 0.0;
        this.simTime = 0.0;
        this.eventCount = 0;
        this.running = false;
        this.paused = false;
        this.stopped = false;
    }

    // ---------- Configuration ----------

    public void setObserver(SimulationObserver observer) { this.observer = observer; }
    public void setObserverInterval(long interval) { this.observerInterval = interval; }
    public void setMaxTime(double maxTime) { this.maxTime = maxTime; }
    public void setMaxEvents(long maxEvents) { this.maxEvents = maxEvents; }

    public Lattice lattice() { return lattice; }
    public double simTime() { return simTime; }
    public long eventCount() { return eventCount; }
    public double totalPropensity() { return totalPropensity; }

    // ---------- Control ----------

    /** Starts or resumes the simulation. */
    public void start() {
        if (!running) {
            running = true;
            stopped = false;
            paused = false;
            if (eventCount == 0) {
                initializePropensities();
            }
            runLoop();
        } else if (paused) {
            paused = false;
        }
    }

    /** Pauses the simulation (may continue later). */
    public void pause() {
        paused = true;
    }

    /** Stops the simulation permanently. */
    public void stop() {
        stopped = true;
        paused = false;
    }

    public boolean isRunning() { return running && !paused; }
    public boolean isStopped() { return stopped; }

    // ---------- Initialization ----------

    private void initializePropensities() {
        int n = domain.totalVoxels();
        totalPropensity = 0.0;
        for (int v = 0; v < n; v++) {
            double rp = ReactionPropensity.totalReaction(lattice, v, params, isOuterShell[v]);
            double dp = diffusionManager.totalDiffusionPropensity(lattice, v);
            lattice.setReactionPropensity(v, rp);
            lattice.setDiffusionPropensity(v, dp);
            double vp = rp + dp;
            lattice.setVoxelPropensity(v, vp);
            totalPropensity += vp;
        }
        nextObserverEvent = observerInterval;
    }

    // ---------- Main Simulation Loop ----------

    private void runLoop() {
        int n = domain.totalVoxels();

        while (!stopped) {
            // Check pause
            while (paused && !stopped) {
                try { Thread.sleep(50); } catch (InterruptedException e) { Thread.currentThread().interrupt(); return; }
            }
            if (stopped) break;

            // Check termination conditions
            if ((maxTime > 0 && simTime >= maxTime) || eventCount >= maxEvents) {
                if (observer != null) observer.onComplete(eventCount, simTime);
                running = false;
                break;
            }

            // Step 1: Compute time to next event
            if (totalPropensity <= 0) {
                // No possible events — wait for spontaneous production
                // In practice, k2 (B production) and k4/k5 (H production/degradation)
                // ensure totalPropensity > 0 once molecules exist
                totalPropensity = 0.0;
                for (int v = 0; v < n; v++) {
                    double rp = ReactionPropensity.totalReaction(lattice, v, params, isOuterShell[v]);
                    double dp = diffusionManager.totalDiffusionPropensity(lattice, v);
                    double vp = rp + dp;
                    lattice.setReactionPropensity(v, rp);
                    lattice.setDiffusionPropensity(v, dp);
                    lattice.setVoxelPropensity(v, vp);
                    totalPropensity += vp;
                }
                if (totalPropensity <= 0) break;
            }

            double tau = -Math.log(rng.nextDouble()) / totalPropensity;
            simTime += tau;

            // Step 2: Select voxel (rejection sampling)
            int selectedVoxel = selectVoxel(n);
            if (selectedVoxel < 0) continue;

            double voxelProp = lattice.getVoxelPropensity(selectedVoxel);
            if (voxelProp <= 0) continue;

            // Step 3: Select event within voxel
            double eventRand = rng.nextDouble(voxelProp);
            double reactionPart = lattice.getReactionPropensity(selectedVoxel);
            double diffusionPart = lattice.getDiffusionPropensity(selectedVoxel);

            if (eventRand < reactionPart) {
                // Reaction event
                ReactionType reaction = selectReaction(selectedVoxel, eventRand);
                if (reaction != null) {
                    ReactionPropensity.execute(reaction, lattice, selectedVoxel);
                }
            } else {
                // Diffusion event
                double diffRand = eventRand - reactionPart;
                diffusionManager.executeDiffusion(lattice, selectedVoxel, diffRand, rng.nextDouble());
            }

            eventCount++;

            // Step 4: Update propensities in affected voxels
            updateAffectedPropensities(selectedVoxel, n);

            // Step 5: Notify observer
            if (observer != null && eventCount >= nextObserverEvent) {
                observer.onStep(eventCount, simTime, lattice);
                nextObserverEvent = eventCount + observerInterval;
            }
        }
        running = false;
        if (observer != null && !stopped) {
            observer.onComplete(eventCount, simTime);
        }
    }

    /**
     * Selects a voxel using rejection sampling (constant-time expected).
     * Uses the maximum per-voxel propensity as the bounding envelope.
     */
    private int selectVoxel(int numVoxels) {
        // Find max per-voxel propensity
        double maxProp = 0;
        for (int v = 0; v < numVoxels; v++) {
            double p = lattice.getVoxelPropensity(v);
            if (p > maxProp) maxProp = p;
        }
        if (maxProp <= 0) return -1;

        // Rejection sampling
        int maxAttempts = 1000;
        for (int attempt = 0; attempt < maxAttempts; attempt++) {
            int v = rng.nextInt(numVoxels);
            double p = lattice.getVoxelPropensity(v);
            if (rng.nextDouble() * maxProp < p) {
                return v;
            }
        }

        // Fallback: linear scan
        double target = rng.nextDouble(totalPropensity);
        double cum = 0;
        for (int v = 0; v < numVoxels; v++) {
            cum += lattice.getVoxelPropensity(v);
            if (target < cum) return v;
        }
        return numVoxels - 1;
    }

    /**
     * Selects a specific reaction type from the per-voxel reaction propensity breakdown.
     */
    private ReactionType selectReaction(int voxelIdx, double rand) {
        double cum = 0;
        for (ReactionType type : ReactionPropensity.REACTION_EVENTS) {
            double p = ReactionPropensity.compute(type, lattice, voxelIdx, params, isOuterShell[voxelIdx]);
            cum += p;
            if (rand < cum) return type;
        }
        // Fallback: return last reaction with non-zero propensity
        for (int i = ReactionPropensity.REACTION_EVENTS.length - 1; i >= 0; i--) {
            ReactionType t = ReactionPropensity.REACTION_EVENTS[i];
            if (ReactionPropensity.compute(t, lattice, voxelIdx, params, isOuterShell[voxelIdx]) > 0) {
                return t;
            }
        }
        return null;
    }

    /**
     * Updates propensities for the event voxel and all its neighbors.
     * This is O(1) since each voxel has at most 6 neighbors.
     */
    private void updateAffectedPropensities(int centerVoxel, int numVoxels) {
        // Update center voxel
        updateVoxelPropensity(centerVoxel);

        // Update all neighbors
        int numNeighbors = neighborTable.neighborCountOf(centerVoxel);
        int[] neighbors = neighborTable.neighborsOf(centerVoxel);
        for (int s = 0; s < numNeighbors; s++) {
            updateVoxelPropensity(neighbors[s]);
        }
    }

    private void updateVoxelPropensity(int v) {
        double rp = ReactionPropensity.totalReaction(lattice, v, params, isOuterShell[v]);
        double dp = diffusionManager.totalDiffusionPropensity(lattice, v);
        double oldVp = lattice.getVoxelPropensity(v);
        double newVp = rp + dp;
        lattice.setReactionPropensity(v, rp);
        lattice.setDiffusionPropensity(v, dp);
        lattice.setVoxelPropensity(v, newVp);
        totalPropensity += (newVp - oldVp);
    }

    /** Runs the simulation synchronously on the current thread (for headless/batch mode). */
    public void runSync() {
        if (eventCount == 0) initializePropensities();
        running = true;
        stopped = false;
        paused = false;
        runLoop();
    }
}
