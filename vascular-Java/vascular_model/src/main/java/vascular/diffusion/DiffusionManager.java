package vascular.diffusion;

import vascular.model.*;
import vascular.reaction.ReactionType;

/**
 * Manages diffusion propensities and executes diffusion events.
 *
 * <p>Each molecule of species S in a voxel diffuses with propensity D_S.
 * The destination is uniformly chosen among the voxel's neighbors.
 *
 * <p>For H with P-mediated transport (R8), the longitudinal movement
 * probability is boosted based on the local P concentration per Eq 3 of Hearn (2019):
 *
 *   P(longitudinal) = (2*DH + a*[P]) / (6*DH + a*[P])    (in 3D)
 *
 * In 2D, there are no longitudinal directions, so P has no effect on H diffusion.
 */
public class DiffusionManager {

    private final DomainMask domain;
    private final NeighborTable neighborTable;
    private final ModelParameters params;
    private final boolean is3D;

    public DiffusionManager(DomainMask domain, NeighborTable neighborTable, ModelParameters params) {
        this.domain = domain;
        this.neighborTable = neighborTable;
        this.params = params;
        this.is3D = domain.dimensions()[2] > 1;
    }

    /**
     * Computes the total diffusion propensity for all species at a voxel.
     */
    public double totalDiffusionPropensity(Lattice lattice, int voxelIdx) {
        double total = 0.0;

        // H diffusion (with P-mediated transport boost)
        int H = lattice.get(Species.H, voxelIdx);
        if (H > 0 && params.DH() > 0) {
            total += params.DH() * H;
        }

        // B diffusion
        int B = lattice.get(Species.B, voxelIdx);
        if (B > 0 && params.DB() > 0) {
            total += params.DB() * B;
        }

        // M diffusion
        int M = lattice.get(Species.M, voxelIdx);
        if (M > 0 && params.DM() > 0) {
            total += params.DM() * M;
        }

        // P does NOT diffuse (DP = 0 per paper)
        return total;
    }

    /**
     * Selects and executes a specific diffusion event at a voxel.
     *
     * <p>The event selection works as follows:
     * 1. Compute per-species total diffusion propensities
     * 2. Select a species proportionally
     * 3. Select a neighbor direction
     * 4. Move one molecule
     *
     * @param lattice   the current lattice state
     * @param voxelIdx  the source voxel
     * @param random    uniform random number in [0, totalPropensity)
     * @param random2   second uniform random number for neighbor selection
     * @return information about which species moved and where, or null if no diffusion possible
     */
    public DiffusionEvent executeDiffusion(
            Lattice lattice, int voxelIdx, double random, double random2) {

        int H = lattice.get(Species.H, voxelIdx);
        int B = lattice.get(Species.B, voxelIdx);
        int M = lattice.get(Species.M, voxelIdx);

        double hProp = (H > 0 && params.DH() > 0) ? params.DH() * H : 0;
        double bProp = (B > 0 && params.DB() > 0) ? params.DB() * B : 0;
        double mProp = (M > 0 && params.DM() > 0) ? params.DM() * M : 0;
        double total = hProp + bProp + mProp;

        if (total <= 0 || H + B + M == 0) return null;

        // Select which species diffuses
        double cumulative = 0;
        Species movingSpecies;
        if (random < hProp) {
            movingSpecies = Species.H;
        } else if (random < hProp + bProp) {
            movingSpecies = Species.B;
        } else {
            movingSpecies = Species.M;
        }

        // Select neighbor direction
        int numNeighbors = neighborTable.neighborCountOf(voxelIdx);
        if (numNeighbors == 0) return null;

        int targetVoxel;
        boolean isLongitudinal;

        if (movingSpecies == Species.H && is3D && params.a() > 0) {
            // P-mediated transport: select direction with biased probabilities
            int P = lattice.get(Species.P, voxelIdx);
            targetVoxel = selectNeighborWithPTransport(voxelIdx, P, numNeighbors, random2);
        } else {
            // Uniform among all neighbors
            int slot = (int) (random2 * numNeighbors);
            if (slot >= numNeighbors) slot = numNeighbors - 1;
            targetVoxel = neighborTable.neighborsOf(voxelIdx)[slot];
        }

        // Determine if this is a longitudinal movement
        // Find which slot corresponds to targetVoxel
        isLongitudinal = false;
        int[] neighbors = neighborTable.neighborsOf(voxelIdx);
        for (int s = 0; s < numNeighbors; s++) {
            if (neighbors[s] == targetVoxel) {
                isLongitudinal = neighborTable.isLongitudinal(voxelIdx, s);
                break;
            }
        }

        // Execute the move
        lattice.adjust(movingSpecies, voxelIdx, -1);
        lattice.adjust(movingSpecies, targetVoxel, 1);

        return new DiffusionEvent(movingSpecies, voxelIdx, targetVoxel, isLongitudinal);
    }

    /**
     * Selects a neighbor direction for H with P-mediated transport bias.
     *
     * <p>Per Eq 3 of the paper:
     *   P(longitudinal) = (2DH + a*[P]) / (6DH + a*[P])
     *   P(lateral, per direction) = DH / (6DH + a*[P])
     *
     * <p>So each of the 2 longitudinal directions gets weight (2DH + a*[P])/(6DH + a*[P]) / 2
     * and each of the 4 lateral directions gets weight DH/(6DH + a*[P]).
     *
     * <p>For simplicity, we assign relative weights:
     *   longitudinal neighbor: weight = (2DH + a*[P]) / 2
     *   lateral neighbor:     weight = DH
     * Then normalize.
     */
    private int selectNeighborWithPTransport(int voxelIdx, int P, int numNeighbors, double random) {
        int[] neighbors = neighborTable.neighborsOf(voxelIdx);
        double DH = params.DH();
        double a = params.a();

        double longWeight = (2.0 * DH + a * P) / 2.0;
        double latWeight = DH;

        // Compute total weight
        double totalWeight = 0;
        double[] weights = new double[numNeighbors];
        for (int s = 0; s < numNeighbors; s++) {
            weights[s] = neighborTable.isLongitudinal(voxelIdx, s) ? longWeight : latWeight;
            totalWeight += weights[s];
        }

        // Select based on cumulative weights
        double target = random * totalWeight;
        double cum = 0;
        for (int s = 0; s < numNeighbors; s++) {
            cum += weights[s];
            if (target < cum) {
                return neighbors[s];
            }
        }
        // Fallback: last neighbor
        return neighbors[numNeighbors - 1];
    }

    public NeighborTable neighborTable() { return neighborTable; }
}
