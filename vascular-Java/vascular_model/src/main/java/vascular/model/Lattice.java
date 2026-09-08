package vascular.model;

import java.util.Arrays;

/**
 * Holds the molecule counts for all species across all voxels in the domain.
 *
 * <p>Uses flat primitive arrays for cache-friendly memory layout.
 * Each voxel stores one int per species, plus per-voxel diffusion/reaction propensities.
 * Designed for high-throughput stochastic simulation with billions of events.
 */
public class Lattice {

    private final DomainMask domain;
    private final int totalVoxels;

    // Molecule counts: counts[species.index()][voxelLinearIndex]
    private final int[][] counts;

    // Per-voxel total propensity (sum of all reaction + diffusion propensities)
    private final double[] voxelPropensity;

    // Per-voxel propensity breakdown by event category
    private final double[] reactionPropensity;  // per voxel
    private final double[] diffusionPropensity; // per voxel

    /** Creates a lattice initialized to zero molecules. */
    public Lattice(DomainMask domain) {
        this.domain = domain;
        this.totalVoxels = domain.totalVoxels();
        this.counts = new int[Species.COUNT][totalVoxels];
        this.voxelPropensity = new double[totalVoxels];
        this.reactionPropensity = new double[totalVoxels];
        this.diffusionPropensity = new double[totalVoxels];
    }

    // ---------- Molecule count access ----------

    /** Returns the count of a species at a specific voxel. */
    public int get(Species species, int voxelIdx) {
        return counts[species.index()][voxelIdx];
    }

    /** Returns the count of a species at a grid position (via domain mapping). */
    public int get(Species species, int x, int y, int z) {
        int idx = domain.toLinear(x, y, z);
        return idx >= 0 ? counts[species.index()][idx] : 0;
    }

    /** Returns a read-only view of counts for a species across all voxels. */
    public int[] getCounts(Species species) {
        return counts[species.index()];
    }

    /** Sets the count of a species at a specific voxel. */
    public void set(Species species, int voxelIdx, int count) {
        if (count < 0) count = 0;
        counts[species.index()][voxelIdx] = count;
    }

    /** Increments a molecule count (delta can be negative). Clamps to 0. */
    public void adjust(Species species, int voxelIdx, int delta) {
        int idx = species.index();
        counts[idx][voxelIdx] = Math.max(0, counts[idx][voxelIdx] + delta);
    }

    /** Total number of a species across all voxels. */
    public int totalCount(Species species) {
        int sum = 0;
        for (int c : counts[species.index()]) sum += c;
        return sum;
    }

    // ---------- Propensity access ----------

    public double getVoxelPropensity(int voxelIdx) {
        return voxelPropensity[voxelIdx];
    }

    public void setVoxelPropensity(int voxelIdx, double p) {
        voxelPropensity[voxelIdx] = p;
    }

    public double getReactionPropensity(int voxelIdx) {
        return reactionPropensity[voxelIdx];
    }

    public void setReactionPropensity(int voxelIdx, double p) {
        reactionPropensity[voxelIdx] = p;
    }

    public double getDiffusionPropensity(int voxelIdx) {
        return diffusionPropensity[voxelIdx];
    }

    public void setDiffusionPropensity(int voxelIdx, double p) {
        diffusionPropensity[voxelIdx] = p;
    }

    public double[] getVoxelPropensities() { return voxelPropensity; }
    public double[] getReactionPropensities() { return reactionPropensity; }
    public double[] getDiffusionPropensities() { return diffusionPropensity; }

    // ---------- Domain access ----------

    public DomainMask domain() { return domain; }
    public int totalVoxels() { return totalVoxels; }

    /** Maximum H concentration across all voxels (for color scaling). */
    public int maxCount(Species species) {
        int max = 0;
        for (int c : counts[species.index()]) {
            if (c > max) max = c;
        }
        return max;
    }

    /** Creates a deep copy. */
    public Lattice copy() {
        Lattice other = new Lattice(domain);
        for (int s = 0; s < Species.COUNT; s++) {
            System.arraycopy(counts[s], 0, other.counts[s], 0, totalVoxels);
        }
        System.arraycopy(voxelPropensity, 0, other.voxelPropensity, 0, totalVoxels);
        System.arraycopy(reactionPropensity, 0, other.reactionPropensity, 0, totalVoxels);
        System.arraycopy(diffusionPropensity, 0, other.diffusionPropensity, 0, totalVoxels);
        return other;
    }

    /** Resets all molecule counts to zero. */
    public void reset() {
        for (int s = 0; s < Species.COUNT; s++) {
            Arrays.fill(counts[s], 0);
        }
        Arrays.fill(voxelPropensity, 0);
        Arrays.fill(reactionPropensity, 0);
        Arrays.fill(diffusionPropensity, 0);
    }
}
