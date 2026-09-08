package vascular.model;

/**
 * Defines which voxels are inside the simulation domain.
 * Supports 2D (circle cross-section) and 3D (right circular cylinder) geometries.
 *
 * <p>Z dimension is 1 for 2D simulations (single slice).
 */
public sealed interface DomainMask permits CircleMask, CylinderMask {

    /** Returns true if the voxel at (x, y, z) is inside the domain. */
    boolean isInside(int x, int y, int z);

    /** Total number of voxels inside the domain. */
    int totalVoxels();

    /** Grid dimensions: [width, height, depth]. */
    int[] dimensions();

    /** Domain radius in voxel units. */
    double radius();

    /** Center coordinates: [cx, cy]. */
    double[] center();

    /**
     * Maps a voxel index to its linear array index (0-based).
     * Returns -1 if the voxel is outside the domain.
     */
    int toLinear(int x, int y, int z);

    /**
     * Creates a VoxelIndex from a linear index.
     */
    VoxelIndex fromLinear(int linear);
}
