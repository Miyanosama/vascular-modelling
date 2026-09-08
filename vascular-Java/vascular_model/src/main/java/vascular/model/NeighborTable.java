package vascular.model;

/**
 * Precomputed neighbor indices for diffusion on the voxel lattice.
 *
 * <p>For each voxel, stores the linear indices of its adjacent voxels.
 * In 3D: up to 6 neighbors (x±1, y±1, z±1).
 * In 2D: up to 4 neighbors (x±1, y±1).
 * Reflective (Neumann) boundary: voxels at domain edges simply have fewer neighbors.
 *
 * <p>Also classifies neighbors as longitudinal (z-direction) or lateral for
 * P-mediated transport calculations.
 */
public class NeighborTable {

    private final int totalVoxels;
    private final int maxNeighbors; // 6 for 3D, 4 for 2D

    // neighbors[v][0..neighborCount[v]-1] = linear indices of adjacent voxels
    private final int[][] neighbors;

    // Number of neighbors for each voxel
    private final int[] neighborCount;

    // Whether each neighbor is in the longitudinal direction (z-axis)
    private final boolean[][] isLongitudinal;

    /**
     * Builds the neighbor table for a given domain.
     */
    public NeighborTable(DomainMask domain) {
        this.totalVoxels = domain.totalVoxels();
        int[] dims = domain.dimensions();
        boolean is3D = dims[2] > 1;
        this.maxNeighbors = is3D ? 6 : 4;

        this.neighbors = new int[totalVoxels][maxNeighbors];
        this.neighborCount = new int[totalVoxels];
        this.isLongitudinal = new boolean[totalVoxels][maxNeighbors];

        // Direction vectors: right, left, down, up, forward, back
        int[][] dirs3D = {{1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1}};
        int[][] dirs2D = {{1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}};
        int[][] dirs = is3D ? dirs3D : dirs2D;

        for (int v = 0; v < totalVoxels; v++) {
            VoxelIndex vi = domain.fromLinear(v);
            int count = 0;

            for (int d = 0; d < dirs.length; d++) {
                int nx = vi.x() + dirs[d][0];
                int ny = vi.y() + dirs[d][1];
                int nz = vi.z() + dirs[d][2];

                if (domain.isInside(nx, ny, nz)) {
                    int nidx = domain.toLinear(nx, ny, nz);
                    if (nidx >= 0) {
                        neighbors[v][count] = nidx;
                        isLongitudinal[v][count] = is3D && (d >= 4); // z-directions are last 2 in dirs3D
                        count++;
                    }
                }
            }
            neighborCount[v] = count;
        }
    }

    /** Returns the neighbor indices for a voxel. */
    public int[] neighborsOf(int voxelIdx) {
        return neighbors[voxelIdx];
    }

    /** Number of neighbors for a voxel. */
    public int neighborCountOf(int voxelIdx) {
        return neighborCount[voxelIdx];
    }

    /** Whether the n-th neighbor of a voxel is in the longitudinal (z) direction. */
    public boolean isLongitudinal(int voxelIdx, int neighborSlot) {
        return isLongitudinal[voxelIdx][neighborSlot];
    }

    public int totalVoxels() { return totalVoxels; }
    public int maxNeighbors() { return maxNeighbors; }
}
