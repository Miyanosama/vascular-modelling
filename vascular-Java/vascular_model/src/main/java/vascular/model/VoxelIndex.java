package vascular.model;

/**
 * Immutable 3D index identifying a voxel in the lattice.
 * For 2D simulations, z is always 0.
 */
public record VoxelIndex(int x, int y, int z) {

    /**
     * Returns the linearized index for a lattice of given dimensions.
     */
    public int toLinear(int width, int height) {
        return (z * height + y) * width + x;
    }

    /**
     * Creates a VoxelIndex from a linear index.
     */
    public static VoxelIndex fromLinear(int linear, int width, int height) {
        int x = linear % width;
        int y = (linear / width) % height;
        int z = linear / (width * height);
        return new VoxelIndex(x, y, z);
    }

    @Override
    public String toString() {
        return "(" + x + ", " + y + ", " + z + ")";
    }
}
