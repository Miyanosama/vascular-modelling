package vascular.model;

/**
 * 3D right circular cylinder domain representing a plant stem segment.
 * The cylinder axis is along the z-direction.
 */
public final class CylinderMask implements DomainMask {

    private final int width, height, depth;
    private final double cx, cy, radius;
    private final int totalVoxels;
    private final int[] linearMap;   // grid index -> linear index, or -1 if outside
    private final int[] reverseMap;  // linear index -> grid index, O(1) lookup

    /**
     * @param diameter Cylinder diameter in voxels.
     * @param depth    Cylinder height (z-dimension) in voxels.
     */
    public CylinderMask(int diameter, int depth) {
        this.width = diameter;
        this.height = diameter;
        this.depth = depth;
        this.cx = (diameter - 1) / 2.0;
        this.cy = (diameter - 1) / 2.0;
        this.radius = diameter / 2.0;
        int totalGrid = width * height * depth;
        this.linearMap = new int[totalGrid];

        int count = 0;
        for (int z = 0; z < depth; z++) {
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    double dx = x - cx;
                    double dy = y - cy;
                    int idx = (z * height + y) * width + x;
                    if (dx * dx + dy * dy <= radius * radius) {
                        linearMap[idx] = count++;
                    } else {
                        linearMap[idx] = -1;
                    }
                }
            }
        }
        this.totalVoxels = count;

        // Build reverse lookup
        this.reverseMap = new int[totalVoxels];
        for (int z = 0; z < depth; z++) {
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    int idx = (z * height + y) * width + x;
                    int lin = linearMap[idx];
                    if (lin >= 0) {
                        reverseMap[lin] = idx;
                    }
                }
            }
        }
    }

    @Override
    public boolean isInside(int x, int y, int z) {
        if (x < 0 || x >= width || y < 0 || y >= height || z < 0 || z >= depth) return false;
        double dx = x - cx;
        double dy = y - cy;
        return dx * dx + dy * dy <= radius * radius;
    }

    @Override
    public int totalVoxels() { return totalVoxels; }

    @Override
    public int[] dimensions() { return new int[]{width, height, depth}; }

    @Override
    public double radius() { return radius; }

    @Override
    public double[] center() { return new double[]{cx, cy}; }

    @Override
    public int toLinear(int x, int y, int z) {
        if (x < 0 || x >= width || y < 0 || y >= height || z < 0 || z >= depth) return -1;
        return linearMap[(z * height + y) * width + x];
    }

    @Override
    public VoxelIndex fromLinear(int linear) {
        if (linear < 0 || linear >= totalVoxels) {
            throw new IndexOutOfBoundsException("Linear index " + linear + " out of [0, " + totalVoxels + ")");
        }
        int gridIdx = reverseMap[linear];
        int x = gridIdx % width;
        int y = (gridIdx / width) % height;
        int z = gridIdx / (width * height);
        return new VoxelIndex(x, y, z);
    }

    /** Grid width (x-dimension). */
    public int width() { return width; }
    /** Grid height (y-dimension). */
    public int height() { return height; }
    /** Grid depth (z-dimension, longitudinal). */
    public int depth() { return depth; }
}
