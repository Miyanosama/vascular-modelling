package vascular.model;

/**
 * 2D circular domain representing a stem cross-section.
 * The circle is centered in the grid with a single z-slice (depth=1).
 */
public final class CircleMask implements DomainMask {

    private final int width, height;
    private final double cx, cy, radius;
    private final int totalVoxels;
    private final int[] linearMap;   // [x + y*width] -> linear index, or -1 if outside
    private final int[] reverseMap;  // linear index -> grid index (y*width + x), O(1) lookup

    /**
     * @param diameter Total grid diameter in voxels (radius = diameter/2).
     */
    public CircleMask(int diameter) {
        this.width = diameter;
        this.height = diameter;
        this.cx = (diameter - 1) / 2.0;
        this.cy = (diameter - 1) / 2.0;
        this.radius = diameter / 2.0;
        this.linearMap = new int[width * height];

        int count = 0;
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                double dx = x - cx;
                double dy = y - cy;
                if (dx * dx + dy * dy <= radius * radius) {
                    linearMap[y * width + x] = count++;
                } else {
                    linearMap[y * width + x] = -1;
                }
            }
        }
        this.totalVoxels = count;

        // Build reverse lookup
        this.reverseMap = new int[totalVoxels];
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int lin = linearMap[y * width + x];
                if (lin >= 0) {
                    reverseMap[lin] = y * width + x;
                }
            }
        }
    }

    @Override
    public boolean isInside(int x, int y, int z) {
        if (z != 0 || x < 0 || x >= width || y < 0 || y >= height) return false;
        double dx = x - cx;
        double dy = y - cy;
        return dx * dx + dy * dy <= radius * radius;
    }

    @Override
    public int totalVoxels() { return totalVoxels; }

    @Override
    public int[] dimensions() { return new int[]{width, height, 1}; }

    @Override
    public double radius() { return radius; }

    @Override
    public double[] center() { return new double[]{cx, cy}; }

    @Override
    public int toLinear(int x, int y, int z) {
        if (z != 0 || x < 0 || x >= width || y < 0 || y >= height) return -1;
        return linearMap[y * width + x];
    }

    @Override
    public VoxelIndex fromLinear(int linear) {
        if (linear < 0 || linear >= totalVoxels) {
            throw new IndexOutOfBoundsException("Linear index " + linear + " out of [0, " + totalVoxels + ")");
        }
        int gridIdx = reverseMap[linear];
        int x = gridIdx % width;
        int y = gridIdx / width;
        return new VoxelIndex(x, y, 0);
    }

    /** Grid width. */
    public int width() { return width; }
    /** Grid height. */
    public int height() { return height; }
}
