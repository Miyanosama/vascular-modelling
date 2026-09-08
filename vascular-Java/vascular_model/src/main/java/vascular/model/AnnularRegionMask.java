package vascular.model;

/**
 * Identifies voxels in the outer annular region of a cylindrical/circular domain.
 * Used to constrain M (inhibitor) production to the outer 4% of the stem radius,
 * matching the biological hypothesis that inhibitor is produced in the cortex.
 */
public final class AnnularRegionMask {

    private final DomainMask domain;
    private final double innerRadiusFraction; // e.g., 0.96 for outer 4%
    private final double innerRadiusSq;
    private final double radiusSq;
    private final double cx, cy;

    /**
     * @param domain              The underlying domain mask.
     * @param innerRadiusFraction Fraction of domain radius defining the inner boundary
     *                            of the annulus (e.g., 0.96 means outer 4%).
     */
    public AnnularRegionMask(DomainMask domain, double innerRadiusFraction) {
        this.domain = domain;
        this.innerRadiusFraction = innerRadiusFraction;
        double r = domain.radius();
        this.radiusSq = r * r;
        double innerR = r * innerRadiusFraction;
        this.innerRadiusSq = innerR * innerR;
        double[] center = domain.center();
        this.cx = center[0];
        this.cy = center[1];
    }

    /**
     * Returns true if the voxel is inside the annular region
     * (i.e., between innerRadiusFraction*radius and radius).
     */
    public boolean isInOuterShell(int x, int y, int z) {
        if (!domain.isInside(x, y, z)) return false;
        double dx = x - cx;
        double dy = y - cy;
        double distSq = dx * dx + dy * dy;
        return distSq > innerRadiusSq && distSq <= radiusSq;
    }

    /** Returns the fraction of voxels in the outer shell (for diagnostic purposes). */
    public double shellFraction() {
        return 1.0 - innerRadiusFraction * innerRadiusFraction;
    }
}
