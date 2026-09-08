package vascular.analysis;

import vascular.model.*;

import java.util.ArrayList;
import java.util.List;

/**
 * Analyzes Turing patterns in simulation results.
 * Detects vascular bundles (regions of high [H]), computes their sizes and spatial statistics.
 */
public class PatternAnalyzer {

    private final int bundleThreshold; // minimum H count to qualify as part of a bundle

    public PatternAnalyzer(int bundleThreshold) {
        this.bundleThreshold = bundleThreshold;
    }

    /**
     * Detects connected regions of high [H] (potential vascular bundles)
     * in a 2D cross-section at z=0.
     *
     * @return list of bundles, each described by [centerX, centerY, area, maxConcentration]
     */
    public List<BundleInfo> detectBundles2D(Lattice lattice, int zSlice) {
        DomainMask domain = lattice.domain();
        int[] dims = domain.dimensions();
        int width = dims[0];
        int height = dims[1];

        boolean[][] visited = new boolean[height][width];
        List<BundleInfo> bundles = new ArrayList<>();

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                if (!domain.isInside(x, y, zSlice) || visited[y][x]) continue;
                int count = lattice.get(Species.H, x, y, zSlice);
                if (count < bundleThreshold) continue;

                // Flood-fill to find connected region
                BundleInfo bundle = floodFill2D(lattice, domain, visited, x, y, zSlice);
                if (bundle != null) {
                    bundles.add(bundle);
                }
            }
        }

        return bundles;
    }

    private BundleInfo floodFill2D(Lattice lattice, DomainMask domain, boolean[][] visited,
                                   int startX, int startY, int z) {
        int[] dims = domain.dimensions();
        int width = dims[0];
        int height = dims[1];

        List<int[]> stack = new ArrayList<>();
        stack.add(new int[]{startX, startY});

        int area = 0;
        int maxConcentration = 0;
        double sumX = 0, sumY = 0;

        int[][] dirs = {{1,0}, {-1,0}, {0,1}, {0,-1}};

        while (!stack.isEmpty()) {
            int[] p = stack.removeLast();
            int x = p[0], y = p[1];

            if (x < 0 || x >= width || y < 0 || y >= height) continue;
            if (!domain.isInside(x, y, z)) continue;
            if (visited[y][x]) continue;

            int count = lattice.get(Species.H, x, y, z);
            if (count < bundleThreshold) continue;

            visited[y][x] = true;
            area++;
            sumX += x;
            sumY += y;
            if (count > maxConcentration) maxConcentration = count;

            for (int[] dir : dirs) {
                stack.add(new int[]{x + dir[0], y + dir[1]});
            }
        }

        if (area == 0) return null;

        double cx = sumX / area;
        double cy = sumY / area;
        return new BundleInfo(cx, cy, area, maxConcentration);
    }

    /**
     * Estimates the characteristic wavelength of the Turing pattern using
     * the average nearest-neighbor distance between bundle centers.
     */
    public double estimateWavelength(List<BundleInfo> bundles) {
        if (bundles.size() < 2) return Double.NaN;

        double sumMinDist = 0;
        int count = 0;

        for (int i = 0; i < bundles.size(); i++) {
            BundleInfo bi = bundles.get(i);
            double minDist = Double.POSITIVE_INFINITY;
            for (int j = 0; j < bundles.size(); j++) {
                if (i == j) continue;
                BundleInfo bj = bundles.get(j);
                double dx = bi.centerX - bj.centerX;
                double dy = bi.centerY - bj.centerY;
                double dist = Math.sqrt(dx * dx + dy * dy);
                if (dist < minDist) minDist = dist;
            }
            if (Double.isFinite(minDist)) {
                sumMinDist += minDist;
                count++;
            }
        }

        return count > 0 ? sumMinDist / count : Double.NaN;
    }

    /**
     * Computes the radial concentration profile of a species.
     * Returns the average concentration at each integer distance from the domain center.
     */
    public double[] radialProfile(Lattice lattice, Species species, int zSlice) {
        DomainMask domain = lattice.domain();
        double[] center = domain.center();
        double cx = center[0], cy = center[1];
        int maxR = (int) Math.ceil(domain.radius());

        double[] sum = new double[maxR + 1];
        int[] count = new int[maxR + 1];

        int[] dims = domain.dimensions();
        for (int y = 0; y < dims[1]; y++) {
            for (int x = 0; x < dims[0]; x++) {
                if (!domain.isInside(x, y, zSlice)) continue;
                double dx = x - cx;
                double dy = y - cy;
                int r = (int) Math.sqrt(dx * dx + dy * dy);
                if (r <= maxR) {
                    sum[r] += lattice.get(species, x, y, zSlice);
                    count[r]++;
                }
            }
        }

        double[] profile = new double[maxR + 1];
        for (int r = 0; r <= maxR; r++) {
            profile[r] = count[r] > 0 ? sum[r] / count[r] : 0;
        }
        return profile;
    }

    /** Information about a detected vascular bundle. */
    public record BundleInfo(double centerX, double centerY, int area, int maxConcentration) {
        @Override
        public String toString() {
            return String.format("Bundle at (%.1f, %.1f) area=%d maxH=%d",
                    centerX, centerY, area, maxConcentration);
        }
    }
}
