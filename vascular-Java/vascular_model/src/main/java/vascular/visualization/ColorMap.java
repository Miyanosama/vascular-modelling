package vascular.visualization;

import java.awt.Color;

/**
 * Perceptually uniform colormaps for concentration visualization.
 * Uses the Viridis colormap (default) or custom hot-cold schemes.
 */
public enum ColorMap {
    VIRIDIS,
    MAGMA,
    INFERNO,
    PLASMA,
    BLUE_RED;

    /**
     * Maps a normalized value [0, 1] to an RGB color.
     * Values outside [0, 1] are clamped.
     */
    public Color map(double value) {
        double v = Math.max(0.0, Math.min(1.0, value));
        return switch (this) {
            case VIRIDIS -> viridis(v);
            case MAGMA -> magma(v);
            case INFERNO -> inferno(v);
            case PLASMA -> plasma(v);
            case BLUE_RED -> blueRed(v);
        };
    }

    /**
     * Maps a raw molecule count to a color using logarithmic scaling.
     * @param count molecule count
     * @param maxCount maximum count for normalization (maps to 1.0)
     */
    public Color mapLog(int count, int maxCount) {
        if (count <= 0) return Color.BLACK;
        if (maxCount <= 0) return Color.BLACK;
        double v = Math.log1p(count) / Math.log1p(maxCount);
        return map(v);
    }

    /** Maps a raw molecule count to a color using linear scaling. */
    public Color mapLinear(int count, int maxCount) {
        if (count <= 0) return Color.BLACK;
        if (maxCount <= 0) return Color.BLACK;
        return map((double) count / maxCount);
    }

    // ---- Viridis ----
    private static Color viridis(double t) {
        double r = viridisR(t);
        double g = viridisG(t);
        double b = viridisB(t);
        return new Color(clamp(r), clamp(g), clamp(b));
    }
    private static double viridisR(double t) {
        if (t < 0.25) return 0.267;
        if (t < 0.5) return 0.267 + (0.427 - 0.267) * (t - 0.25) / 0.25;
        if (t < 0.75) return 0.427 + (0.746 - 0.427) * (t - 0.5) / 0.25;
        return 0.746 + (0.993 - 0.746) * (t - 0.75) / 0.25;
    }
    private static double viridisG(double t) {
        if (t < 0.25) return 0.004 + (0.286 - 0.004) * t / 0.25;
        if (t < 0.5) return 0.286 + (0.543 - 0.286) * (t - 0.25) / 0.25;
        if (t < 0.75) return 0.543 + (0.741 - 0.543) * (t - 0.5) / 0.25;
        return 0.741 + (0.993 - 0.741) * (t - 0.75) / 0.25;
    }
    private static double viridisB(double t) {
        if (t < 0.25) return 0.329 + (0.565 - 0.329) * t / 0.25;
        if (t < 0.5) return 0.565 + (0.559 - 0.565) * (t - 0.25) / 0.25;
        if (t < 0.75) return 0.559 + (0.383 - 0.559) * (t - 0.5) / 0.25;
        return 0.383 + (0.104 - 0.383) * (t - 0.75) / 0.25;
    }

    // ---- Magma ----
    private static Color magma(double t) {
        double r, g, b;
        if (t < 0.5) {
            double s = t / 0.5;
            r = 0.001 + 0.575 * s;
            g = 0.001 + 0.137 * s;
            b = 0.120 + 0.374 * s;
        } else {
            double s = (t - 0.5) / 0.5;
            r = 0.575 + (0.988 - 0.575) * s;
            g = 0.137 + (0.753 - 0.137) * s;
            b = 0.494 + (0.418 - 0.494) * s;
        }
        return new Color(clamp(r), clamp(g), clamp(b));
    }

    // ---- Inferno ----
    private static Color inferno(double t) {
        double r, g, b;
        if (t < 0.5) {
            double s = t / 0.5;
            r = 0.001 + 0.807 * s;
            g = 0.001 + 0.119 * s;
            b = 0.006 + 0.350 * s;
        } else {
            double s = (t - 0.5) / 0.5;
            r = 0.807 + (0.988 - 0.807) * s;
            g = 0.119 + (0.685 - 0.119) * s;
            b = 0.356 + (0.075 - 0.356) * s;
        }
        return new Color(clamp(r), clamp(g), clamp(b));
    }

    // ---- Plasma ----
    private static Color plasma(double t) {
        double r, g, b;
        if (t < 0.5) {
            double s = t / 0.5;
            r = 0.050 + 0.728 * s;
            g = 0.030 + 0.183 * s;
            b = 0.527 + 0.059 * s;
        } else {
            double s = (t - 0.5) / 0.5;
            r = 0.778 + (0.940 - 0.778) * s;
            g = 0.213 + (0.888 - 0.213) * s;
            b = 0.586 + (0.491 - 0.586) * s;
        }
        return new Color(clamp(r), clamp(g), clamp(b));
    }

    // ---- Blue-Red diverging ----
    private static Color blueRed(double t) {
        if (t < 0.5) {
            double s = t * 2;
            return new Color(clamp(s), clamp(s), clamp(1.0));
        } else {
            double s = (t - 0.5) * 2;
            return new Color(clamp(1.0), clamp(1.0 - s), clamp(1.0 - s));
        }
    }

    private static float clamp(double v) {
        return (float) Math.max(0.0, Math.min(1.0, v));
    }
}
