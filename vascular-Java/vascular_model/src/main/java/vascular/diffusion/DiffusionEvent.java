package vascular.diffusion;

import vascular.model.Species;

/**
 * Records the result of a diffusion event.
 */
public record DiffusionEvent(
        Species species,
        int sourceVoxel,
        int targetVoxel,
        boolean isLongitudinal
) {
    @Override
    public String toString() {
        String dir = isLongitudinal ? "∥" : "⟂";
        return species + " " + sourceVoxel + "→" + targetVoxel + " [" + dir + "]";
    }
}
