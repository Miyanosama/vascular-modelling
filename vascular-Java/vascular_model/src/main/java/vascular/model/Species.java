package vascular.model;

/**
 * Molecular species in the vascular patterning model.
 *
 * H = activator/regulatory molecule (promotes vascular differentiation)
 * B = substrate (consumed during H autocatalysis)
 * P = efflux carrier (facilitates longitudinal H transport)
 * M = inhibitor (suppresses H and B, produced at periphery)
 */
public enum Species {
    H(0),
    B(1),
    P(2),
    M(3);

    private final int index;

    Species(int index) {
        this.index = index;
    }

    /** Zero-based index for array lookups. */
    public int index() {
        return index;
    }

    /** Number of species in the model. */
    public static final int COUNT = values().length;
}
