package vascular.model;

/**
 * All rate constants and diffusion coefficients for the vascular patterning model.
 * Uses a builder pattern for convenient withXxx() modification.
 *
 * <p>Reaction rates k1–k12 correspond to reactions X1–X12 in Table 1 of Hearn (2019).
 * Diffusion coefficients DH, DB, DM control passive diffusion; DP is always 0.
 * Affinity parameter 'a' controls P-mediated longitudinal H transport strength.
 * Voxel volume V is used to scale reaction propensities (usually 1.0).
 */
public record ModelParameters(
        double k1,   // autocatalysis: 2H + B -> 3H
        double k2,   // B production: ∅ -> B
        double k3,   // B degradation: B -> ∅
        double k4,   // H degradation: H -> ∅
        double k5,   // H production: ∅ -> H
        double k6,   // P production: H -> H + P
        double k7,   // P degradation: P -> ∅
        double k8,   // not used directly (P-mediated transport handled separately)
        double k9,   // M production: ∅ -> M (outer 4% of cylinder only)
        double k10,  // M degradation: M -> ∅
        double k11,  // M + H -> ∅ (annihilation)
        double k12,  // M + B -> ∅ (annihilation)
        double DH,   // H diffusion coefficient
        double DB,   // B diffusion coefficient
        double DP,   // P diffusion coefficient (always 0)
        double DM,   // M diffusion coefficient
        double a,    // H affinity for P (P-mediated transport strength)
        double V     // voxel volume (usually 1.0)
) {

    /** Returns a copy with k1 replaced. */
    public ModelParameters withK1(double v)  { return new ModelParameters(v, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, DH, DB, DP, DM, a, V); }
    /** Returns a copy with k2 replaced. */
    public ModelParameters withK2(double v)  { return new ModelParameters(k1, v, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, DH, DB, DP, DM, a, V); }
    /** Returns a copy with k3 replaced. */
    public ModelParameters withK3(double v)  { return new ModelParameters(k1, k2, v, k4, k5, k6, k7, k8, k9, k10, k11, k12, DH, DB, DP, DM, a, V); }
    /** Returns a copy with k4 replaced. */
    public ModelParameters withK4(double v)  { return new ModelParameters(k1, k2, k3, v, k5, k6, k7, k8, k9, k10, k11, k12, DH, DB, DP, DM, a, V); }
    /** Returns a copy with k5 replaced. */
    public ModelParameters withK5(double v)  { return new ModelParameters(k1, k2, k3, k4, v, k6, k7, k8, k9, k10, k11, k12, DH, DB, DP, DM, a, V); }
    /** Returns a copy with k6 replaced. */
    public ModelParameters withK6(double v)  { return new ModelParameters(k1, k2, k3, k4, k5, v, k7, k8, k9, k10, k11, k12, DH, DB, DP, DM, a, V); }
    /** Returns a copy with k7 replaced. */
    public ModelParameters withK7(double v)  { return new ModelParameters(k1, k2, k3, k4, k5, k6, v, k8, k9, k10, k11, k12, DH, DB, DP, DM, a, V); }
    /** Returns a copy with DH replaced. */
    public ModelParameters withDH(double v)  { return new ModelParameters(k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, v, DB, DP, DM, a, V); }
    /** Returns a copy with DB replaced. */
    public ModelParameters withDB(double v)  { return new ModelParameters(k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, DH, v, DP, DM, a, V); }
    /** Returns a copy with DM replaced. */
    public ModelParameters withDM(double v)  { return new ModelParameters(k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, DH, DB, DP, v, a, V); }
    /** Returns a copy with k9 replaced. */
    public ModelParameters withK9(double v)  { return new ModelParameters(k1, k2, k3, k4, k5, k6, k7, k8, v, k10, k11, k12, DH, DB, DP, DM, a, V); }
    /** Returns a copy with k10 replaced. */
    public ModelParameters withK10(double v) { return new ModelParameters(k1, k2, k3, k4, k5, k6, k7, k8, k9, v, k11, k12, DH, DB, DP, DM, a, V); }
    /** Returns a copy with k11 replaced. */
    public ModelParameters withK11(double v) { return new ModelParameters(k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, v, k12, DH, DB, DP, DM, a, V); }
    /** Returns a copy with k12 replaced. */
    public ModelParameters withK12(double v) { return new ModelParameters(k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, v, DH, DB, DP, DM, a, V); }

    public double DBoverDH() { return DH > 0 ? DB / DH : Double.POSITIVE_INFINITY; }

    /**
     * Returns a new ModelParameters with specified fields overridden.
     * Accepts keys like "k1", "k2", ..., "k12", "DH", "DB", "DP", "DM", "a", "V".
     * Unknown keys are silently ignored.
     */
    public ModelParameters applyOverrides(java.util.Map<String, Double> overrides) {
        if (overrides == null || overrides.isEmpty()) return this;
        ModelParameters p = this;
        for (var entry : overrides.entrySet()) {
            p = switch (entry.getKey().toLowerCase()) {
                case "k1"  -> p.withK1(entry.getValue());
                case "k2"  -> p.withK2(entry.getValue());
                case "k3"  -> p.withK3(entry.getValue());
                case "k4"  -> p.withK4(entry.getValue());
                case "k5"  -> p.withK5(entry.getValue());
                case "k6"  -> p.withK6(entry.getValue());
                case "k7"  -> p.withK7(entry.getValue());
                case "k9"  -> p.withK9(entry.getValue());
                case "k10" -> p.withK10(entry.getValue());
                case "k11" -> p.withK11(entry.getValue());
                case "k12" -> p.withK12(entry.getValue());
                case "dh"  -> p.withDH(entry.getValue());
                case "db"  -> p.withDB(entry.getValue());
                case "dm"  -> p.withDM(entry.getValue());
                default    -> p;
            };
        }
        return p;
    }

    /**
     * Default HB (two-molecule) model parameters.
     * Based on Gray-Scott-Schnakenberg dynamics producing Turing spots.
     */
    public static ModelParameters defaultHB() {
        return new ModelParameters(
                0.001,   // k1: autocatalysis (scaled for stochastic integer counts)
                0.5,     // k2: B production (must dominate autocatalysis consumption)
                0.01,    // k3: B degradation
                0.03,    // k4: H degradation
                0.05,    // k5: H production (basal level, small)
                0.0,     // k6: P production (off)
                0.0,     // k7: P degradation (off)
                0.0,     // k8: (unused)
                0.0,     // k9: M production (off)
                0.0,     // k10: M degradation (off)
                0.0,     // k11: M-H annihilation (off)
                0.0,     // k12: M-B annihilation (off)
                0.05,    // DH: H diffusion
                5.0,     // DB: B diffusion (DB/DH = 100)
                0.0,     // DP: P diffusion (always 0)
                0.0,     // DM: M diffusion (off)
                0.0,     // a: H affinity for P (off)
                1.0      // V: voxel volume
        );
    }

    /**
     * Default HBP (three-molecule) model parameters.
     * Adds P-mediated longitudinal transport.
     */
    public static ModelParameters defaultHBP() {
        return defaultHB()
                .withK6(0.01)   // P production
                .withK7(0.01)   // P degradation
                .withK6(0.005); // use moderate P production
    }

    /**
     * Default HBPM (four-molecule) model parameters.
     * Adds M inhibitor for SVB localization.
     */
    public static ModelParameters defaultHBPM() {
        return new ModelParameters(
                1.0,     // k1
                0.05,    // k2
                0.05,    // k3
                0.02,    // k4
                0.02,    // k5
                0.005,   // k6: P production
                0.01,    // k7: P degradation
                0.0,     // k8: (unused)
                0.005,   // k9: M production (outer 4%)
                0.01,    // k10: M degradation
                1.0,     // k11: M+H annihilation
                1.0,     // k12: M+B annihilation
                0.0375,  // DH
                7.5,     // DB
                0.0,     // DP
                5.0,     // DM
                5.0,     // a: H affinity for P
                1.0      // V
        );
    }

    /** Builder for programmatic construction. */
    public static final class Builder {
        private double k1 = 1.0, k2 = 0.05, k3 = 0.05, k4 = 0.02, k5 = 0.02;
        private double k6, k7, k8, k9, k10, k11, k12;
        private double DH = 0.0375, DB = 7.5, DP, DM, a, V = 1.0;

        public Builder k1(double v)  { k1 = v; return this; }
        public Builder k2(double v)  { k2 = v; return this; }
        public Builder k3(double v)  { k3 = v; return this; }
        public Builder k4(double v)  { k4 = v; return this; }
        public Builder k5(double v)  { k5 = v; return this; }
        public Builder k6(double v)  { k6 = v; return this; }
        public Builder k7(double v)  { k7 = v; return this; }
        public Builder k8(double v)  { k8 = v; return this; }
        public Builder k9(double v)  { k9 = v; return this; }
        public Builder k10(double v) { k10 = v; return this; }
        public Builder k11(double v) { k11 = v; return this; }
        public Builder k12(double v) { k12 = v; return this; }
        public Builder DH(double v)  { DH = v; return this; }
        public Builder DB(double v)  { DB = v; return this; }
        public Builder DP(double v)  { DP = v; return this; }
        public Builder DM(double v)  { DM = v; return this; }
        public Builder a(double v)   { a = v; return this; }
        public Builder V(double v)   { V = v; return this; }

        public ModelParameters build() {
            return new ModelParameters(k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, DH, DB, DP, DM, a, V);
        }
    }
}
