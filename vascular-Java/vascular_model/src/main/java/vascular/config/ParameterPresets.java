package vascular.config;

import vascular.model.ModelParameters;
import vascular.model.CircleMask;
import vascular.model.CylinderMask;
import vascular.model.DomainMask;

/**
 * Pre-configured simulation configurations matching the paper's parameter studies.
 */
public final class ParameterPresets {

    private ParameterPresets() {}

    // ---------- 2D HB model configurations ----------

    /** Default HB model: classic Turing spot pattern. */
    public static SimulationConfig hbDefault2D() {
        return new SimulationConfig(
                "HB-Default-2D",
                "HB",
                ModelParameters.defaultHB(),
                new CircleMask(60),
                0,      // maxTime (0 = unlimited)
                30_000_000, // maxEvents (patterns emerge ~15M)
                500_000,    // observerInterval
                42      // seed
        );
    }

    /** 2D HB model with DH=0.016, varying DB/DH ratios. */
    public static SimulationConfig hb2D_DH_low(double dbOverDH) {
        double DH = 0.016;
        return new SimulationConfig(
                "HB-DH=0.016-DB/DH=" + (int) dbOverDH,
                "HB",
                ModelParameters.defaultHB().withDH(DH).withDB(DH * dbOverDH),
                new CircleMask(60),
                0, 50_000_000, 500_000, 42
        );
    }

    /** 2D HB model with DH=0.05, varying DB/DH ratios. */
    public static SimulationConfig hb2D_DH_high(double dbOverDH) {
        double DH = 0.05;
        return new SimulationConfig(
                "HB-DH=0.05-DB/DH=" + (int) dbOverDH,
                "HB",
                ModelParameters.defaultHB().withDH(DH).withDB(DH * dbOverDH),
                new CircleMask(60),
                0, 50_000_000, 500_000, 43
        );
    }

    // ---------- 3D HBP model configurations ----------

    /** Default HBP model with P-mediated longitudinal transport. */
    public static SimulationConfig hbpDefault3D() {
        ModelParameters p = new ModelParameters(
                1.0, 0.05, 0.05, 0.02, 0.02,  // k1-k5
                0.005, 0.01, 0.0,               // k6-k8: P production/degradation
                0.0, 0.0, 0.0, 0.0,            // k9-k12: M off
                0.0375, 7.5, 0.0, 0.0, 5.0, 1.0 // diffusion + affinity
        );
        return new SimulationConfig(
                "HBP-Default-3D",
                "HBP",
                p,
                new CylinderMask(20, 40),
                0, 500_000_000, 10_000_000, 100
        );
    }

    // ---------- 3D HBPM model configurations ----------

    /** Default HBPM model with all 4 molecules. */
    public static SimulationConfig hbpmDefault3D() {
        return new SimulationConfig(
                "HBPM-Default-3D",
                "HBPM",
                ModelParameters.defaultHBPM(),
                new CylinderMask(20, 40),
                0, 500_000_000, 10_000_000, 200
        );
    }

    /** HBPM model with high M production (strong peripheral suppression). */
    public static SimulationConfig hbpmHighM() {
        return new SimulationConfig(
                "HBPM-HighM-3D",
                "HBPM",
                ModelParameters.defaultHBPM().withK9(0.01).withDM(10.0),
                new CylinderMask(20, 40),
                0, 500_000_000, 10_000_000, 201
        );
    }

    /** HBPM model with low M diffusion (sharper SVB zone). */
    public static SimulationConfig hbpmLowDM() {
        return new SimulationConfig(
                "HBPM-LowDM-3D",
                "HBPM",
                ModelParameters.defaultHBPM().withDM(1.0),
                new CylinderMask(20, 40),
                0, 500_000_000, 10_000_000, 202
        );
    }

    // ---------- Size effect studies ----------

    /** Large diameter stem (2D). */
    public static SimulationConfig largeStem2D() {
        return new SimulationConfig(
                "HB-LargeStem-2D",
                "HB",
                ModelParameters.defaultHB(),
                new CircleMask(80),
                0, 100_000_000, 1_000_000, 300
        );
    }

    /** Medium diameter stem (2D). */
    public static SimulationConfig mediumStem2D() {
        return new SimulationConfig(
                "HB-MediumStem-2D",
                "HB",
                ModelParameters.defaultHB(),
                new CircleMask(60), // 75% of large
                0, 100_000_000, 1_000_000, 301
        );
    }

    /** Small diameter stem (2D) — SVBs suppressed. */
    public static SimulationConfig smallStem2D() {
        return new SimulationConfig(
                "HB-SmallStem-2D",
                "HB",
                ModelParameters.defaultHB(),
                new CircleMask(30), // 37% of large
                0, 100_000_000, 1_000_000, 302
        );
    }

    // ---------- 3D diffusion rate studies ----------

    public static SimulationConfig hbpm3D_DH(double DH, double dbOverDH) {
        return new SimulationConfig(
                "HBPM-DH=" + DH + "-DB/DH=" + (int) dbOverDH,
                "HBPM",
                ModelParameters.defaultHBPM().withDH(DH).withDB(DH * dbOverDH),
                new CylinderMask(20, 40),
                0, 500_000_000, 10_000_000, 400
        );
    }

    /**
     * Returns all preset simulation configurations.
     */
    public static SimulationConfig[] allPresets() {
        return new SimulationConfig[]{
                hbDefault2D(),
                hb2D_DH_low(1000),
                hb2D_DH_high(100),
                hbpDefault3D(),
                hbpmDefault3D(),
                largeStem2D(),
                mediumStem2D(),
                smallStem2D(),
        };
    }
}
