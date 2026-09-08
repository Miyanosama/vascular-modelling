package vascular.reaction;

import vascular.model.Lattice;
import vascular.model.ModelParameters;
import vascular.model.Species;

/**
 * Computes reaction propensities for all 12 reactions in the HBPM model.
 *
 * <p>Propensities follow the stochastic formulation from the paper:
 * zero-order reactions scale by voxel volume V,
 * first-order reactions scale by molecule count,
 * second-order reactions use combinatorial factors.
 */
public final class ReactionPropensity {

    private ReactionPropensity() { /* static utility */ }

    /**
     * Computes the propensity for a specific reaction at a specific voxel.
     *
     * @param type   the reaction type
     * @param lattice the current lattice state
     * @param voxelIdx the linear voxel index
     * @param params model parameters (rate constants, volume)
     * @param isOuterShell whether this voxel is in the outer 4% annulus (for M production)
     * @return the reaction propensity (non-negative)
     */
    public static double compute(
            ReactionType type,
            Lattice lattice,
            int voxelIdx,
            ModelParameters params,
            boolean isOuterShell
    ) {
        int H = lattice.get(Species.H, voxelIdx);
        int B = lattice.get(Species.B, voxelIdx);
        int P = lattice.get(Species.P, voxelIdx);
        int M = lattice.get(Species.M, voxelIdx);
        double V = params.V();

        return switch (type) {
            case AUTOCATALYSIS -> {
                // R1: 2H + B → 3H, rate = k1 * H*(H-1)*B / V^2
                if (H >= 2 && B >= 1) {
                    yield params.k1() * H * (H - 1.0) * B / (V * V);
                }
                yield 0.0;
            }
            case B_PRODUCTION -> {
                // R2: ∅ → B, rate = k2 * V
                yield params.k2() * V;
            }
            case B_DEGRADATION -> {
                // R3: B → ∅, rate = k3 * B
                yield params.k3() * B;
            }
            case H_DEGRADATION -> {
                // R4: H → ∅, rate = k4 * H
                yield params.k4() * H;
            }
            case H_PRODUCTION -> {
                // R5: ∅ → H, rate = k5 * V
                yield params.k5() * V;
            }
            case P_PRODUCTION -> {
                // R6: H → H + P, rate = k6 * H
                yield params.k6() * H;
            }
            case P_DEGRADATION -> {
                // R7: P → ∅, rate = k7 * P
                yield params.k7() * P;
            }
            case P_TRANSPORT -> {
                // R8: P-mediated longitudinal H transport
                // Handled by the diffusion system; this is a placeholder for the event pool
                // The actual propensity is computed in DiffusionManager
                yield 0.0;
            }
            case M_PRODUCTION -> {
                // R9: ∅ → M, rate = k9 * V (only in outer 4% annulus)
                yield isOuterShell ? params.k9() * V : 0.0;
            }
            case M_DEGRADATION -> {
                // R10: M → ∅, rate = k10 * M
                yield params.k10() * M;
            }
            case MH_ANNIHILATION -> {
                // R11: M + H → ∅, rate = k11 * M * H / V
                if (M >= 1 && H >= 1) {
                    yield params.k11() * M * H / V;
                }
                yield 0.0;
            }
            case MB_ANNIHILATION -> {
                // R12: M + B → ∅, rate = k12 * M * B / V
                if (M >= 1 && B >= 1) {
                    yield params.k12() * M * B / V;
                }
                yield 0.0;
            }
        };
    }

    /**
     * Computes total reaction propensity for a voxel (sum of all 12 reactions).
     */
    public static double totalReaction(
            Lattice lattice,
            int voxelIdx,
            ModelParameters params,
            boolean isOuterShell
    ) {
        double total = 0.0;
        for (ReactionType type : ReactionType.ALL) {
            total += compute(type, lattice, voxelIdx, params, isOuterShell);
        }
        return total;
    }

    /**
     * Executes a reaction by modifying molecule counts at the given voxel.
     * Returns the delta for each species (for propensity recalculation).
     */
    public static int[] execute(
            ReactionType type,
            Lattice lattice,
            int voxelIdx
    ) {
        return switch (type) {
            case AUTOCATALYSIS -> {
                // 2H + B → 3H: net +1 H, -1 B
                lattice.adjust(Species.H, voxelIdx, 1);
                lattice.adjust(Species.B, voxelIdx, -1);
                yield new int[]{1, -1, 0, 0}; // H, B, P, M deltas
            }
            case B_PRODUCTION -> {
                // ∅ → B: +1 B
                lattice.adjust(Species.B, voxelIdx, 1);
                yield new int[]{0, 1, 0, 0};
            }
            case B_DEGRADATION -> {
                // B → ∅: -1 B
                lattice.adjust(Species.B, voxelIdx, -1);
                yield new int[]{0, -1, 0, 0};
            }
            case H_DEGRADATION -> {
                // H → ∅: -1 H
                lattice.adjust(Species.H, voxelIdx, -1);
                yield new int[]{-1, 0, 0, 0};
            }
            case H_PRODUCTION -> {
                // ∅ → H: +1 H
                lattice.adjust(Species.H, voxelIdx, 1);
                yield new int[]{1, 0, 0, 0};
            }
            case P_PRODUCTION -> {
                // H → H + P: +1 P (H unchanged)
                lattice.adjust(Species.P, voxelIdx, 1);
                yield new int[]{0, 0, 1, 0};
            }
            case P_DEGRADATION -> {
                // P → ∅: -1 P
                lattice.adjust(Species.P, voxelIdx, -1);
                yield new int[]{0, 0, -1, 0};
            }
            case P_TRANSPORT -> {
                // Handled by diffusion system — no direct molecule count change here
                yield new int[]{0, 0, 0, 0};
            }
            case M_PRODUCTION -> {
                // ∅ → M: +1 M
                lattice.adjust(Species.M, voxelIdx, 1);
                yield new int[]{0, 0, 0, 1};
            }
            case M_DEGRADATION -> {
                // M → ∅: -1 M
                lattice.adjust(Species.M, voxelIdx, -1);
                yield new int[]{0, 0, 0, -1};
            }
            case MH_ANNIHILATION -> {
                // M + H → ∅: -1 M, -1 H
                lattice.adjust(Species.M, voxelIdx, -1);
                lattice.adjust(Species.H, voxelIdx, -1);
                yield new int[]{-1, 0, 0, -1};
            }
            case MB_ANNIHILATION -> {
                // M + B → ∅: -1 M, -1 B
                lattice.adjust(Species.M, voxelIdx, -1);
                lattice.adjust(Species.B, voxelIdx, -1);
                yield new int[]{0, -1, 0, -1};
            }
        };
    }

    /** Array of all reaction types excluding P_TRANSPORT (which is handled by diffusion). */
    public static final ReactionType[] REACTION_EVENTS = {
            ReactionType.AUTOCATALYSIS,
            ReactionType.B_PRODUCTION,
            ReactionType.B_DEGRADATION,
            ReactionType.H_DEGRADATION,
            ReactionType.H_PRODUCTION,
            ReactionType.P_PRODUCTION,
            ReactionType.P_DEGRADATION,
            ReactionType.M_PRODUCTION,
            ReactionType.M_DEGRADATION,
            ReactionType.MH_ANNIHILATION,
            ReactionType.MB_ANNIHILATION
    };

    /** Number of non-diffusion reaction types. */
    public static final int REACTION_COUNT = REACTION_EVENTS.length;
}
