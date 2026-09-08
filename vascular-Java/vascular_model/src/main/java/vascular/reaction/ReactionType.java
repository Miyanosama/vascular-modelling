package vascular.reaction;

/**
 * The 12 reaction types in the HBPM vascular patterning model.
 * Each corresponds to a row in Table 1 of Hearn (2019).
 */
public enum ReactionType {

    /** R1: 2H + B → 3H (autocatalysis, core Gray-Scott-Schnakenberg) */
    AUTOCATALYSIS(1, true, true, false),

    /** R2: ∅ → B (constitutive B production) */
    B_PRODUCTION(2, false, true, false),

    /** R3: B → ∅ (B degradation) */
    B_DEGRADATION(3, false, true, false),

    /** R4: H → ∅ (H degradation) */
    H_DEGRADATION(4, true, false, false),

    /** R5: ∅ → H (constitutive H production) */
    H_PRODUCTION(5, true, false, false),

    /** R6: H → H + P (P production, induced by H) */
    P_PRODUCTION(6, true, false, true),

    /** R7: P → ∅ (P degradation) */
    P_DEGRADATION(7, false, false, true),

    /** R8: P-mediated longitudinal H transport (handled by diffusion system) */
    P_TRANSPORT(8, true, false, true),

    /** R9: ∅ → M (M production, constrained to outer 4% annulus) */
    M_PRODUCTION(9, false, false, false),

    /** R10: M → ∅ (M degradation) */
    M_DEGRADATION(10, false, false, false),

    /** R11: M + H → ∅ (M-H mutual annihilation) */
    MH_ANNIHILATION(11, true, false, false),

    /** R12: M + B → ∅ (M-B mutual annihilation) */
    MB_ANNIHILATION(12, false, true, false);

    private final int reactionNumber;
    private final boolean affectsH;
    private final boolean affectsB;
    private final boolean affectsP;

    ReactionType(int reactionNumber, boolean affectsH, boolean affectsB, boolean affectsP) {
        this.reactionNumber = reactionNumber;
        this.affectsH = affectsH;
        this.affectsB = affectsB;
        this.affectsP = affectsP;
    }

    public int reactionNumber() { return reactionNumber; }
    public boolean affectsH() { return affectsH; }
    public boolean affectsB() { return affectsB; }
    public boolean affectsP() { return affectsP; }

    /** All reaction types in order X1-X12. */
    public static final ReactionType[] ALL = values();
}
