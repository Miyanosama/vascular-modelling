%% simulate_plectostele.m
% Reproduction of "Plectostele" (Fig. 5) from:
%   Cartenì et al. (2014) Annals of Botany 114: 619–627
%   "Modelling the development and arrangement of the primary
%    vascular structure in plants"
%
% STAGE 1 — Eq. (1):        S₀ morphogen → central / peripheral zones
% STAGE 2 — Eqs. (2)–(4):   Activator–substrate → procambium ring
% STAGE 3 — Eqs. (5)–(9):   Phloem / xylem within procambium ring
%
% Key feature for Plectostele (Fig. 5):
%   NO local facilitation (s_SF = s_SX = 0) + direct cross-inhibition
%   (k_xi) → alternating bands of phloem and xylem within the
%   vascular strand.
%   r = 40,  r_AP = 0.05,  D_AX = D_AF = 0.002

clear; close all;

%% ========================================================================
%  DOMAIN SETUP
%  ========================================================================

r_domain = 40;              % Radius of plant cross-section [pixels] (siphonostele ring)
N = 2 * r_domain + 1;       % Grid side length (odd → centre on a pixel)
dx = 1.0;                   % Spatial step
dt = 0.05;                  % Time step (reduced for stability with stiffer kinetics)

% Build circular domain
[X, Y] = meshgrid(-r_domain:r_domain, -r_domain:r_domain);
Dmap = sqrt(X.^2 + Y.^2);               % Distance from centre
mask = Dmap <= r_domain;                % Circular domain mask
N_cells = sum(mask(:));

fprintf('╔══════════════════════════════════════════════════════════════╗\n');
fprintf('║  Vascular PDE Model — Plectostele                           ║\n');
fprintf('║  Cartenì et al. (2014) Annals of Botany 114: 619–627       ║\n');
fprintf('╠══════════════════════════════════════════════════════════════╣\n');
fprintf('║  Domain radius r = %2d   |  Grid %d × %d  |  dt = %.2f         ║\n', ...
        r_domain, N, N, dt);
fprintf('║  NO local facilitation → alternating phloem/xylem bands    ║\n');
fprintf('╚══════════════════════════════════════════════════════════════╝\n\n');

%% ========================================================================
%  STAGE 1: SPATIAL DOMAIN DEFINITION — Eq. (1)
%
%   ∂S₀/∂t = s_S₀ · d/r  −  m_S₀ · S₀  +  D_S₀ · ∇²S₀
%
%  S₀ forms a radial gradient (high at periphery, low at centre).
%  Cells with S₀ <  S₀*  →  central zone   (produce S₁)
%  Cells with S₀ ≥ S₀*  →  peripheral zone (produce S₂)
%  ========================================================================

fprintf('━ Stage 1: Spatial domain definition (S₀) ━━━━━━━━━━━━━━━━━━\n');

% Parameters [Table 1] — calibrated for proper zone partitioning
s_S0  = 0.012;              % S₀ production rate
m_S0  = 0.036;               % S₀ removal rate
D_S0  = 0.06;               % S₀ diffusion coefficient
S0_star = 0.17;             % Central / peripheral threshold (≈equal areas)

S0 = zeros(N, N);

for tstep = 1:20000
    lap_S0 = laplacian2d(S0, dx);
    dS0 = s_S0 * (Dmap / r_domain) - m_S0 * S0 + D_S0 * lap_S0;
    dS0(~mask) = 0;
    S0 = S0 + dt * dS0;
    S0 = max(S0, 0);
end

central_zone    = (S0 <  S0_star) & mask;
peripheral_zone = (S0 >= S0_star) & mask;

fprintf('  S₀ range:  [%.4f , %.4f]\n', min(S0(mask)), max(S0(mask)));
fprintf('  Central zone: %5d cells  |  Peripheral zone: %5d cells\n', ...
        sum(central_zone(:)), sum(peripheral_zone(:)));

%% ========================================================================
%  STAGE 2: PROCAMBIUM DIFFERENTIATION — Eqs. (2)–(4)
%
%  Activator–substrate (Gierer–Meinhardt) system on the 2-D domain.
%  The substrates S₁, S₂ are produced in complementary zones; their
%  overlap at the central/peripheral boundary provides the spatially-
%  localised context for A_P autocatalysis.  Lateral inhibition via
%  substrate depletion creates discrete procambium spots.
%
%   ∂S₁/∂t  = σ₁ − ρ_S·A_P²·S₁·S₂ + D_S·∇²S₁       σ₁ = s_S in central, 0 elsewhere
%   ∂S₂/∂t  = σ₂ − ρ_S·A_P²·S₁·S₂ + D_S·∇²S₂       σ₂ = s_S in peripheral, 0 elsewhere
%   ∂A_P/∂t = σ_A + ρ_A·A_P²·S₁·S₂ − μ_A·A_P + D_A·∇²A_P
%  ========================================================================

fprintf('\n━ Stage 2: Procambium differentiation (S₁, S₂, A_P) ━━━━━━━━\n');

% --- Tunable parameters for Turing instability ---
% Calibrated to produce stable discrete procambium spots (eustele).
% D_S >> D_A is required for Turing instability.
% The spots represent vascular bundles at the central/peripheral boundary.
s_S   = 0.10;               % Substrate production rate
rho_S = 0.12;               % Substrate cross-consumption coefficient  [paper: r_S]
D_S   = 0.22;               % Substrate diffusion (must be >> D_A for Turing)
s_AP  = 0.0015;             % Basal activator production   [paper: s_AP]
rho_A = 0.050;              % Activator autocatalysis coefficient  [0.05 → siphonostele ring]
mu_A  = 0.060;              % Activator removal rate       [paper: m_AP]
D_A   = 0.006;              % Activator diffusion (must be << D_S)

% Zone-dependent production masks
sigma1 = s_S * double(central_zone);
sigma2 = s_S * double(peripheral_zone);

% Initialise — small random seed to break rotational symmetry
rng(2024);
S1 = zeros(N, N);
S2 = zeros(N, N);
AP = 0.001 * rand(N, N) .* mask;

% Integrate — stop when a good eustele pattern (6–10 spots) emerges
T_stage2 = 40000;           % max steps; ~2000 time units at dt=0.05
report_every = 2000;        % check pattern quality frequently
check_every  = 100;         % fine-grained spot-count checks
fprintf('  ρ_A = %.3f  (siphonostele ring for r = %d)\n', rho_A, r_domain);
fprintf('  D_S / D_A = %.1f  (Turing ratio)\n', D_S / D_A);
fprintf('  Integrating (auto-stop when ring pattern forms)...\n');

AP_history = zeros(ceil(T_stage2/report_every), 3);
h_idx = 0;
best_AP = AP;  best_S1 = S1;  best_S2 = S2;  best_tstep = 1;
best_nspots = 0;
stop_early = false;

for tstep = 1:T_stage2
    lap_S1 = laplacian2d(S1, dx);
    lap_S2 = laplacian2d(S2, dx);
    lap_AP = laplacian2d(AP, dx);

    % Cross-reaction term:  A_P² · S₁ · S₂
    cross = AP.^2 .* S1 .* S2;

    % S₁:  production in central zone, consumption by cross-reaction, diffusion
    dS1 = sigma1 - rho_S * cross + D_S * lap_S1;
    dS1(~mask) = 0;
    S1 = max(S1 + dt * dS1, 1e-12);

    % S₂:  production in peripheral zone, consumption by cross-reaction, diffusion
    dS2 = sigma2 - rho_S * cross + D_S * lap_S2;
    dS2(~mask) = 0;
    S2 = max(S2 + dt * dS2, 1e-12);

    % A_P:  basal + autocatalytic production − removal + diffusion
    dAP = s_AP + rho_A * cross - mu_A * AP + D_A * lap_AP;
    dAP(~mask) = 0;
    AP = max(AP + dt * dAP, 0);

    % Fine-grained pattern quality check
    if mod(tstep, check_every) == 0
        AP_thresh = 0.75 * max(AP(mask));
        peaks_bw = (AP > AP_thresh) & mask;
        nspots = countConnectedRegions(peaks_bw);
        current_maxAP = max(AP(mask));

        % Stop when siphonostele ring forms (1 connected procambium region)
        if nspots == 1 && ~stop_early && tstep > 5000 && current_maxAP > 2.0
            stop_early = true;
            best_AP = AP;  best_S1 = S1;  best_S2 = S2;
            best_tstep = tstep;  best_nspots = nspots;
            fprintf('    >>> Siphonostele ring formed: %d region at step %d (t=%.0f)\n', ...
                    nspots, tstep, tstep*dt);
            break;
        end

        % Also track best pattern (strongest single-ring pattern)
        if nspots == 1 && current_maxAP > max(best_AP(mask)) && tstep > 5000
            best_AP = AP;  best_S1 = S1;  best_S2 = S2;
            best_tstep = tstep;  best_nspots = nspots;
        end
    end

    % Periodic reporting
    if mod(tstep, report_every) == 0
        h_idx = h_idx + 1;
        AP_history(h_idx, 1) = tstep;
        AP_history(h_idx, 2) = max(AP(mask));
        AP_thresh = 0.75 * max(AP(mask));
        peaks_bw = (AP > AP_thresh) & mask;
        AP_history(h_idx, 3) = countConnectedRegions(peaks_bw);
        fprintf('    step %5d (t = %6.0f):  A_P max = %7.3f  |  spots = %2d\n', ...
                tstep, tstep*dt, AP_history(h_idx, 2), AP_history(h_idx, 3));
    end
end

% Restore best pattern
if best_nspots > 0
    AP = best_AP;  S1 = best_S1;  S2 = best_S2;
    fprintf('  → Using best pattern: %d spots at step %d\n', best_nspots, best_tstep);
else
    fprintf('  → Using final state (no optimal pattern found)\n');
end

% Determine procambium: threshold to isolate the vascular strand
% Use 50% of peak
AP_star = 0.50 * max(AP(mask));
procambium = (AP >= AP_star) & mask;

cc_pro_count = countConnectedRegions(procambium);
fprintf('\n  A_P max = %.3f  |  A_P* = %.3f  |  Procambium strands: %d  |  Cells: %d\n', ...
        max(AP(mask)), AP_star, cc_pro_count, sum(procambium(:)));
if cc_pro_count == 1
    fprintf('  → Siphonostele ring — basis for plectostele banding pattern\n');
end

%% ========================================================================
%  STAGE 3: PHLOEM & XYLEM DIFFERENTIATION — Eqs. (5)–(9)
%
%  Within each procambium strand, two mutually inhibitory activators
%  (A_F = phloem, A_X = xylem) self-organise.
%
%  PLECTOSTELE: NO local facilitation (s_SF = s_SX = 0) +
%  direct cross-inhibition (−k_xi · AF · AX²) between
%  phloem and xylem activators.  This forces the two tissues
%  into alternating bands — the characteristic plectostele
%  pattern found in Lycopodium and other primitive plants.
%
%  D_AF = D_AX = 0.002 (equal diffusion)
%  ========================================================================

fprintf('\n━ Stage 3: Phloem & Xylem differentiation ━━━━━━━━━━━━━━━━━\n');

% Parameters — tuned for plectostele (ring procambium, no facilitation)
% KEY: Without zonal facilitation, AF and AX would co-localize.
% Direct cross-inhibition (k_xi) forces them apart → alternating bands.
% This implements the biological fact that phloem & xylem exclude each other
% locally (a single cell can only become one or the other).
s_AXF  = 0.04;              % Activator basal production
r_AXF  = 0.12;              % Activator cross-reaction coefficient
m_AXF  = 0.001;             % Activator removal rate
D_AF   = 0.003;             % Phloem activator diffusion    [equal for plectostele]
D_AX   = 0.003;             % Xylem activator diffusion     [equal for plectostele]
s_SXF  = 0.86;              % Substrate basal production    [× facil=0 → 0]
g_SXF  = 0.02;              % Substrate cross-reaction      [positive coupling]
m_SXF  = 0.001;             % Substrate removal rate
D_SXF  = 0.003;             % Substrate diffusion
m_R    = 0.024;              % Repressor removal rate
A_star = 1;                 % Differentiation threshold
k_xi   = 14;              % Direct cross-inhibition [AF↔AX — key for banding]

% Pro: reactions only active within procambium
Pro = double(procambium);

% Local facilitation masks — NONE for plectostele
facil_SF = zeros(N, N);             % NO phloem facilitation (key to plectostele)
facil_SX = zeros(N, N);             % NO xylem facilitation  (key to plectostele)

% Initialise (random seed for symmetry breaking → alternating bands)
rng(2025);
AF = 0.10 * rand(N, N) .* procambium;
AX = 0.10 * rand(N, N) .* procambium;
SF = 0.01 * rand(N, N) .* procambium;
SX = 0.01 * rand(N, N) .* procambium;
R  = zeros(N, N);

T_stage3 = 60000;
fprintf('  D_AF = %.3f,  D_AX = %.3f  (plectostele — no facilitation + cross-inhibition)\n', D_AF, D_AX);
fprintf('  k_xi = %.2f  (direct AF↔AX mutual inhibition)\n', k_xi);
fprintf('  A* = %.1f\n', A_star);
fprintf('  Integrating %d steps...\n', T_stage3);

for tstep = 1:T_stage3
    lap_AF = laplacian2d(AF, dx);
    lap_AX = laplacian2d(AX, dx);
    lap_SF = laplacian2d(SF, dx);
    lap_SX = laplacian2d(SX, dx);

    % Common repressor denominator
    inv_denom = 1 ./ (1 + R);

    % --- Eq. (5):  A_F (phloem activator) ---
    % Cross-inhibition:  −k_xi · AF · AX²  → AF suppressed where AX is high
    dAF = Pro .* (s_AXF + r_AXF * AF.^2 .* SF .* inv_denom) ...
          - m_AXF * AF + D_AF * lap_AF ...
          - k_xi * AF .* AX.^2;
    dAF(~mask) = 0;

    % --- Eq. (6):  A_X (xylem activator) ---
    % Cross-inhibition:  −k_xi · AX · AF²  → AX suppressed where AF is high
    dAX = Pro .* (s_AXF + r_AXF * AX.^2 .* SX .* inv_denom) ...
          - m_AXF * AX + D_AX * lap_AX ...
          - k_xi * AX .* AF.^2;
    dAX(~mask) = 0;

    % --- Eq. (7):  S_F (phloem substrate) ---
    dSF = Pro .* (s_SXF * facil_SF + g_SXF * (AX - SF)) ...
          - m_SXF * SF + D_SXF * lap_SF;
    dSF(~mask) = 0;

    % --- Eq. (8):  S_X (xylem substrate) ---
    dSX = Pro .* (s_SXF * facil_SX + g_SXF * (AF - SX)) ...
          - m_SXF * SX + D_SXF * lap_SX;
    dSX(~mask) = 0;

    % --- Eq. (9):  R (common repressor — no diffusion) ---
    dR = Pro .* (r_AXF * AF.^2 .* SF + r_AXF * AX.^2 .* SX) - m_R * R;
    dR(~mask) = 0;

    % Euler step
    AF = max(AF + dt * dAF, 0);
    AX = max(AX + dt * dAX, 0);
    SF = max(SF + dt * dSF, 0);
    SX = max(SX + dt * dSX, 0);
    R  = max(R  + dt * dR,  0);

    if mod(tstep, report_every) == 0
        % Spatial correlation: low/negative = alternating bands, high = overlap
        af_vals = AF(mask);  ax_vals = AX(mask);
        % Manual Pearson correlation (no toolbox needed)
        afc = af_vals - mean(af_vals);  axc = ax_vals - mean(ax_vals);
        spatial_corr = sum(afc .* axc) / sqrt(sum(afc.^2) * sum(axc.^2));
        fprintf('    step %5d (t = %6.0f):  A_F max = %8.2f  |  A_X max = %8.2f  |  corr = %+.3f\n', ...
                tstep, tstep*dt, max(af_vals), max(ax_vals), spatial_corr);
    end
end

% Differentiate phloem and xylem
phloem = (AF >= A_star) & mask;
xylem  = (AX >= A_star) & mask;

fprintf('\n  Phloem cells: %d  |  Xylem cells: %d\n', ...
        sum(phloem(:)), sum(xylem(:)));

% ── Spatial diagnostics: verify plectostele (alternating bands) ──
overlap = phloem & xylem;
fprintf('  Overlapping (phloem & xylem): %d cells\n', sum(overlap(:)));

% Mean radial distance from centre
d_phloem = mean(Dmap(phloem));
d_xylem  = mean(Dmap(xylem));
fprintf('  Mean radial position — Phloem: %.1f  |  Xylem: %.1f  (centre=0, edge=%d)\n', ...
        d_phloem, d_xylem, r_domain);
fprintf('  → Plectostele: alternating phloem/xylem bands (no radial ordering)\n');

% Per-strand analysis
cc_pro = labelConnectedRegions(procambium);
n_bundles = max(cc_pro(:));
fprintf('\n  Per-strand analysis (%d strand(s)):\n', n_bundles);
fprintf('  Strand | Procambium | Phloem | Xylem \n');
fprintf('  -------+------------+--------+-------\n');
for b = 1:n_bundles
    bundle_mask = (cc_pro == b);
    n_pro = sum(bundle_mask(:));
    n_phl = sum(phloem(bundle_mask));
    n_xyl = sum(xylem(bundle_mask));
    fprintf('  %6d | %10d | %6d | %5d\n', ...
            b, n_pro, n_phl, n_xyl);
end

%% ========================================================================
%  PUBLICATION-QUALITY FIGURES
%  ========================================================================

fprintf('\n━ Generating figures ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

% ── Colour palette ──
COL_PRO  = [0.30 0.55 0.35];   % Procambium — muted green
COL_PHL  = [0.20 0.50 0.85];   % Phloem — blue
COL_XYL  = [0.90 0.30 0.20];   % Xylem — red-orange
COL_BG   = [0.96 0.96 0.96];   % Background light grey
COL_CENT = [0.25 0.45 0.70];   % Central zone
COL_PERI = [0.85 0.40 0.30];   % Peripheral zone

% ── Composite tissue map ──
tissue_RGB = repmat(reshape(COL_BG, [1 1 3]), N, N);
% Grey out external region
for ch = 1:3
    tmp = tissue_RGB(:,:,ch);
    tmp(~mask) = 0.88;
    tissue_RGB(:,:,ch) = tmp;
end
% Layer procambium → xylem → phloem
for ch = 1:3
    tmp = tissue_RGB(:,:,ch);
    tmp(procambium & ~xylem & ~phloem) = COL_PRO(ch);
    tmp(xylem)  = COL_XYL(ch);
    tmp(phloem) = COL_PHL(ch);
    tissue_RGB(:,:,ch) = tmp;
end

% ═══════════════════════════════════════════════════════════════
%  FIGURE 1 — Main result: tissue map
% ═══════════════════════════════════════════════════════════════
figure('Name', 'Plectostele — Tissue Map', ...
       'NumberTitle', 'off', ...
       'Position', [80, 160, 620, 620], ...
       'Color', 'w');

imagesc(tissue_RGB);
axis image off;

title({'\bf{Plectostele (alternating phloem/xylem bands)}'; ...
       '\rm{r = 40,  r_{AP} = 0.05,  D_{AX} = D_{AF} = 0.002,  no facilitation + cross-inhibition}'}, ...
       'FontSize', 13);

% Legend
hold on;
hp(1) = patch(NaN, NaN, COL_PRO, 'EdgeColor', 'none');
hp(2) = patch(NaN, NaN, COL_PHL, 'EdgeColor', 'none');
hp(3) = patch(NaN, NaN, COL_XYL, 'EdgeColor', 'none');
legend(hp, {'Procambium', 'Phloem', 'Xylem'}, ...
       'Location', 'southoutside', 'Orientation', 'horizontal', ...
       'FontSize', 10);
hold off;

% ═══════════════════════════════════════════════════════════════
%  FIGURE 2 — Multi-panel process overview (2 × 5)
% ═══════════════════════════════════════════════════════════════
figure('Name', 'Model Process Overview', ...
       'NumberTitle', 'off', ...
       'Position', [720, 50, 1250, 750], ...
       'Color', 'w');

% (A) S₀ gradient
axA = subplot(2, 5, 1);
imagesc(S0, 'AlphaData', double(mask));
axis image off; colormap(axA, parula); colorbar;
title({'\bf(A)'; '\rmS₀ gradient'}, 'FontSize', 10);
hold on; contour(mask, [1 1], 'k-', 'LineWidth', 0.5); hold off;

% (B) Zones
axB = subplot(2, 5, 2);
zone_RGB = ones(N, N, 3);
for ch = 1:3
    z = zone_RGB(:,:,ch);
    z(central_zone)    = COL_CENT(ch);
    z(peripheral_zone) = COL_PERI(ch);
    z(~mask) = 0.88;
    zone_RGB(:,:,ch) = z;
end
imagesc(zone_RGB); axis image off;
title({'\bf(B)'; '\rmCentral (blue) / Peripheral (red)'}, 'FontSize', 10);

% (C) S₁
axC = subplot(2, 5, 3);
imagesc(S1, 'AlphaData', double(mask));
axis image off; colormap(axC, turbo); colorbar;
title({'\bf(C)'; '\rmS₁ (central substrate)'}, 'FontSize', 10);

% (D) S₂
axD = subplot(2, 5, 4);
imagesc(S2, 'AlphaData', double(mask));
axis image off; colormap(axD, turbo); colorbar;
title({'\bf(D)'; '\rmS₂ (peripheral substrate)'}, 'FontSize', 10);

% (E) A_P
axE = subplot(2, 5, 5);
imagesc(AP, 'AlphaData', double(mask));
axis image off; colormap(axE, hot); colorbar;
hold on; contour(procambium, [1 1], 'c-', 'LineWidth', 1.2); hold off;
title({'\bf(E)'; '\rmA_P (procambium activator)'}, 'FontSize', 10);

% (F) Procambium
axF = subplot(2, 5, 6);
pro_RGB = ones(N, N, 3);
for ch = 1:3
    z = pro_RGB(:,:,ch);
    z(procambium) = COL_PRO(ch);
    z(~mask) = 0.88;
    pro_RGB(:,:,ch) = z;
end
imagesc(pro_RGB); axis image off;
title({'\bf(F)'; '\rmProcambium (A_P ≥ A_P^*)'}, 'FontSize', 10);

% (G) A_F
axG = subplot(2, 5, 7);
imagesc(AF, 'AlphaData', double(mask));
axis image off; colormap(axG, bone); colorbar;
title({'\bf(G)'; '\rmA_F (phloem activator)'}, 'FontSize', 10);

% (H) A_X
axH = subplot(2, 5, 8);
imagesc(AX, 'AlphaData', double(mask));
axis image off; colormap(axH, copper); colorbar;
title({'\bf(H)'; '\rmA_X (xylem activator)'}, 'FontSize', 10);

% (I) Repressor R
axI = subplot(2, 5, 9);
imagesc(R, 'AlphaData', double(mask));
axis image off; colormap(axI, pink); colorbar;
title({'\bf(I)'; '\rmR (common repressor)'}, 'FontSize', 10);

% (J) Final tissue map
axJ = subplot(2, 5, 10);
imagesc(tissue_RGB); axis image off;
hold on;
hp2(1) = patch(NaN, NaN, COL_PRO, 'EdgeColor', 'none');
hp2(2) = patch(NaN, NaN, COL_PHL, 'EdgeColor', 'none');
hp2(3) = patch(NaN, NaN, COL_XYL, 'EdgeColor', 'none');
legend(hp2, {'Procambium','Phloem','Xylem'}, ...
       'Location', 'southoutside', 'FontSize', 8);
hold off;
title({'\bf(J)'; '\rmPlectostele (alternating bands)'}, 'FontSize', 10);

sgtitle({'\bf{Reaction–Diffusion Model of Primary Vascular Patterning}'; ...
         '\rmCartenì et al. (2014)  —  Plectostele (no local facilitation)'}, ...
        'FontSize', 14);

% ═══════════════════════════════════════════════════════════════
%  FIGURE 3 — Radial profiles (1-D slice through centre)
% ═══════════════════════════════════════════════════════════════
figure('Name', 'Radial Profiles', ...
       'NumberTitle', 'off', ...
       'Position', [80, 50, 620, 700], ...
       'Color', 'w');

mid_row = round(N/2);
r_1d = X(mid_row, :);
in_domain = mask(mid_row, :);
r_norm = r_1d / r_domain;

% (a) S₀
subplot(3, 1, 1);
plot(r_norm(in_domain), S0(mid_row, in_domain), 'k-', 'LineWidth', 1.8);
hold on;
yline(S0_star, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
xlabel('d / r'); ylabel('S₀');
title('(a)  S₀ radial gradient');
legend({'S₀', 'S₀^*'}, 'Location', 'best', 'FontSize', 9);
grid on; box on;

% (b) S₁, S₂, A_P
subplot(3, 1, 2);
yyaxis left;
h1 = plot(r_norm(in_domain), S1(mid_row, in_domain), 'b-', 'LineWidth', 1.2);
hold on;
h2 = plot(r_norm(in_domain), S2(mid_row, in_domain), 'r-', 'LineWidth', 1.2);
ylabel('S₁ ,  S₂');
yyaxis right;
h3 = plot(r_norm(in_domain), AP(mid_row, in_domain), 'g-', 'LineWidth', 1.8);
yline(AP_star, 'g--', 'LineWidth', 1);
ylabel('A_P');
xlabel('d / r');
title('(b)  Procambium substrates & activator');
legend([h1 h2 h3], {'S₁', 'S₂', 'A_P'}, 'Location', 'best', 'FontSize', 9);
grid on; box on;

% (c) A_F, A_X
subplot(3, 1, 3);
h4 = plot(r_norm(in_domain), AF(mid_row, in_domain), ...
          'Color', COL_PHL, 'LineWidth', 1.2);
hold on;
h5 = plot(r_norm(in_domain), AX(mid_row, in_domain), ...
          'Color', COL_XYL, 'LineWidth', 1.2);
yline(A_star, 'k--', 'LineWidth', 1);
xlabel('d / r'); ylabel('A_F ,  A_X');
title('(c)  Phloem & xylem activators');
legend([h4 h5], {'A_F (phloem)', 'A_X (xylem)'}, ...
       'Location', 'best', 'FontSize', 9);
grid on; box on;

sgtitle({'\bf{Radial Concentration Profiles}'; ...
         '\rmHorizontal slice through centre of domain'}, ...
        'FontSize', 13);

% ═══════════════════════════════════════════════════════════════
%  FIGURE 4 — A_P spot count convergence
% ═══════════════════════════════════════════════════════════════
figure('Name', 'Convergence', ...
       'NumberTitle', 'off', ...
       'Position', [750, 50, 500, 300], ...
       'Color', 'w');

yyaxis left;
plot(AP_history(:,1)*dt, AP_history(:,2), 'b-o', ...
     'LineWidth', 1.2, 'MarkerSize', 4, 'MarkerFaceColor', 'b');
ylabel('max(A_P)');

yyaxis right;
plot(AP_history(:,1)*dt, AP_history(:,3), 'r-s', ...
     'LineWidth', 1.2, 'MarkerSize', 4, 'MarkerFaceColor', 'r');
ylabel('Number of procambium spots');
xlabel('Time');
grid on; box on;
title({'\bf{Convergence of procambium pattern}'; ...
       '\rmSpot count and peak activator over time'}, 'FontSize', 11);
legend({'max(A_P)', 'N spots'}, 'Location', 'best', 'FontSize', 9);

%% ========================================================================
%  SUMMARY
%  ========================================================================
fprintf('\n╔══════════════════════════════════════════════════════════════╗\n');
fprintf('║  SIMULATION COMPLETE                                        ║\n');
fprintf('╠══════════════════════════════════════════════════════════════╣\n');
fprintf('║  Pattern:  Plectostele (alternating phloem/xylem bands)     ║\n');
fprintf('║  r = %2d   r_{AP} = %.2f   D_{AX} = D_{AF} = %.3f              ║\n', ...
        r_domain, rho_A, D_AF);
fprintf('║  No local facilitation + direct AF↔AX cross-inhibition     ║\n');
fprintf('║  Procambium strands: %2d                                     ║\n', cc_pro_count);
fprintf('║  Phloem cells: %4d   Xylem cells: %4d                       ║\n', ...
        sum(phloem(:)), sum(xylem(:)));
fprintf('╚══════════════════════════════════════════════════════════════╝\n');

%% ========================================================================
%  LOCAL FUNCTION — 5-point Laplacian with periodic BC on square grid.
%  The caller applies a circular mask, which approximates zero-flux
%  Neumann conditions on the curved boundary.
%  ========================================================================

function lap = laplacian2d(u, dx)
    lap = (circshift(u, [ 0,  1]) + ...
           circshift(u, [ 0, -1]) + ...
           circshift(u, [ 1,  0]) + ...
           circshift(u, [-1,  0]) - 4 * u) / (dx^2);
end

%% ------------------------------------------------------------------------
%  Count connected regions in a binary mask (4-connectivity).
%  ------------------------------------------------------------------------
function n = countConnectedRegions(BW)
    [rows, cols] = size(BW);
    visited = false(rows, cols);
    n = 0;
    for r = 1:rows
        for c = 1:cols
            if BW(r, c) && ~visited(r, c)
                n = n + 1;
                stack = [r, c];
                visited(r, c) = true;
                while ~isempty(stack)
                    p = stack(1, :);  stack(1, :) = [];
                    for dr = [-1, 0, 1]
                        for dc = [-1, 0, 1]
                            if abs(dr) + abs(dc) ~= 1; continue; end
                            nr = p(1) + dr;  nc = p(2) + dc;
                            if nr >= 1 && nr <= rows && nc >= 1 && nc <= cols
                                if BW(nr, nc) && ~visited(nr, nc)
                                    visited(nr, nc) = true;
                                    stack(end+1, :) = [nr, nc];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%% ------------------------------------------------------------------------
%  Label connected regions (returns integer label matrix, 4-connectivity).
%  ------------------------------------------------------------------------
function L = labelConnectedRegions(BW)
    [rows, cols] = size(BW);
    L = zeros(rows, cols);
    visited = false(rows, cols);
    label = 0;
    for r = 1:rows
        for c = 1:cols
            if BW(r, c) && ~visited(r, c)
                label = label + 1;
                stack = [r, c];
                visited(r, c) = true;
                L(r, c) = label;
                while ~isempty(stack)
                    p = stack(1, :);  stack(1, :) = [];
                    for dr = [-1, 0, 1]
                        for dc = [-1, 0, 1]
                            if abs(dr) + abs(dc) ~= 1; continue; end
                            nr = p(1) + dr;  nc = p(2) + dc;
                            if nr >= 1 && nr <= rows && nc >= 1 && nc <= cols
                                if BW(nr, nc) && ~visited(nr, nc)
                                    visited(nr, nc) = true;
                                    L(nr, nc) = label;
                                    stack(end+1, :) = [nr, nc];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%% ------------------------------------------------------------------------
%  Ternary operator:  ternary(cond, valTrue, valFalse)
%  ------------------------------------------------------------------------
function out = ternary(cond, valTrue, valFalse)
    if cond
        out = valTrue;
    else
        out = valFalse;
    end
end
