%% runVascularModel.m
% Top-level script to simulate the Carteni et al. (2014) primary vascular
% structure model. Runs the full spatial PDE simulation and visualizes results.
%
% Reference:
%   Carteni F, Giannino F, Schweingruber FH, Mazzoleni S (2014)
%   Modelling the development and arrangement of the primary vascular
%   structure in plants. Annals of Botany 114(4): 619-627.

clear; close all;
addpath(fileparts(mfilename('fullpath')));

%% ---- Configuration ----
% Domain
N  = 40;   r = 18;   h = 2*r/N;   h2 = h*h;

% Module 1: Morphogen S0 gradient
p.sigma_S0 = 0.012;  p.mu_S0 = 0.015;  p.D_S0 = 0.8;
p.S0_star  = 0.35;

% Module 2: Procambium (Turing activator-substrate)
p.sigma_S  = 0.1;    p.k_S    = 5;      p.rho_S  = 0.3;   p.D_S  = 1.0;
p.sigma_AP = 0.001;  p.rho_AP = 0.5;    p.mu_AP  = 0.015; p.D_AP = 0.005;
p.AP_star  = 0.8;

% Module 3: Phloem/Xylem toggle-switch
p.sigma_AF = 0.0005; p.rho_AF = 0.5;   p.mu_AF  = 0.015; p.D_AF = 0.0002;
p.sigma_AX = 0.0005; p.rho_AX = 0.5;   p.mu_AX  = 0.015; p.D_AX = 0.0002;
p.k_R      = 1;      p.rho_R  = 0.02;  p.mu_R   = 0.02;  p.D_R  = 0.001;
p.K_AF     = 1.5;    p.K_AX   = 1.5;

% Simulation
tSpan = [0, 600];

%% ---- Build spatial domain ----
fprintf('========================================\n');
fprintf('Carteni et al. (2014) Vascular Model\n');
fprintf('========================================\n');

[X, Y] = meshgrid(linspace(-r, r, N), linspace(-r, r, N));
Dmat = sqrt(X.^2 + Y.^2);
mask = Dmat <= r;
nTot  = N * N;
nActive = sum(mask(:));
maskVec = mask(:);
fprintf('Grid: %dx%d (%d active cells, %d ODEs)\n', N, N, nActive, 7*nTot);
fprintf('Time span: [%.0f, %.0f]\n', tSpan(1), tSpan(2));

% S0 steady-state gradient
S0_ss = (p.sigma_S0 / p.mu_S0) * (Dmat / r);
central = S0_ss < p.S0_star;
peripheral = ~central;
fprintf('Central zone: %d cells, Peripheral: %d cells\n', ...
        sum(central(mask)), sum(peripheral(mask)));

%% ---- Build Laplacian ----
[I_vals, J_vals, V_vals] = deal([]);
for i = 1:N
    for j = 1:N
        k = (i-1)*N + j;
        if ~mask(i,j)
            I_vals(end+1)=k; J_vals(end+1)=k; V_vals(end+1)=0;
            continue;
        end
        ndiag = 0;
        if i>1&&mask(i-1,j), nk=(i-2)*N+j; I_vals(end+1)=k; J_vals(end+1)=nk; V_vals(end+1)=1/h2; ndiag=ndiag-1/h2; end
        if i<N&&mask(i+1,j), nk=i*N+j; I_vals(end+1)=k; J_vals(end+1)=nk; V_vals(end+1)=1/h2; ndiag=ndiag-1/h2; end
        if j>1&&mask(i,j-1), nk=(i-1)*N+(j-1); I_vals(end+1)=k; J_vals(end+1)=nk; V_vals(end+1)=1/h2; ndiag=ndiag-1/h2; end
        if j<N&&mask(i,j+1), nk=(i-1)*N+(j+1); I_vals(end+1)=k; J_vals(end+1)=nk; V_vals(end+1)=1/h2; ndiag=ndiag-1/h2; end
        I_vals(end+1)=k; J_vals(end+1)=k; V_vals(end+1)=ndiag;
    end
end
Lap = sparse(I_vals, J_vals, V_vals, nTot, nTot);

%% ---- Initial conditions ----
rng(42);
S0_0 = S0_ss;
S1_0 = 1.0 * central .* mask;
S2_0 = 1.0 * peripheral .* mask;
AP_0 = (0.01 + 0.05 * randn(N, N)) .* mask;
AP_0 = max(AP_0, 0);
AF_0 = zeros(N, N);
AX_0 = zeros(N, N);
R_0  = zeros(N, N);

% Spatial bias for phloem/xylem differentiation
rng(5678);
biasField = 1 + 0.4 * (rand(N, N) - 0.5);
biasField(~mask) = 1;
tmpAF = p.sigma_AF * biasField;  sigmaAFvec = tmpAF(:);
tmpAX = p.sigma_AX * (2 - biasField);  sigmaAXvec = tmpAX(:);

y0 = [S0_0(:); S1_0(:); S2_0(:); AP_0(:); AF_0(:); AX_0(:); R_0(:)];
centralVec = central(:);
periphVec  = peripheral(:);

%% ---- Build JPattern ----
fprintf('Building sparse Jacobian pattern...\n');
nVars = 7;
L_pat = spones(Lap); L_pat(~maskVec,:) = 0;
I_jp = []; J_jp = [];
for rv = 1:nVars
    for cv = 1:nVars
        if rv == cv, pat = L_pat; else, pat = speye(nTot); pat(~maskVec,:) = 0; end
        [pi, pj] = find(pat);
        I_jp = [I_jp; (rv-1)*nTot + pi];
        J_jp = [J_jp; (cv-1)*nTot + pj];
    end
end
JPattern = spones(sparse(I_jp, J_jp, 1, nVars*nTot, nVars*nTot) + speye(nVars*nTot));
fprintf('  %d non-zeros (%.3f%% dense)\n', nnz(JPattern), 100*nnz(JPattern)/(nVars*nTot)^2);

%% ---- Solve ----
fprintf('Solving...\n'); tic;
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'JPattern', JPattern);
[t, y] = ode15s(@(t,y) vascularRHS(t,y, nTot, p, centralVec, periphVec, Lap, maskVec, sigmaAFvec, sigmaAXvec), ...
                tSpan, y0, opts);
fprintf('Done in %.1f sec (%d time points).\n', toc, length(t));

%% ---- Results ----
off = @(k) (k-1)*nTot + (1:nTot);
AP = reshape(y(end, off(4)), [N, N]);
AF = reshape(y(end, off(5)), [N, N]);
AX = reshape(y(end, off(6)), [N, N]);

fprintf('\nAP: Min=%.2f Max=%.2f Median=%.3f Std=%.2f\n', ...
    min(AP(mask)), max(AP(mask)), median(AP(mask)), std(AP(mask)));

% Count AP peaks
np = 0; afp = []; axp = []; ratios = [];
for i = 2:N-1
    for j = 2:N-1
        if mask(i,j) && AP(i,j) > p.AP_star ...
           && AP(i,j) > AP(i-1,j) && AP(i,j) > AP(i+1,j) ...
           && AP(i,j) > AP(i,j-1) && AP(i,j) > AP(i,j+1)
            np = np + 1;
            afp(np) = AF(i,j); axp(np) = AX(i,j);
            ratios(np) = max(afp(np), axp(np)) / max(min(afp(np), axp(np)), 1e-10);
        end
    end
end

fprintf('AP local maxima: %d\n', np);
fprintf('AF range at peaks: [%.2f, %.2f]\n', min(afp), max(afp));
fprintf('AX range at peaks: [%.2f, %.2f]\n', min(axp), max(axp));
fprintf('AF/AX toggle ratio: median=%.1f, >10:1 = %d/%d\n', ...
    median(ratios), sum(ratios > 10), np);

%% ---- Visualize ----
APp = AP; APp(~mask) = NaN;
AFp = AF; AFp(~mask) = NaN;
AXp = AX; AXp(~mask) = NaN;

figure('Position', [50, 50, 1400, 500], 'Color', 'w');
subplot(1,4,1); imagesc(X(1,:), Y(:,1), APp); axis equal tight; colorbar;
set(gca, 'YDir', 'normal'); title(sprintf('A_P (procambium)\\n%d peaks', np));
subplot(1,4,2); imagesc(X(1,:), Y(:,1), AFp); axis equal tight; colorbar;
set(gca, 'YDir', 'normal'); title('A_F (phloem)');
subplot(1,4,3); imagesc(X(1,:), Y(:,1), AXp); axis equal tight; colorbar;
set(gca, 'YDir', 'normal'); title('A_X (xylem)');

overlay = zeros(N, N, 3);
AFn = AF ./ max(AF(mask), [], 'all');
AXn = AX ./ max(AX(mask), [], 'all');
overlay(:,:,1) = AXn; overlay(:,:,2) = AFn;
overlay(repmat(~mask, [1,1,3])) = 1;
subplot(1,4,4); image(X(1,:), Y(:,1), overlay); axis equal tight;
set(gca, 'YDir', 'normal'); title('Overlay: Xylem(red)+Phloem(green)');
sgtitle(sprintf('Carteni et al. (2014) Vascular Pattern (t=%.0f)', t(end)));

fprintf('\nDone. Figures displayed in MATLAB desktop.\n');

%% ---- RHS function ----
function dydt = vascularRHS(~, y, nTot, p, centralVec, periphVec, Lap, maskVec, sigmaAFvec, sigmaAXvec)
    off_n = @(k) (k-1)*nTot + (1:nTot);
    S0_v = y(off_n(1)); S1_v = y(off_n(2)); S2_v = y(off_n(3));
    AP_v = y(off_n(4)); AF_v = y(off_n(5)); AX_v = y(off_n(6)); R_v = y(off_n(7));

    dS0_v = zeros(nTot, 1);  % S0 frozen at steady state

    denom_S = 1 + p.k_S .* AP_v;
    cross_S = p.rho_S .* (AP_v.^2) .* S1_v .* S2_v;
    dS1_v = p.sigma_S .* centralVec .* (1 - S1_v ./ max(denom_S, 1e-6)) ...
            - cross_S + p.D_S .* (Lap * S1_v);
    dS2_v = p.sigma_S .* periphVec .* (1 - S2_v ./ max(denom_S, 1e-6)) ...
            - cross_S + p.D_S .* (Lap * S2_v);
    dAP_v = p.sigma_AP + p.rho_AP .* (AP_v.^2) .* S1_v .* S2_v ...
            - p.mu_AP .* AP_v + p.D_AP .* (Lap * AP_v);

    pcMask = double(AP_v > p.AP_star);
    auto_AF = p.rho_AF .* (AF_v.^2) ./ (p.K_AF^2 + AF_v.^2);
    auto_AX = p.rho_AX .* (AX_v.^2) ./ (p.K_AX^2 + AX_v.^2);
    denom_R = max(1 + p.k_R .* R_v, 1e-6);

    dAF_v = pcMask .* (sigmaAFvec + auto_AF ./ denom_R - p.mu_AF .* AF_v) ...
            + p.D_AF .* (Lap * AF_v);
    dAX_v = pcMask .* (sigmaAXvec + auto_AX ./ denom_R - p.mu_AX .* AX_v) ...
            + p.D_AX .* (Lap * AX_v);
    dR_v  = pcMask .* (p.rho_R .* (AF_v + AX_v) - p.mu_R .* R_v) ...
            + p.D_R .* (Lap * R_v);

    dS0_v(~maskVec)=0; dS1_v(~maskVec)=0; dS2_v(~maskVec)=0;
    dAP_v(~maskVec)=0; dAF_v(~maskVec)=0; dAX_v(~maskVec)=0; dR_v(~maskVec)=0;

    dydt = [dS0_v; dS1_v; dS2_v; dAP_v; dAF_v; dAX_v; dR_v];
end
