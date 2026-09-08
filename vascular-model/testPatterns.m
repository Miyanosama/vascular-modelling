%% testPatterns.m — Test script for vascular pattern formation
clear; close all;
addpath(fileparts(mfilename('fullpath')));

%% Domain setup
N = 40;  r = 18;  h = 2*r/N;  h2 = h*h;
[X, Y] = meshgrid(linspace(-r, r, N), linspace(-r, r, N));
Dmat = sqrt(X.^2 + Y.^2);
mask = Dmat <= r;
nTot = N * N;
maskVec = mask(:);

%% S0 steady-state gradient (pre-computed)
sigma_S0 = 0.012;  mu_S0 = 0.015;  D_S0 = 0.8;
S0_ss = (sigma_S0 / mu_S0) * (Dmat / r);
S0_star = 0.35;  % threshold for central vs peripheral
central = S0_ss < S0_star;
peripheral = ~central;
fprintf('Domain: %d total, %d active, %d central, %d peripheral\n', ...
        N*N, sum(maskVec), sum(central(mask)), sum(peripheral(mask)));

%% Build Laplacian
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

%% Module 2 parameters (Turing-tuned)
p.sigma_S  = 0.1;
p.k_S      = 5;
p.rho_S    = 0.3;
p.D_S      = 1.0;
p.sigma_AP = 0.001;
p.rho_AP   = 0.5;
p.mu_AP    = 0.015;
p.D_AP     = 0.005;
p.AP_star  = 0.8;

%% Module 3 parameters — strong autocatalysis, weak repressor for toggle
p.sigma_AF = 0.0005; p.rho_AF = 0.5;  p.mu_AF = 0.015; p.D_AF = 0.0002;
p.sigma_AX = 0.0005; p.rho_AX = 0.5;  p.mu_AX = 0.015; p.D_AX = 0.0002;
p.k_R      = 1;      p.rho_R  = 0.02; p.mu_R  = 0.02;  p.D_R  = 0.001;
p.K_AF     = 1.5;    p.K_AX    = 1.5;

%% Pack into vectors for ODE
centralVec = central(:);
periphVec = peripheral(:);
dovrVec = Dmat(:) / r;

%% Initial conditions
rng(1234);
S0_0 = S0_ss;
S1_0 = 1.0 * central .* mask;
S2_0 = 1.0 * peripheral .* mask;
AP_0 = (0.01 + 0.05 * randn(N, N)) .* mask;
AP_0 = max(AP_0, 0);
% Initial AF/AX = 0 everywhere (only active at procambium via sigma terms)
AF_0 = zeros(N, N);
AX_0 = zeros(N, N);
R_0  = zeros(N, N);

% Create spatial bias field: sigma_AF_bias varies +-20% across domain
% This provides persistent asymmetry for phloem vs xylem differentiation
rng(5678);
biasField = 1 + 0.4 * (rand(N, N) - 0.5);  % range 0.8 to 1.2
biasField(~mask) = 1;
sigma_AF_field = p.sigma_AF * biasField;
sigma_AX_field = p.sigma_AX * (2 - biasField);  % complementary bias

y0 = [S0_0(:); S1_0(:); S2_0(:); AP_0(:); AF_0(:); AX_0(:); R_0(:)];

%% Build JPattern
nVars = 7;
L_pat = spones(Lap);
L_pat(~maskVec, :) = 0;
I_jp = []; J_jp = [];
for rv = 1:nVars
    for cv = 1:nVars
        if rv == cv
            pat = L_pat;
        else
            pat = speye(nTot);
            pat(~maskVec, :) = 0;
        end
        [pi, pj] = find(pat);
        I_jp = [I_jp; (rv-1)*nTot + pi];
        J_jp = [J_jp; (cv-1)*nTot + pj];
    end
end
JPattern = sparse(I_jp, J_jp, 1, nVars*nTot, nVars*nTot);
JPattern = spones(JPattern + speye(nVars*nTot));
fprintf('JPattern: %d nnz (%.3f%% dense)\n', nnz(JPattern), 100*nnz(JPattern)/(nVars*nTot)^2);

%% Solve
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'JPattern', JPattern);
fprintf('Solving...\n'); tic;
sigmaAFvec = sigma_AF_field(:);
sigmaAXvec = sigma_AX_field(:);

[t, y] = ode15s(@(t,y) vascularRHS(t,y, nTot, p, centralVec, periphVec, Lap, maskVec, sigmaAFvec, sigmaAXvec), ...
                [0, 600], y0, opts);
fprintf('Done in %.1f sec (%d points)\n', toc, length(t));

%% Results
off = @(k) (k-1)*nTot + (1:nTot);
AP = reshape(y(end, off(4)), [N, N]);
AF = reshape(y(end, off(5)), [N, N]);
AX = reshape(y(end, off(6)), [N, N]);
S0 = reshape(y(end, off(1)), [N, N]);
S1 = reshape(y(end, off(2)), [N, N]);
S2 = reshape(y(end, off(3)), [N, N]);
R  = reshape(y(end, off(7)), [N, N]);

fprintf('\nAP: Min=%.3f, Max=%.3f, Median=%.3f, Std=%.3f\n', ...
    min(AP(mask)), max(AP(mask)), median(AP(mask)), std(AP(mask)));
fprintf('AP > %.1f: %d cells\n', p.AP_star, sum(AP(mask) > p.AP_star));

% Count local maxima
peaks = 0;
for i = 2:N-1
    for j = 2:N-1
        if mask(i,j) && AP(i,j) > p.AP_star ...
           && AP(i,j) > AP(i-1,j) && AP(i,j) > AP(i+1,j) ...
           && AP(i,j) > AP(i,j-1) && AP(i,j) > AP(i,j+1)
            peaks = peaks + 1;
        end
    end
end
fprintf('AP local maxima: %d\n', peaks);

%% Plot
APp = AP; APp(~mask) = NaN;
AFp = AF; AFp(~mask) = NaN;
AXp = AX; AXp(~mask) = NaN;

figure('Position', [50, 50, 1400, 500], 'Color', 'w');
subplot(1,3,1); imagesc(X(1,:), Y(:,1), APp); axis equal tight; colorbar;
set(gca, 'YDir', 'normal'); title('A_P (procambium)');
subplot(1,3,2); imagesc(X(1,:), Y(:,1), AFp); axis equal tight; colorbar;
set(gca, 'YDir', 'normal'); title('A_F (phloem)');
subplot(1,3,3); imagesc(X(1,:), Y(:,1), AXp); axis equal tight; colorbar;
set(gca, 'YDir', 'normal'); title('A_X (xylem)');
sgtitle(sprintf('Vascular Pattern (t=%.0f, %d local AP maxima)', t(end), peaks));

%% RHS function
function dydt = vascularRHS(~, y, nTot, p, centralVec, periphVec, Lap, maskVec, sigmaAFvec, sigmaAXvec)
    off_n = @(k) (k-1)*nTot + (1:nTot);
    S0_v = y(off_n(1)); S1_v = y(off_n(2)); S2_v = y(off_n(3));
    AP_v = y(off_n(4)); AF_v = y(off_n(5)); AX_v = y(off_n(6)); R_v = y(off_n(7));

    % S0: frozen at steady state
    dS0_v = zeros(nTot, 1);

    % Module 2: Procambium
    denom_S = 1 + p.k_S .* AP_v;
    cross_S = p.rho_S .* (AP_v.^2) .* S1_v .* S2_v;
    dS1_v = p.sigma_S .* centralVec .* (1 - S1_v ./ max(denom_S, 1e-6)) ...
            - cross_S + p.D_S .* (Lap * S1_v);
    dS2_v = p.sigma_S .* periphVec .* (1 - S2_v ./ max(denom_S, 1e-6)) ...
            - cross_S + p.D_S .* (Lap * S2_v);
    dAP_v = p.sigma_AP + p.rho_AP .* (AP_v.^2) .* S1_v .* S2_v ...
            - p.mu_AP .* AP_v + p.D_AP .* (Lap * AP_v);

    % Module 3: Phloem/Xylem toggle-switch
    % Hill-saturated autocatalysis with common repressor R
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

    % Zero outside domain
    dS0_v(~maskVec)=0; dS1_v(~maskVec)=0; dS2_v(~maskVec)=0;
    dAP_v(~maskVec)=0; dAF_v(~maskVec)=0; dAX_v(~maskVec)=0; dR_v(~maskVec)=0;

    dydt = [dS0_v; dS1_v; dS2_v; dAP_v; dAF_v; dAX_v; dR_v];
end
