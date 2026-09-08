function [t, y, grid] = simulateVascularPattern(params)
% SIMULATEVASCULARPATTERN  Simulate the Carteni et al. (2014) vascular model.
%
%   [t, y, grid] = simulateVascularPattern(params)
%
%   Solves the reaction-diffusion PDE system on a 2D circular domain using
%   the method of lines (finite differences + ode15s with sparse Jacobian).

    N  = params.gridSize;
    r  = params.radius;
    h  = params.h;
    h2 = h * h;

    % ---- Build spatial grid ----
    [X, Y] = meshgrid(linspace(-r, r, N), linspace(-r, r, N));
    D = sqrt(X.^2 + Y.^2);
    d_over_r = D / r;
    mask = D <= r;
    nActive = sum(mask(:));

    % Store grid info for output
    grid.X = X;
    grid.Y = Y;
    grid.D = D;
    grid.mask = mask;
    grid.nActive = nActive;

    % ---- State vector layout ----
    % [S0(:); S1(:); S2(:); AP(:); AF(:); AX(:); R(:)], each block N^2 elements
    nTot  = N * N;
    nVars = 7;         % number of species
    nEqs  = nVars * nTot;

    % Block offsets for each species
    off_S0 = 0;  off_S1 = 1;  off_S2 = 2;
    off_AP = 3;  off_AF = 4;  off_AX = 5;  off_R  = 6;

    % Build sparse Laplacian
    Lap = buildLaplacian(N, h2, mask);

    % Build block-diagonal Laplacian for all species
    % L_block = kron(diag(D_values), Lap) but done efficiently
    Dvals = [params.D_S0, params.D_S, params.D_S, ...
             params.D_AP, params.D_AF, params.D_AX, params.D_R];
    L_blocks = cell(nVars, 1);
    for v = 1:nVars
        L_blocks{v} = Dvals(v) * Lap;
    end
    Lfull = blkdiag(L_blocks{:});

    % ---- Build Jacobian sparsity pattern ----
    fprintf('Building sparse Jacobian pattern...\n');
    JPattern = buildJPattern(N, mask, nVars, Lap);
    nnzJ = nnz(JPattern);
    fprintf('  Jacobian: %d x %d with %d non-zeros (%.2f%% dense)\n', ...
            nEqs, nEqs, nnzJ, 100*nnzJ/nEqs^2);

    % ---- Initial conditions ----
    rng(42);
    S0_0 = 0.1 * mask .* (1 - d_over_r);
    S1_0 = params.sigma_S * 0.1 * ones(N, N) .* mask;
    S2_0 = params.sigma_S * 0.1 * ones(N, N) .* mask;
    AP_0 = params.sigma_AP * 0.1 * mask + 0.001 * randn(N, N) .* mask;
    AF_0 = zeros(N, N);
    AX_0 = zeros(N, N);
    R_0  = zeros(N, N);

    y0 = [S0_0(:); S1_0(:); S2_0(:); AP_0(:); AF_0(:); AX_0(:); R_0(:)];

    % ---- Create mask vector for speed ----
    maskVec = mask(:);
    dovrVec = d_over_r(:);

    % ---- Solve ----
    fprintf('Solving on %dx%d grid (%d active cells, %d ODEs)...\n', ...
            N, N, nActive, nEqs);
    tic;

    opts = odeset('RelTol', params.relTol, 'AbsTol', params.absTol, ...
                  'JPattern', JPattern, ...
                  'Vectorized', 'off');

    [t, y] = ode15s(@rhsODE, params.tSpan, y0, opts);

    elapsed = toc;
    fprintf('Simulation complete in %.1f sec (%d time points).\n', ...
            elapsed, length(t));

    % ====================================================================
    % Nested ODE right-hand side (all computed as vectors for speed)
    % ====================================================================
    function dydt = rhsODE(~, y)
        % Unpack
        S0 = y(off_S0*nTot + (1:nTot));
        S1 = y(off_S1*nTot + (1:nTot));
        S2 = y(off_S2*nTot + (1:nTot));
        AP = y(off_AP*nTot + (1:nTot));
        AF = y(off_AF*nTot + (1:nTot));
        AX = y(off_AX*nTot + (1:nTot));
        R  = y(off_R *nTot + (1:nTot));

        % ---- Module 1: Morphogen S0 ----
        dS0 = params.sigma_S0 .* dovrVec - params.mu_S0 .* S0 + Dvals(1) .* (Lap * S0);
        dS0(~maskVec) = 0;

        % ---- Module 2: Procambium ----
        central    = S0 < params.S0_star;
        peripheral = ~central;

        denom_S = 1 + params.k_S .* AP;
        cross_S = params.rho_S .* (AP.^2) .* S1 .* S2;

        dS1 = params.sigma_S .* central .* (1 - S1 ./ denom_S) - cross_S + Dvals(2) .* (Lap * S1);
        dS2 = params.sigma_S .* peripheral .* (1 - S2 ./ denom_S) - cross_S + Dvals(3) .* (Lap * S2);
        dAP = params.sigma_AP + params.rho_AP .* (AP.^2) .* S1 .* S2 ...
              - params.mu_AP .* AP + Dvals(4) .* (Lap * AP);

        dS1(~maskVec) = 0; dS2(~maskVec) = 0; dAP(~maskVec) = 0;

        % ---- Module 3: Phloem/Xylem ----
        pcMask = AP > params.AP_star;
        pcMask_f = double(pcMask);

        denom_R = 1 + params.k_R .* R;
        denom_R(denom_R < 1e-10) = 1e-10;

        dAF = Dvals(5) .* (Lap * AF) ...
              + pcMask_f .* (params.sigma_AF + params.rho_AF .* (AF.^2) ./ denom_R ...
                             - params.mu_AF .* AF);
        dAX = Dvals(6) .* (Lap * AX) ...
              + pcMask_f .* (params.sigma_AX + params.rho_AX .* (AX.^2) ./ denom_R ...
                             - params.mu_AX .* AX);
        dR  = Dvals(7) .* (Lap * R) ...
              + pcMask_f .* (params.rho_R .* (AF + AX) - params.mu_R .* R);

        dAF(~maskVec) = 0; dAX(~maskVec) = 0; dR(~maskVec) = 0;

        % Pack
        dydt = [dS0; dS1; dS2; dAP; dAF; dAX; dR];
    end
end

function Lap = buildLaplacian(N, h2, mask)
% BUILD LAPLACIAN  Sparse 2D Laplacian with zero-flux Neumann BC.
    nTot = N * N;
    [ii, jj] = ndgrid(1:N, 1:N);
    ii = ii(:); jj = jj(:);

    I = []; J = []; V = [];

    for k = 1:nTot
        i = ii(k); j = jj(k);
        if ~mask(i, j)
            I(end+1) = k; J(end+1) = k; V(end+1) = 0; %#ok<*AGROW>
            continue;
        end
        ndiag = 0;
        % North
        if i > 1 && mask(i-1, j)
            nk = (i-2)*N + j;
            I(end+1) = k; J(end+1) = nk; V(end+1) = 1/h2;
            ndiag = ndiag - 1/h2;
        end
        % South
        if i < N && mask(i+1, j)
            nk = i*N + j;
            I(end+1) = k; J(end+1) = nk; V(end+1) = 1/h2;
            ndiag = ndiag - 1/h2;
        end
        % West
        if j > 1 && mask(i, j-1)
            nk = (i-1)*N + (j-1);
            I(end+1) = k; J(end+1) = nk; V(end+1) = 1/h2;
            ndiag = ndiag - 1/h2;
        end
        % East
        if j < N && mask(i, j+1)
            nk = (i-1)*N + (j+1);
            I(end+1) = k; J(end+1) = nk; V(end+1) = 1/h2;
            ndiag = ndiag - 1/h2;
        end
        I(end+1) = k; J(end+1) = k; V(end+1) = ndiag;
    end
    Lap = sparse(I, J, V, nTot, nTot);
end

function Jpat = buildJPattern(N, mask, nVars, Lap)
% BUILDJPATTERN  Sparse Jacobian pattern for ode15s JPattern option.
%
%   Reaction-diffusion system: each species at grid point (i,j) depends on:
%   1. All species at same grid point  (reaction terms)    → nVars x nVars block
%   2. Same species at neighbour points (diffusion terms)  → nVars diag blocks of Lap pattern

    nTot = N * N;
    nEqs = nVars * nTot;

    % Start with block-diagonal: all species interact at same grid point
    [L_i, L_j, ~] = find(Lap);
    L_pattern = sparse(L_i, L_j, 1, nTot, nTot);  % binary pattern of Lap

    % Each species interacts with itself and neighbours spatially
    % PLUS all other species at the SAME grid point

    % Build the composite pattern
    % For each block (r,c) in the nVars x nVars block matrix:
    [I_blocks, J_blocks] = ndgrid(1:nVars, 1:nVars);

    I_all = []; J_all = [];

    for b = 1:numel(I_blocks)
        r = I_blocks(b);  % row species
        c = J_blocks(b);  % col species

        rowOffset = (r-1) * nTot;
        colOffset = (c-1) * nTot;

        if r == c
            % Same species: includes both reaction (diagonal) and diffusion (off-diagonal)
            pat = L_pattern;
        else
            % Different species: only reaction at same grid point (diagonal only)
            pat = speye(nTot);
        end

        % Only keep active (masked) rows
        pat(~mask(:), :) = 0;

        [pi, pj, ~] = find(pat);
        I_all = [I_all; rowOffset + pi]; %#ok<*AGROW>
        J_all = [J_all; colOffset + pj];
    end

    Jpat = sparse(I_all, J_all, 1, nEqs, nEqs);

    % Ensure diagonal entries are present (ode15s requires them)
    diagEntries = (1:nEqs)';
    Jpat = Jpat + sparse(diagEntries, diagEntries, 1, nEqs, nEqs);
    Jpat = spones(Jpat);  % binary
end
