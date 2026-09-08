function model = buildSimbioVascularModel(params)
% BUILDSIMBIOVASCULARMODEL  Build SimBiology model of the vascular patterning network.
%
%   model = buildSimbioVascularModel(params)
%
%   Builds a SimBiology model representing the Carteni et al. (2014)
%   reaction network at a single spatial point (well-mixed compartment).
%   This captures the core biological interactions — morphogen gradient,
%   procambium differentiation, and phloem/xylem toggle-switch — without
%   the spatial diffusion terms needed for full pattern formation.
%
%   For the full spatial simulation, use runVascularModel() instead.

    if nargin < 1
        % Default parameters matching the working spatial model
        params.sigma_S0 = 0.012; params.mu_S0 = 0.015;
        params.S0_star  = 0.35;
        params.sigma_S  = 0.1;   params.k_S    = 5;    params.rho_S  = 0.3;
        params.sigma_AP = 0.001; params.rho_AP = 0.5;  params.mu_AP = 0.015;
        params.AP_star  = 0.8;
        params.sigma_AF = 0.0005; params.rho_AF = 0.5;  params.mu_AF = 0.015;
        params.sigma_AX = 0.0005; params.rho_AX = 0.5;  params.mu_AX = 0.015;
        params.k_R      = 1;      params.rho_R  = 0.02; params.mu_R  = 0.02;
        params.K_AF     = 1.5;    params.K_AX   = 1.5;
        params.isCentral = true;  % true = central domain (S1 active), false = peripheral
    end

    % Add helper scripts
    skillDir = fullfile(getenv('USERPROFILE'), '.matlab', 'agentic-toolkits', ...
                        'matlab', 'skills-catalog', 'computational-biology', ...
                        'matlab-build-simbiology-model', 'scripts');
    if exist(skillDir, 'dir'), addpath(skillDir); end

    % ---- Create model ----
    fprintf('Building SimBiology Vascular Pattern Model...\n');
    model = sbiomodel('VascularPattern');
    fprintf('Model UUID: %s\n', model.uuid);

    % ---- Compartment ----
    comp = addcompartment(model, 'Tissue', 1);
    comp.Constant = true;
    comp.Note = 'Represents a single spatial point in the stem/root cross-section';

    % ---- Species (7 total) ----
    % Initial values set to approximate pre-pattern state
    sp_S0 = addspecies(comp, 'S0', 0.5);       % morphogen (mid-gradient)
    sp_S1 = addspecies(comp, 'S1', 1.0);       % central substrate
    sp_S2 = addspecies(comp, 'S2', 0.1);       % peripheral substrate
    sp_AP = addspecies(comp, 'AP', 0.01);      % procambium activator
    sp_AF = addspecies(comp, 'AF', 0.0);       % phloem activator
    sp_AX = addspecies(comp, 'AX', 0.0);       % xylem activator
    sp_R  = addspecies(comp, 'R',  0.0);       % common repressor

    % All species use arbitrary units (dimensionless model)

    % ---- Parameters ----
    p_sigma_S0 = addparameter(model, 'sigma_S0', params.sigma_S0);
    p_mu_S0    = addparameter(model, 'mu_S0', params.mu_S0);
    p_sigma_S  = addparameter(model, 'sigma_S', params.sigma_S);
    p_k_S      = addparameter(model, 'k_S', params.k_S);
    p_rho_S    = addparameter(model, 'rho_S', params.rho_S);
    p_sigma_AP = addparameter(model, 'sigma_AP', params.sigma_AP);
    p_rho_AP   = addparameter(model, 'rho_AP', params.rho_AP);
    p_mu_AP    = addparameter(model, 'mu_AP', params.mu_AP);
    p_AP_star  = addparameter(model, 'AP_star', params.AP_star);

    p_sigma_AF = addparameter(model, 'sigma_AF', params.sigma_AF);
    p_rho_AF   = addparameter(model, 'rho_AF', params.rho_AF);
    p_mu_AF    = addparameter(model, 'mu_AF', params.mu_AF);
    p_K_AF     = addparameter(model, 'K_AF', params.K_AF);

    p_sigma_AX = addparameter(model, 'sigma_AX', params.sigma_AX);
    p_rho_AX   = addparameter(model, 'rho_AX', params.rho_AX);
    p_mu_AX    = addparameter(model, 'mu_AX', params.mu_AX);
    p_K_AX     = addparameter(model, 'K_AX', params.K_AX);

    p_k_R      = addparameter(model, 'k_R', params.k_R);
    p_rho_R    = addparameter(model, 'rho_R', params.rho_R);
    p_mu_R     = addparameter(model, 'mu_R', params.mu_R);

    % ---- Module 1: S0 morphogen production/decay ----
    rxn = addreaction(model, 'null -> Tissue.S0');
    rxn.ReactionRate = 'sigma_S0';
    rxn.Name = 'S0_production';

    rxn = addreaction(model, 'Tissue.S0 -> null');
    rxn.ReactionRate = 'mu_S0 * Tissue.S0';
    rxn.Name = 'S0_decay';

    % ---- Module 2: Substrate dynamics ----
    % S1 production (central domain — conditionally active)
    if params.isCentral
        rxn = addreaction(model, 'null -> Tissue.S1');
        rxn.ReactionRate = 'sigma_S * (1 - Tissue.S1 / (1 + k_S * Tissue.AP))';
        rxn.Name = 'S1_production';
    end

    % S2 production (peripheral domain)
    if ~params.isCentral
        rxn = addreaction(model, 'null -> Tissue.S2');
        rxn.ReactionRate = 'sigma_S * (1 - Tissue.S2 / (1 + k_S * Tissue.AP))';
        rxn.Name = 'S2_production';
    end

    % Cross-reaction consumption
    rxn = addreaction(model, 'Tissue.S1 + Tissue.S2 -> null');
    rxn.ReactionRate = 'rho_S * Tissue.AP^2 * Tissue.S1 * Tissue.S2';
    rxn.Name = 'Substrate_cross_consumption';

    % ---- Module 2: Procambium activator ----
    rxn = addreaction(model, 'null -> Tissue.AP');
    rxn.ReactionRate = 'sigma_AP + rho_AP * Tissue.AP^2 * Tissue.S1 * Tissue.S2';
    rxn.Name = 'AP_production';

    rxn = addreaction(model, 'Tissue.AP -> null');
    rxn.ReactionRate = 'mu_AP * Tissue.AP';
    rxn.Name = 'AP_decay';

    % ---- Module 3: Phloem/Xylem toggle-switch ----
    % Use a smooth Hill-function gate for procambium-dependent activation
    % pcActive ≈ 0 when AP < AP_star, ≈ 1 when AP > AP_star
    p_nHill = addparameter(model, 'nHill', 10);
    p_pcActive = addparameter(model, 'pcActive', 0);
    p_pcActive.Constant = false;
    addrule(model, 'pcActive = Tissue.AP^nHill / (AP_star^nHill + Tissue.AP^nHill)', 'repeatedAssignment');

    % AF: autocatalytic with Hill saturation, repressed by R
    rxn = addreaction(model, 'null -> Tissue.AF');
    rxn.ReactionRate = 'pcActive * (sigma_AF + rho_AF * Tissue.AF^2 / (K_AF^2 + Tissue.AF^2) / (1 + k_R * Tissue.R))';
    rxn.Name = 'AF_production';

    rxn = addreaction(model, 'Tissue.AF -> null');
    rxn.ReactionRate = 'mu_AF * Tissue.AF';
    rxn.Name = 'AF_decay';

    % AX: autocatalytic with Hill saturation, repressed by R
    rxn = addreaction(model, 'null -> Tissue.AX');
    rxn.ReactionRate = 'pcActive * (sigma_AX + rho_AX * Tissue.AX^2 / (K_AX^2 + Tissue.AX^2) / (1 + k_R * Tissue.R))';
    rxn.Name = 'AX_production';

    rxn = addreaction(model, 'Tissue.AX -> null');
    rxn.ReactionRate = 'mu_AX * Tissue.AX';
    rxn.Name = 'AX_decay';

    % Repressor R: produced by AF + AX
    rxn = addreaction(model, 'null -> Tissue.R');
    rxn.ReactionRate = 'pcActive * rho_R * (Tissue.AF + Tissue.AX)';
    rxn.Name = 'R_production';

    rxn = addreaction(model, 'Tissue.R -> null');
    rxn.ReactionRate = 'mu_R * Tissue.R';
    rxn.Name = 'R_decay';

    % ---- Observables ----
    addobservable(model, 'PhloemDominance', 'Tissue.AF ./ (Tissue.AF + Tissue.AX + 1e-10)');
    addobservable(model, 'TotalVascular', 'Tissue.AF + Tissue.AX');

    % ---- Simulation config ----
    cs = getconfigset(model, 'active');
    cs.StopTime = 100;
    cs.SolverType = 'ode15s';
    cs.RuntimeOptions.StatesToLog = 'all';

    % ---- Verify ----
    fprintf('Verifying model...\n');
    verify(model);
    fprintf('Model built successfully.\n');
    fprintf('  Species: %d\n', numel(model.Species));
    fprintf('  Reactions: %d\n', numel(model.Reactions));
    fprintf('  Parameters: %d\n', numel(model.Parameters));
    fprintf('  Rules: %d\n', numel(model.Rules));
    fprintf('  To simulate: [t, x, names] = sbiosimulate(model);\n');
end
