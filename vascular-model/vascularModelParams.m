function params = vascularModelParams(varargin)
% VASCULARMODELPARAMS  Define parameters for Carteni et al. (2014) vascular model.
%
%   params = vascularModelParams()           returns default parameters
%   params = vascularModelParams('Name',Value,...)  override specific parameters
%
%   Reference: Carteni F, Giannino F, Schweingruber FH, Mazzoleni S (2014)
%   Modelling the development and arrangement of the primary vascular
%   structure in plants. Annals of Botany 114(4): 619-627.

    % =========================================================================
    % Default parameters from the paper
    % =========================================================================

    % ---- Spatial domain ----
    params.gridSize  = 50;          % N x N grid points
    params.radius    = 20;          % domain radius (spatial units, r in paper)
    params.h         = [];          % grid spacing (computed from radius/gridSize)
    params.domainType = 'circular'; % 'circular' or 'square'

    % ---- Module 1: Morphogen S0 gradient (Eq 1) ----
    params.sigma_S0  = 0.012;       % S0 basic production rate
    params.mu_S0     = 0.015;       % S0 consumption rate
    params.D_S0      = 0.8;         % S0 diffusion coefficient
    params.S0_star   = 0.5;         % threshold for central/peripheral zones

    % ---- Module 2: Procambium differentiation (Eqs 2-4) ----
    % Substrates S1 (central), S2 (peripheral), Activator A_P
    params.sigma_S   = 0.04;        % S basic production rate
    params.k_S       = 20;          % S production saturation constant
    params.rho_S     = 0.08;        % S cross-reaction coefficient
    params.D_S       = 0.5;         % S diffusion coefficient

    params.sigma_AP  = 0.001;       % A_P basic production rate
    params.rho_AP    = 0.03;        % A_P cross-reaction coefficient (0.03 or 0.05)
    params.mu_AP     = 0.02;        % A_P removal rate
    params.D_AP      = 0.02;        % A_P diffusion coefficient
    params.AP_star   = 0.5;         % threshold for procambium differentiation

    % ---- Module 3: Phloem/Xylem differentiation (Eqs 5-9) ----
    % Activators A_F (phloem), A_X (xylem), Repressor R
    params.sigma_AF  = 0.01;        % A_F basic production rate
    params.sigma_AX  = 0.01;        % A_X basic production rate
    params.rho_AF    = 0.1;         % A_F autocatalytic strength
    params.rho_AX    = 0.1;         % A_X autocatalytic strength
    params.k_R       = 1.0;         % repressor saturation constant
    params.mu_AF     = 0.01;        % A_F removal rate
    params.mu_AX     = 0.01;        % A_X removal rate
    params.D_AF      = 0.002;       % A_F diffusion coefficient (range 0.001-0.003)
    params.D_AX      = 0.002;       % A_X diffusion coefficient (range 0.001-0.003)

    params.rho_R     = 0.5;         % repressor production rate (mu_SXF in paper)
    params.mu_R      = 0.5;         % repressor removal rate
    params.D_R       = 0.02;        % repressor diffusion coefficient

    params.A_star    = 30;          % threshold for phloem/xylem differentiation

    % ---- Simulation settings ----
    params.tSpan     = [0, 500];    % simulation time span
    params.solver    = 'ode15s';    % stiff solver for reaction-diffusion
    params.relTol    = 1e-6;        % relative tolerance
    params.absTol    = 1e-8;        % absolute tolerance

    % ---- Override with user-specified values ----
    for i = 1:2:numel(varargin)
        if isfield(params, varargin{i})
            params.(varargin{i}) = varargin{i+1};
        else
            warning('vascularModelParams:UnknownParameter', ...
                    'Unknown parameter: %s', varargin{i});
        end
    end

    % ---- Compute derived quantities ----
    params.h = 2 * params.radius / params.gridSize;

end
