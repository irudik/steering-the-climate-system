%% Initialize model parameters

% Climate parameters
params.delta = .0138; % Base decay (unitless)
params.alpha = 5.35; % Forcing parameter (W/m^2)
params.s = 0.809; % Transformed climate sensitivity (deg C / [W/m^2])
params.mpre = 596.4; % Preindustrial CO2 (GtC)

% Golosov decay parameters: must map decadal parameters into our setting
params.psi_l = 0.2; % permanent fraction of emissions (unitless)
params.psi_0 = 0.393; % Fraction of non-permanent emissions that decays geometrically (unitless)
params.psi = 0.0228/10; % Geometric decay rate (unitless)

% Abatement cost parameters
params.a_0 = 1.17; 
params.a_1 = 2;
params.a_2 = 2.8;
params.sigma_0 = 0.13;
params.g_sig = -0.0073;
params.g_psi = -0.005;
params.delta_sig = 0.003;

% Initial conditions
params.initial_co2 = 808.9; % (GtC)
params.initial_temp = 0.7307; % (deg C)
params.initial_emissions = 9.9662; % (GtC)
params.initial_gdp = 84.6902; % (trillions of dollars)
params.initial_co2_perm = 684; % (GtC)
params.initial_co2_geo = 118; % (GtC)

% Computation parameters
params.initial_time = 0; % Default: 0
params.time_horizon = 400; % Terminal time to simulate the model beyond the target
params.grid_res = .005; % Resolution of the returned solution mesh in units of years

% Search over a longer time horizon if using the GHKT decay structure
if ~logic.ghkt 
    params.tau = 200; % (years)
else
    params.tau = 1000; % (years)
end

% If GHKT, do not vary parameters and do not use non-stationary emissions,
% if non-stationary emissions do not vary parameters
if logic.ghkt
    logic.vary_params = 0;
    logic.stat_ems = 1;
    disp('Golosov decay structure can only run with stationary emissions and base parameterization. Switching to stationary emissions and base parameters.');
elseif ~logic.stat_ems
    logic.vary_params = 0;
    disp('Non-stationary emissions can only run with base parameterization. Switching to base parameterization.');
end

% Initialize high, medium and low parameter values for paper results
Phi.base = .0091; % Inertia
R.base = .055; % Discount rate

% If varying inertia and discount, initialize a structure with the
% alternative specifications
if logic.vary_params
    Phi.h = Phi.base*2;
    Phi.l = Phi.base/2;
    R.l = .014;
end

fields_phi = fieldnames(Phi);
fields_r = fieldnames(R);