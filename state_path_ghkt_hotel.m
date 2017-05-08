function Dlambda = state_path_ghkt_hotel(t,lambda)

global params timespace it_count logic

% Initialize abatement
abatement = ((params.initial_emissions^params.a_2/params.sigma_0/params.a_0/params.initial_gdp).*...
    (params.psi_l*lambda(1)+lambda(3)*params.psi_0*(1-params.psi_l)))^(1/(params.a_2-1));

%% GHKT Hotelling model transition equations when negative emissions
% constraint does not bind

% Equations of motion are represented as: row contents = d/dt(shadow_cost
% or state_var)
Dlambda = [(params.r)*lambda(1);                                           % CO2 permanent co-state, lambda_M1
    params.psi_l*(params.initial_emissions-abatement);                     % CO2 permanent transition, M1 
    (params.r+params.psi)*lambda(3);                                       % CO2 geometric co-state, lambda_M2
    params.psi_0*(1-params.psi_l)*(params.initial_emissions-abatement)-...
    (params.psi)*lambda(4)];                                               % CO2 geometric transition, M2
