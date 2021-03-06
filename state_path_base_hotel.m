function Dlambda = state_path_base_hotel(t,lambda)

global params timespace it_count logic

% Initialize emissions and abatement
emissions = params.initial_emissions;
abatement = ( lambda(1).*emissions.^params.a_2 /...
    (params.a_0*params.initial_gdp*params.sigma_0) ).^(1/(params.a_2-1)); 

%% Base Hotelling model transition equations

% Equations of motion are represented as: row contents = d/dt(shadow_cost
% or state_var)
Dlambda = [(params.r+params.delta)*lambda(1)                          % CO2 co-state, lambda_M 
    emissions-abatement-params.delta*(lambda(2)-params.mpre)];        % CO2 transition, M

