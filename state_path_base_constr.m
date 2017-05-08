function Dlambda = state_path_base_constr(t,lambda)

global params timespace it_count logic

% Initialize emissions and abatement
emissions = params.initial_emissions;
abatement = ( lambda(2).*emissions.^params.a_2 / (params.a_0*params.initial_gdp*params.sigma_0) ).^(1/(params.a_2-1)); 

%% Base model transition equations

% Equations of motion are represented as: row contents = d/dt(shadow_cost
% or state_var)
Dlambda = [(params.r+params.phi)*lambda(1)                          % Temperature co-state, lambda_T 
    (params.r+params.delta)*lambda(2) ...
        - params.phi*params.s*params.alpha*lambda(1)./lambda(4)     % CO2 co-state, lambda_M 
    params.phi*...
        (params.s*params.alpha*log(lambda(4)/params.mpre)-lambda(3))% Temperature transition, T 
    emissions-abatement-params.delta*(lambda(4)-params.mpre)];      % CO2 transition, M
Dlambda(5) = Dlambda(2);