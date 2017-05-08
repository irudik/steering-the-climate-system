function Dlambda = state_path_constr(t,lambda)

global params timespace it_count logic

% Initialize abatement
abatement = params.initial_emissions;

%% GHKT inertia model transition equations when negative emissions
% constraint does not bind

% Equations of motion are represented as: row contents = d/dt(shadow_cost
% or state_var)
Dlambda = [(params.r+params.phi)*lambda(1);                                % Temperature co-state, lambda_T
    (params.r)*lambda(2)-...
    params.phi*params.s*params.alpha*lambda(1)./(lambda(4)+lambda(6));     % CO2 permanent co-state, lambda_M1
    params.phi*(params.s*params.alpha*...
    log((lambda(4)+lambda(6))/params.mpre)-lambda(3));                     % Temperature transition, T 
    params.psi_l*(params.initial_emissions-abatement)                      % CO2 permanent transition, M1 
    (params.r+params.psi)*lambda(5)-...
    params.phi*params.s*params.alpha*lambda(1)./(lambda(4)+lambda(6));     % CO2 geometric co-state, lambda_M2
    params.psi_0*(1-params.psi_l)*(params.initial_emissions-abatement)-...
    (params.psi)*lambda(6)];                                               % CO2 geometric transition, M2

Dlambda(7) = params.psi_l*Dlambda(2) + params.psi_0*(1-params.psi_l)*Dlambda(5);  % shadow cost of constraint
