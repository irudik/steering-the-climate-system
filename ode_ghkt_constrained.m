function [t,y,te,ye,ie] = ode_ghkt_constrained(input)

%% Simulate (in reverse) the state and costate transitions for the
% GHKT model while the negative emissions constraint binds

global params timespace it_count logic

% Seventh element is negative emission shadow cost
input(7) = 0;  

% Create mesh of times to discretize the problem and evaluate the ODE system
timespace = params.tau:-params.grid_res:params.initial_time;

% Set tolerances
options_ode = odeset('AbsTol', 1e-12, 'RelTol', 1e-12, 'Events', @event_mu0);

% Simulate the system of equations backwards given the terminal conditions
if logic.hotelling
    [t,y,te,ye,ie] = ode45(@state_path_ghkt_constr_hotel, timespace, input, options_ode);
else
    [t,y,te,ye,ie] = ode45(@state_path_ghkt_constr, timespace, input, options_ode);
end
