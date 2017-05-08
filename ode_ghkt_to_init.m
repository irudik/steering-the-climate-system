function [t,y,te,ye] = ode_ghkt_to_init(input)

%% Simulate (in reverse) the state and costate transitions for the
% GHKT models before the negative emissions constraint binds

global params timespace it_count logic

% Omit constraint shadow cost variable
input(7) = []; 


% Create mesh of times to discretize the problem and evaluate the ODE system
timespace = params.tau:-params.grid_res:params.initial_time;

% Set tolerances
options_ode = odeset('AbsTol', 1e-12, 'RelTol', 1e-12, 'Events', @event_M1);%, 'NonNegative', [6]);

% Simulate the system of equations backwards given the terminal conditions
if logic.hotelling
    input = input(1:4);
    [t,y,te,ye] = ode45(@state_path_ghkt_hotel, timespace, input, options_ode);
else
    [t,y,te,ye] = ode45(@state_path_ghkt, timespace, input, options_ode);
end
