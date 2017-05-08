function [value,isterminal,direction] = event_T0(t,y)

%% Stop ode solver for event
% Event to hit initial temperature

global initial_temp

% Event for temperature equal to initial temperature
value = y(3) - initial_temp; % stop when reach initial temp

isterminal = 1; % terminate at a zero of value
direction = 0; 


end