function [value,isterminal,direction] = event_t0(t,y)

%% Stop ode solver for event
% Event when hit initial time of 0

global params

% Event for hitting initial time
value = t - params.initial_time; % stop when t=0

isterminal = 1; % terminate at a zero of value
direction = 0;




end