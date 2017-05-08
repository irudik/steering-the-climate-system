function [value,isterminal,direction] = event_M1(t,y)

%% Stop ode solver for event
% Event for hitting initial CO2 stock

global params logic


value = y(2) - params.initial_co2; % stop when reach initial CO2 stock
isterminal = 1; % terminate at a zero of value
direction = 0;


end