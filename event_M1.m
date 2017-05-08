function [value,isterminal,direction] = event_M1(t,y)

%% Stop ode solver for event
% Event for hitting initial CO2 stock or initial time

global params logic

if ~logic.stat_ems
    % Event for hitting initial CO2 stock
    value = y(4) - params.initial_co2; % stop when reach initial M
    isterminal = 1; % terminate at a zero of value
    direction = 0;
    
    % Event for hitting initial time
    value(2) = t; % stop when reach initial time!
    isterminal(2) = 1; % terminate at a zero of value
    direction(2) = 0;
    
else
    
    % Event for hitting initial permanent CO2 stock
    value = y(4) - params.initial_co2_perm; % stop when reach initial M1
    
    isterminal = 1; % terminate at a zero of value
    direction = 0;
end

end