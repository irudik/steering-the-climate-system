function [value,isterminal,direction] = event_mu0(t,y)

%% Stop ode solver for event
% Event when negative emissions constraint no longer binds or hit initial
% CO2
global params logic

if logic.stat_ems
    
    value = y(7); % stop when reach mu = 0
    isterminal = 1; % terminate at a zero of value
    direction = 0; % terminate only at a zero at which is decreasing
    
    % Also tell it to stop if finds time at which get to initial value
    value(2) = y(4) - params.initial_co2_perm; % stop when reach initial M1
    isterminal(2) = 1; % terminate at a zero of value
    direction(2) = 0;
    
else
    
    value = y(5); % stop when reach mu=0
    isterminal = 1; % terminate at a zero of value
    direction = 0; % terminate only at a zero at which is decreasing
    
    % Also tell it to stop if finds time at which get to initial value
    value(2) = y(4) - params.initial_co2; % stop when reach initial M1
    isterminal(2) = 1; % terminate at a zero of value
    direction(2) = 0;
    
end


end