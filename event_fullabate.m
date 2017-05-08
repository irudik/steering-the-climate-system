function [value,isterminal,direction] = event_fullabate(t,y)

%% Stop ode solver for event
% Event to hit full abatement or initial CO2

global params logic

% Marginal abatement cost at full abatement
mac_full = params.a_0*params.sigma_0*params.initial_gdp/params.initial_emissions;

if logic.ghkt
    
    % Event for full abatement
    if logic.hotelling
        value(1) = mac_full - params.psi_l*y(1) - params.psi_0*(1-params.psi_l)*y(3); % stop when reach A=E
    else
        value(1) = mac_full - params.psi_l*y(2) - params.psi_0*(1-params.psi_l)*y(5); % stop when reach A=E
    end
    
    isterminal(1) = 1; % terminate at a zero of value
    direction(1) = 0;
    
    % Event for hitting the initial level of the permanent CO2 stock
    if logic.hotelling
        value(2) = y(2) - params.initial_co2_perm; % stop when reach initial M!
    else
        value(2) = y(4) - params.initial_co2_perm; % stop when reach initial M!
    end
    
    isterminal(2) = 1; % terminate at a zero of value
    direction(2) = 0;
    
elseif ~logic.stat_ems
    
    % Initialize non-stationary processes given the 'start time' guess
    emissions = params.initial_emissions*exp(.0068*(t));
    
    % Full abatement
    mac_full = params.a_0*params.sigma_0*params.initial_gdp./emissions;
    
    % Event for full abatement
    value(1) = mac_full - y(2); % stop when reach A=E
    
    isterminal(1) = 1; % terminate at a zero of value
    direction(1) = 0;
    
    % Event for hitting the initial level of the permanent CO2 stock
    if logic.hotelling
        value(2) = y(2) - params.initial_co2;
    else
        value(2) = y(4) - params.initial_co2; % stop when reach initial M!
    end
    
    isterminal(2) = 1; % terminate at a zero of value
    direction(2) = 0;
    
else
    
    % Initialize non-stationary processes given the 'start time' guess
    emissions = params.initial_emissions;
    
    % Full abatement
    mac_full = params.a_0*params.sigma_0*params.initial_gdp./params.initial_emissions;
    
    % Event for full abatement
    value(1) = mac_full - y(2); % stop when reach A=E
    
    isterminal(1) = 1; % terminate at a zero of value
    direction(1) = 0;
    
    % Event for hitting the initial level of the permanent CO2 stock
    if logic.hotelling
        value(2) = y(2) - params.initial_co2;
    else
        value(2) = y(4) - params.initial_co2; % stop when reach initial M!
    end
    
    isterminal(2) = 1; % terminate at a zero of value
    direction(2) = 0;    
    
end

end