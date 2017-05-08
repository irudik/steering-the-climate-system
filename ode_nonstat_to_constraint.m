function [t,y,te,ye,ie] = ode_nonstat_to_constraint(terminal_search)

%% Simulate (in reverse) the state and costate transitions for the
% non-stationary Hotelling model, and the non-stationary inertia model
% after the negative emissions constraint stops binding


global params timespace it_count logic

% Create mesh of times to discretize the problem and evaluate the ODE system
if logic.hotelling
    timespace = terminal_search(1):-params.grid_res:-200;
    nonstat_time = 0;
    options_ode = odeset('RelTol',  1e-12, 'AbsTol', 1e-12,'Events', @event_time0);

else
    timespace = terminal_search(2):-params.grid_res:-200;
    shadow_temp = terminal_search(1);
    nonstat_time = 0;
    options_ode = odeset('RelTol',  1e-12, 'AbsTol', 1e-12, 'Events', @event_fullabate);

end

% Initialize non-stationary processes given the 'start time' guess
emissions = params.initial_emissions*exp(.0068*(timespace));
terminal_abatement_nonstat = emissions(1)-params.delta*(params.terminal_co2-params.mpre);
shadow_co2_nonstat = params.initial_gdp*params.a_0*params.sigma_0*(terminal_abatement_nonstat).^(params.a_2-1)/emissions(1)^params.a_2;

% Simulate the system of equations backwards given the terminal conditions
if logic.hotelling
    [t,z,te,ye,ie] = ode23(@(t,z) state_path_nonstat_hotel(t, z, timespace, emissions), ...
        timespace, [shadow_co2_nonstat params.terminal_co2], options_ode);
    
    % Adjust output matrix to account for the fact that Hotelling runs do
    % not have any temperature states
    y(:,2) = z(:,1);
    y(:,4) = z(:,2);
    
    % Emissions
    y(:,5) =  params.initial_emissions*exp(.0068*t);
    
    % Abatement
    y(:,6)=((y(:,5).^params.a_2).*y(:,2)./(params.a_0*params.initial_gdp*params.sigma_0)).^(1/(params.a_2-1)); %abatement

    % Time
    y(:,7)=t;
    
else
    [t,y,te,ye,ie] = ode23(@(t,y) state_path_nonstat(t, y, timespace, emissions), ...
        timespace, [shadow_temp(1) shadow_co2_nonstat params.terminal_temp params.terminal_co2], options_ode);
    
    % Emissions
    y(:,5) =  params.initial_emissions*exp(.0068*(t));
    
    % Abatement
    y(:,6)=((y(:,5).^params.a_2).*y(:,2)./(params.a_0*params.initial_gdp*params.sigma_0)).^(1/(params.a_2-1)); %abatement
    
    % Time
    y(:,7)=t;
    
end
