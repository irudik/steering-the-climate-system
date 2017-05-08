function [t,y,te,ye,ie] = ode_base_to_constraint(terminal_search)

%% Simulate (in reverse) the state and costate transitions for the base model

global params timespace it_count logic

% Create mesh of times to discretize the problem and evaluate the ODE system
timespace = params.tau:-params.grid_res:params.initial_time;

% Calculate terminal CO2 shadow price
terminal_abatement = params.initial_emissions-params.delta*(params.terminal_co2-params.mpre);
shadow_co2 = params.initial_gdp*params.a_0*params.sigma_0*(terminal_abatement/params.initial_emissions)^(params.a_2-1)/params.initial_emissions;


% Simulate the system of equations backwards given the terminal conditions
if logic.hotelling
    
    % Set tolerances for the ODE solver
    options_ode = odeset('RelTol',  1e-12, 'AbsTol', 1e-12, 'Events', @event_M);
    
    [t,z] = ode23(@state_path_base_hotel, timespace, [shadow_co2 params.terminal_co2], options_ode);
    
    % Adjust output matrix to account for the fact that Hotelling runs do
    % not have any temperature states
    y(:,2) = z(:,1);
    y(:,4) = z(:,2);
    
    % Emissions
    y(:,5) = params.initial_emissions;
    
    % Abatement
    y(:,6) = ((y(:,5).^params.a_2).*y(:,2)/(params.a_0*params.initial_gdp*params.sigma_0)).^(1/(params.a_2-1));
    
    % Time
    y(:,7) = t;
    
else
    
    % Set tolerances for the ODE solver
    options_ode = odeset('RelTol',  1e-10, 'AbsTol', 1e-10, 'Events', @event_fullabate);
    
    [t,y,te,ye,ie] = ode23(@(t,y) state_path_base(t, y), ...
        timespace, [terminal_search shadow_co2 params.terminal_temp params.terminal_co2], options_ode);
    
    % Emissions
    y(:,5) = params.initial_emissions;
    
    % Abatement
    y(:,6) = ((y(:,5).^params.a_2).*y(:,2)/(params.a_0*params.initial_gdp*params.sigma_0)).^(1/(params.a_2-1));
    
    % Time
    y(:,7) = t;
end

