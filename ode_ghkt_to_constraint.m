function [t,y,te,ye,ie] = ode_ghkt_to_constraint(terminal_search)

%% Simulate (in reverse) the state and costate transitions for the
% GHKT model after the negative emissions constraint stops binding

global params timespace logic

% Initialize terminal conditions for Hotelling and non-Hotelling runs
if logic.hotelling
    
    end_co2_perm = terminal_search(1);
    
    % Post-tau time vector, go 800 years to approx integrating to infinity
    inf_time = params.tau:params.grid_res:params.tau+800;
    
    % Post-tau abatement vector
    inf_abate = params.initial_emissions - ...
        exp(-(params.psi*params.psi_l/(params.psi_l+(params.psi_0*(1-params.psi_l))))*(inf_time-params.tau)).*...
        (params.psi/(params.psi_l+params.psi_0*(1-params.psi_l))*(params.terminal_co2-end_co2_perm));
    
    % Exponential part of the integral
    int_1 = exp(-(params.r + params.psi*params.psi_l./(params.psi_l+params.psi_0*(1-params.psi_l))).*(inf_time-params.tau));
    
    % MAC part of the integral
    int_2 = (params.a_0*params.sigma_0*params.initial_gdp./(params.initial_emissions.^params.a_2))*(inf_abate.^(params.a_2-1));
    
    % Full integral from lambda_M1(tau+) condition
    integral = params.psi/(params.psi_l+params.psi_0*(1-params.psi_l))*trapz(inf_time, int_1.*int_2);

    
    % Time tau abatement and MAC
    end_abate = params.initial_emissions - params.psi*(params.terminal_co2-end_co2_perm)/(params.psi_l+params.psi_0*(1-params.psi_l));
    end_mac = (params.a_0*params.sigma_0*params.initial_gdp/(params.initial_emissions^params.a_2))*(end_abate^(params.a_2-1));
    
    % Constraint shadow price
    nu = (params.psi_l*integral - end_mac) / ((params.psi_l+params.psi_0*(1-params.psi_l)));
    disp(['Jump coefficient: ' num2str(nu,4)])
    
    shadow_co2_perm = integral - nu;
    shadow_co2_geo = - nu;
    
    end_co2_geo = params.terminal_co2 - end_co2_perm;
    

else
    
    shadow_temp = terminal_search(1);
    end_co2_perm = terminal_search(2);
    
    % Calculate integral term from appendix to get CO2 shadow prices
    
    % Post-tau time vector, go 800 years to approx integrating to infinity
    inf_time = params.tau:params.grid_res:params.tau+800;
    
    % Post-tau abatement vector
    inf_abate = params.initial_emissions - ...
        exp(-(params.psi*params.psi_l/(params.psi_l+(params.psi_0*(1-params.psi_l))))*(inf_time-params.tau)).*...
        (params.psi/(params.psi_l+params.psi_0*(1-params.psi_l))*(params.terminal_co2-end_co2_perm));
    
    % Exponential part of the integral
    int_1 = exp(-(params.r + params.psi*params.psi_l./(params.psi_l+params.psi_0*(1-params.psi_l))).*(inf_time-params.tau));
    
    % MAC part of the integral
    int_2 = (params.a_0*params.sigma_0*params.initial_gdp./(params.initial_emissions.^params.a_2))*(inf_abate.^(params.a_2-1));
    
    % Full integral from lambda_M1(tau+) condition
    integral = params.psi/(params.psi_l+params.psi_0*(1-params.psi_l))*trapz(inf_time, int_1.*int_2);

    end_co2_geo = params.terminal_co2 - end_co2_perm;
    
    % Time tau abatement and MAC
    end_abate = params.initial_emissions - params.psi*(params.terminal_co2-end_co2_perm)/(params.psi_l+params.psi_0*(1-params.psi_l));
    end_mac = (params.a_0*params.sigma_0*params.initial_gdp/(params.initial_emissions^params.a_2))*(end_abate^(params.a_2-1));
    
    % Constraint shadow price
    nu = (params.psi_l*integral - end_mac) / (params.s*params.phi*params.alpha/params.terminal_co2/(params.psi_l+params.psi_0*(1-params.psi_l)));
    disp(['Jump coefficient: ' num2str(nu,4)])
    
    shadow_co2_perm = integral - nu*(params.s*params.phi*params.alpha/params.terminal_co2);
    shadow_co2_geo = - nu*(params.s*params.phi*params.alpha/params.terminal_co2);
    
    end_temp = params.terminal_temp;
    
end

% Create mesh of times to discretize the problem and evaluate the ODE system
timespace = params.tau:-params.grid_res:params.initial_time;

% Set tolerances
options_ode = odeset('AbsTol', 1e-12, 'RelTol', 1e-12, 'Events', @event_fullabate);

% Simulate the system of equations backwards given the terminal conditions
if logic.hotelling
    terminal_conditions = [shadow_co2_perm; end_co2_perm; shadow_co2_geo; end_co2_geo];
    [t,y,te,ye,ie] = ode45(@state_path_ghkt_hotel, timespace, terminal_conditions, options_ode);
else
    terminal_conditions = [shadow_temp; shadow_co2_perm; end_temp;...
    end_co2_perm; shadow_co2_geo; end_co2_geo];
    [t,y,te,ye,ie] = ode45(@state_path_ghkt, timespace, terminal_conditions, options_ode);
end
