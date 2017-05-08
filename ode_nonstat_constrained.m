function [t,y,te,ye,ie] = ode_nonstat_constrained(input)

%% Simulate (in reverse) the state and costate transitions for the
% non-stationary inertia model while the negative emissions constraint
% binds

global params timespace it_count logic

% Create mesh of times to discretize the problem and evaluate the ODE system
timespace = fliplr(-200:params.grid_res:input(5));

% Initialize non-stationary processes
input(5) = 0;
emissions = params.initial_emissions*exp(.0068*(timespace));

% Set tolerances for the ODE solver
options_ode = odeset('AbsTol', 1e-12, 'RelTol', 1e-12, 'Events', @event_mu0);

% Simulate the system of equations backwards given the terminal conditions
[t,y,te,ye,ie] = ode45(@(t,y) state_path_nonstat_constr(t, y, timespace, emissions), timespace, input, options_ode);

% Emissions
y(:,5) =  params.initial_emissions*exp(.0068*(t));

% Abatement
y(:,6)=((y(:,5).^params.a_2).*y(:,2)./(params.a_0*params.initial_gdp*params.sigma_0)).^(1/(params.a_2-1)); %abatement

% Time
y(:,7)=t;

        