% Replication code for:
% Steering the Climate System: Using Inertia to Lower the Cost of Policy
% by Derek Lemoine and Ivan Rudik
% American Economic Review
% The git repository can be found at:
% http://github.com/irudik/steering-the-climate-system

clear all
clear global
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set-up model and user options

% Set global parameters
global params it_count logic

% Set solver options file name for KNITRO solver
options_file = 'knitro_options_ijr.opt';

% Switch for type of solver, note: code has only been tested with KNITRO
% and current settings in the options file
logic.knitro = 1; % =1 for KNITRO; =0 for fmincon.


%% Switches for different types of runs
% The GHKT carbon decay structure and non-stationary emissions only run for
% the base discount rate and inertia settings
logic.stat_ems = 1; % =1 for stationary emissions
logic.ghkt = 0; % =1 to use GHKT decay structure
logic.hotelling = 0; % =1 to solve the Hotelling (CO2 target) problem
logic.vary_params = 0; % =1 to loop and vary discount rate and inertia
logic.solve_nonstat = 0; % =1 to solve the non-stationary problem, will take a VERY long time; =0 to simulate based on the pre-solved initial guesses

%% Initialize model parameters
run initialize_parameters

% Set temperature targets to loop over
params.first_temp = 2; % Set lowest terminal temperature to loop through
params.final_temp = 2; % Set highest terminal temperature to loop through
params.temp_incr = 0.5; % Temperature increment for looping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve model and simulate trajectories
% Loop through phi
for phis=1:numel(fields_phi)
    
    % Loop through r
    for rs=1:numel(fields_r)
        
        % Only do runs with one parameter changed from base specification
        if ( R.(fields_r{rs}) ~= R.base && Phi.(fields_phi{phis}) ~= Phi.base )
            break;
        end
        
        % Loop through temperature targets
        for temp=params.first_temp:params.temp_incr:params.final_temp
            % Clear output storage matrices between runs
            clear inertia_out hotelling_out terminal_conditions
            
            % Set inertia and discounting for this run
            params.phi = Phi.(fields_phi{phis});
            params.r = R.(fields_r{rs});
            
            % Initial (actually, terminal) conditions
            params.terminal_temp = temp;
            params.terminal_co2 = params.mpre*exp(params.terminal_temp/(params.s*params.alpha));
            
            % Initialize solver guesses and bounds
            run initialize_solver
            
            % Run solver
            disp('Solving model.')
            
            % Base Hotelling runs do not need to search over terminal
            % conditions, so simply simulate the model in reverse until we
            % hit the initial conditions
            if logic.hotelling && logic.stat_ems && ~logic.ghkt
                
                % We can infer the terminal abatement level and terminal
                % CO2 shadow cost from the CO2 target
                terminal_abatement = params.initial_emissions - ...
                    params.delta*(params.terminal_co2 - params.mpre);
                terminal_co2_shadow = params.initial_gdp*params.a_0*params.sigma_0*...
                    (terminal_abatement/params.initial_emissions)^(params.a_2-1)/params.initial_emissions;
                
                % With this information, simulate the model in reverse and
                % recover the simulated state trajectories in
                % 'inertia_out'
                [t,inertia_out] = ode_base_to_constraint(terminal_co2_shadow);

                % Adjust the time vector so time starts at 0
                inertia_out(:,7) = inertia_out(:,7)-inertia_out(end,7);
                
                % Flip to zero to plot below
                exitflag = 0;
                
                % If we aren't solving the base Hotelling model, we need to
                % search over terminal conditions
            else
                
                % Search over terminal conditions for all other model types
                % minimize an objective function of f(x) = 0 subject to a
                % model-specific set of nonlinear equality constraints
                % which require initial conditions to be satified
                if ~logic.stat_ems && ~logic.hotelling && ~logic.solve_nonstat
                    terminal_conditions = terminal_search;
                    exitflag = 0;
                else
                    if logic.knitro
                        [terminal_conditions, ~, exitflag] = ...
                            knitromatlab(@(x) 0,terminal_search,[],[],[],[],l_bounds,u_bounds,nonlincons,[],[],options_file);
                    else
                        [terminal_conditions, ~, exitflag] = ...
                            fmincon(@(x) 0,terminal_search,[],[],[],[],l_bounds,u_bounds,nonlincons);
                    end
                end
            end
            
            % Simulate trajectories forward in time if solver returns a
            % successful exitflag.
            if exitflag > -203
                
                % If using GHKT decay structure
                if logic.ghkt
                    
                    % Simulate the trajectories forward
                    inertia_out = simulate_trajectories_ghkt(terminal_conditions);
                    inertia_out = flipud(inertia_out);
                    
                    % Set initial temperature and infer temperature
                    % trajectory for the Hotelling GHKT model
                    inertia_out(1,3) = params.initial_temp;
                    if logic.hotelling
                        for i = 2:size(inertia_out,1)
                            inertia_out(i,3) = inertia_out(i-1,3) + ...
                                params.grid_res*params.phi*...
                                (params.s*params.alpha*log((inertia_out(i-1,4)+inertia_out(i-1,6))/params.mpre)-inertia_out(i-1,3));
                        end
                    end
                    
                    % Generate time mesh for the simulation
                    inertia_out(:,8) = params.initial_time:params.grid_res:(size(inertia_out,1)-1)*params.grid_res;
                    inertia_out = flipud(inertia_out);
                    
                    % Simulate post-tau and calculate results
                    run calculate_results_ghkt
                    
                    % Else if using non-stationary emissions
                elseif ~logic.stat_ems
                    
                    % Simulate post-tau and calculate results
                    [inertia_out] = simulate_trajectories_nonstat(terminal_conditions);
                    run calculate_results_non_ghkt
                    
                    % Else using base model
                else
                    
                    % If non-hotelling, simulate the model
                    if ~logic.hotelling
                        [inertia_out] = simulate_trajectories_base(terminal_conditions);
                        inertia_out(:,7) = inertia_out(:,7)-inertia_out(end,7);
                    end
                    
                    % Simulate post-tau and calculate results
                    run calculate_results_non_ghkt
                    
                end
                
                % Throw error if shadow cost search fails
            else
                [inertia_out] = ode_solver_reverse_time(terminal_conditions);
                error('Solution not found.');
            end
            
        end
        
    end
end


