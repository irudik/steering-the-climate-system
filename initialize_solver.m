%% Set initial guesses of free terminal conditions for each model type


% GHKT decay model
if logic.ghkt
        
    disp(['Loading initial guesses and beginning run with GHKT carbon decay structure, and a temperature target of ' num2str(params.terminal_temp, 3) '.']);
   
    % GHKT hotelling model
    if logic.hotelling 
        
        % Solver bounds for M1(tau)
        l_bounds = [params.initial_co2_perm];
        u_bounds = [params.terminal_co2];
        
        % Search for for M1(tau)
        terminal_search = [748.816];
    
    % GHKT inertia model
    else  
        
        % Solver bounds for lambda_T(tau) and M1(tau)
        l_bounds = [0 params.initial_co2_perm];
        u_bounds = [inf params.terminal_co2];
        
        % Search for for lambda_T(tau), M1(tau)
        if params.terminal_temp == 2
            terminal_search = [2.4561e+07 787.939];
        elseif params.terminal_temp == 2.5
            terminal_search = [2.15621e+07 853.222];
        elseif params.terminal_temp == 3
            terminal_search = [24552089.5969459 0.5*(params.terminal_co2+params.initial_co2_perm)];
        end
        
    end
    
    % Initialize anonymous function where we will search over the above
    % terminal conditions to try to find trajectories that satisfy our
    % initial conditions
    nonlincons = @(x) parameter_search_ghkt(x);
    
% Non-stationary emissions
elseif ~logic.stat_ems
    
    disp(['Loading initial guesses and beginning run with non-stationary emissions, and a temperature target of ' num2str(params.terminal_temp, 3) '.']);
    
    params.terminal_abatement = params.initial_emissions - ...
        params.delta*(params.terminal_co2-params.mpre);
    
    params.shadow_co2 = params.initial_gdp*params.a_0*params.sigma_0*...
        (params.terminal_abatement/params.initial_emissions)^(params.a_2-1)/params.initial_emissions;
    
    
    
    % Non-stationary Hotelling model: search for lambda_T(tau)
    if logic.hotelling
        if params.terminal_temp == 2
            terminal_search = [38.3236786187447];
        elseif params.terminal_temp == 2.5
            terminal_search = [58.0542];
        elseif params.terminal_temp == 3
            terminal_search = [78.9327];
        end
        l_bounds = 0;
        u_bounds = 500;
    % Non-stationary inertia model: search for lambda_T(tau) and tau   
    else
        if params.terminal_temp == 2
            terminal_search = [2747.59623483090	130.785659574311];
        elseif params.terminal_temp == 2.5
            terminal_search = [2476.67912419411	158.978898290394];
        elseif params.terminal_temp == 3
            terminal_search = [2245.62021062529	185.833980250711];
        end
        l_bounds = [0 -200];
        u_bounds = [Inf 2000];
    end
    
    % Initialize anonymous function where we will search over the above
    % terminal conditions to try to find trajectories that satisfy our
    % initial conditions
    nonlincons = @(x) parameter_search_nonstat(x);
    
% Base model
else
    
    disp(['Loading initial guesses and beginning base specification run, and a temperature target of ' num2str(params.terminal_temp, 3) '.']);
    
    % Base inertia model, base Hotelling does not need a parameter search:
    % search for lambda_T(tau)
    if ~logic.hotelling
        params.tau = 700;
        run load_initial_base
        l_bounds = 0;
        u_bounds = Inf;
    end
    
    % Initialize anonymous function where we will search over the above
    % terminal conditions to try to find trajectories that satisfy our
    % initial conditions
    nonlincons = @(x) parameter_search_base(x);
    
    
end

% Reset iteration counter for reporting results
it_count = 1;