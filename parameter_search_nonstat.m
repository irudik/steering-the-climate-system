function [cineq, ceq] = parameter_search_nonstat(terminal_search)

%% Search for terminal conditions that solve the model

% Non-stationary model: search for lambda_T(tau) for the Hotelling model,
% and lambda_T(tau) and tau for the inertia model such that Hotelling: CO2
% at t0 is equal to initial CO2 Least-cost: CO2 and temperature at t0 are
% equal to initial values, allow for negative emissions constraint to bind


global params timespace it_count logic

cineq = 0; % Automatically satisfy inequality constraints

% Call the solver for the passed initial time for non-stationary processes
[t1,y1,te1,ye1,ie1] = ode_nonstat_to_constraint(terminal_search);

if logic.hotelling
    
    disp(['Non-stationary hotelling run ' num2str(it_count) '.  T0: ' num2str(terminal_search) '.']);
    y = y1;
    
    % Set equality constraint to be M(t0) = initial CO2
    ceq = y(end,4) - params.initial_co2;
    
    % If we happen to select perfectly between two points pick the first one
    if numel(ceq) == 2
        ceq(2) = [];
    end
    
    % Return error in initial emissions/gdp level or note that we have a
    % divergent CO2 trajectory and force a large error to re-adjust the initial
    % guess (this doesn't always work)
    if ~isempty(ceq)
        disp(['True initial CO2: ' num2str(params.initial_co2,4), ' simulated initial CO2: ' num2str(y(end,4),4)]);
        disp(['True initial emissions: ' num2str(params.initial_emissions,4), ' simulated initial emissions: ' num2str(y(end,5),4)]);
    else
        disp('Monotonically increasing CO2, readjusting feasiblity error.');
        ceq = [30000];
    end
    
    it_count = it_count+1; % Advance counter
    
else 
    %% Find last time that negative emission constraint binds

    disp(['Non-stationary least-cost run ' num2str(it_count) '.  lambda_T, T0: ' num2str(terminal_search(1)) ' ' num2str(terminal_search(2)) '.']);
    
    % If we cannot find the first time, change the time horizon
    if isempty(ye1)
        terminal_search(2) = terminal_search(2)*2;
        [t1,y1,te1,ye1,ie1] = ode_nonstat_to_constraint(terminal_search);
        if isempty(ye1)
            terminal_search(2) = terminal_search(2)/4;
            [t1,y1,te1,ye1,ie1] = ode_nonstat_to_constraint(terminal_search);
        end
        if isempty(ye1)
            disp('Never found time at which either constraint binds or hit initial CO2.');
            ceq = [-30000];
            return;
        end
    end
    
    
    if size(ie1,1)>1 % first instant may have had an event
        ie1 = ie1(end,:);
        ye1 = ye1(end,:);
        te1 = te1(end,:);
    end
    
    
    
    if ie1 == 2 % if the negative emissions constraint never binds
        
        % Hit initial time
        disp(['Found initial value with ' num2str(terminal_search(2)-te1) ' years from there until guessed steady-state-like conditions.']);
        y = y1;
        t = t1;
        ye3 = ye1;
        
    else % ie1 == 1, negative emissions constraint is binding so we must
         % simulate the model where abatement is constrained
        
        disp(['Negative emission constraint binds at state ' ...
            mat2str(ye1,4) ', with ' num2str(terminal_search(2)-te1,4)...
            ' years from there until guessed steady-state-like conditions.']);
        
        
        %% Use the simulated results for after the constraint has stopped binding to find the
        % first time that negative emission constraint starts binding, or we hit initial temperature
        ye1 = [ye1 te1];
        [t2,y2,te2,ye2,ie2] = ode_nonstat_constrained(ye1);
        
        % If we cannot find the first time, change the time horizon
        if isempty(ye2)
            terminal_search(2) =  terminal_search(2)*4;
            [t2,y2,te2,ye2,ie2] = ode_nonstat_constrained(ye1);
            if isempty(ye2)
                terminal_search(2) =  terminal_search(2)/16;
                [t2,y2,te2,ye2,ie2] = ode_nonstat_constrained(ye1);
            end
            if isempty(ye2)
                disp('Never found first time at which constraint started binding.');
                ceq = [-30000];
                return;
            end
        end
        
        if size(ie2,1) > 1 % first instant should have had an event
            ie2 = ie2(end,:);
            ye2 = ye2(end,:);
            te2 = te2(end,:);
        end
        
        % New time
        t2=t2+(t1(end)-terminal_search(2));
        
        if ie2 == 2
            
            % found time with initial temperature, so store output and skip next
            % simulations
            disp(['Found initial value with ' ...
                num2str( terminal_search(2)-te2)...
                ' years from there until constraint stopped binding.']);
            
            y = vertcat(y1,y2);
            t = vertcat(t1,t2);
            ye3 = ye2;
            
        else % ie2 == 1, so found first time with full abatement
            
            disp(['Negative emission constraint first binds at state '...
                mat2str(ye2,4) ' and then binds for '...
                num2str( terminal_search(2)-te2,4) ' years.']);
            
            %% Now simulate backwards to the initial time
            ye2(5) = te2;
            [t3,y3,te3,ye3,ie3] = ode_nonstat_to_init(ye2);
            
            if all(y3(end,:)==y3(end-1,:))
                y3=y3(1:end-1,:);
            end
            
            % If we cannot find the first time, change the time horizon
            if isempty(ye3)
                
                terminal_search(2) = terminal_search(2)*4;
                [t3,y3,te3,ye3] = ode_nonstat_to_init(ye2);
                
                if isempty(ye3)
                    
                    terminal_search(2) = terminal_search(2)/16;
                    [t3,y3,te3,ye3] = ode_nonstat_to_init(ye2);
                    
                end
                
                if isempty(ye3)
                    disp('Never found initial value.');
                    ceq = [-30000];
                    return;
                    
                end
                
            end
            
            disp(['Negative emission constraint first binds after '...
                mat2str(terminal_search(2)-te3) ' years.']);
            
            t3=t3+(t2(end)-terminal_search(2));
            
            y = vertcat(y1(:,1:7),y2(:,1:7),y3(:,1:7));
            t = vertcat(t1,t2,t3);
            
        end
        
    end
    
    % Set equality constraint to be M(t0) = initial CO2, T(t0) = initial
    % temp
    ceq = [y(end,3)-params.initial_temp; (y(end,4)-params.initial_co2); y(end,5)-params.initial_emissions];
    
    disp(['True initial temp: ' num2str(params.initial_temp,4), ' simulated initial temp: ' num2str(y(end,3),4)]);
    disp(['True initial CO2: ' num2str(params.initial_co2,4), ' simulated initial CO2: ' num2str(y(end,4),4)]);
    disp(['True initial emissions: ' num2str(params.initial_emissions,4), ' simulated initial emissions: ' num2str(y(end,5),4)]);
    disp(['Error: ' num2str(ceq(1,1)) ' ' num2str(ceq(2,1)) ' ' num2str(ceq(3,1)) '.']);
    % Return error in initial emissions/gdp level and temperature or note that
    % we have a divergent CO2 trajectory and force a large error to re-adjust
    % the initial guess (this doesn't always work)
    if isempty(ceq)
        
        disp('Monotonically increasing CO2, readjusting feasiblity error.');
        ceq = [-30000,-30000];
        
    end
    
    it_count = it_count+1; % Advance counter
    
end

