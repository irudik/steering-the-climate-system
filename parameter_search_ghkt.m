function [cineq, ceq] = parameter_search_ghkt(terminal_search)

%% Search for terminal conditions that solve the model

% GHKT model: search for M1(tau) for the Hotelling model, and M1(tau)
% and lambda_T(tau) for the least-cost model such that
% Hotelling: Both CO2 stocks at t0 are equal to their initial levels
% Least-cost: Both CO2 stocks and temperature at t0 are equal to initial
% values, allow for negative emissions constraint to bind

global params timespace it_count logic

cineq = 0; % Automatically satisfy inequality constraints

if logic.hotelling
    disp(['GHKT hotelling run ' num2str(it_count) '. Terminal conditions: ' mat2str(terminal_search,6) '.']);
else
    disp(['GHKT least-cost run ' num2str(it_count) '. Terminal conditions: ' mat2str(terminal_search,6) '.']);
end

%% Find last time that negative emission constraint binds
[t1,y1,te1,ye1,ie1] = ode_ghkt_to_constraint(terminal_search);

% If we cannot find the last time, change the time horizon
if isempty(ye1)
    params.tau = params.tau*2;
    [t1,y1,te1,ye1,ie1] = ode_ghkt_to_constraint(terminal_search);
    if isempty(ye1)
        params.tau = params.tau/4;
        [t1,y1,te1,ye1,ie1] = ode_ghkt_to_constraint(terminal_search);
    end
    if isempty(ye1)
        disp('Never found time at which either constraint binds or hit initial temperature.');
        ceq = [-30000];
        return;
    end
end

if size(ie1,1)>1 % first instant may have had an event
    ie1 = ie1(end,:);
    ye1 = ye1(end,:);
    te1 = te1(end,:);
end

if ie1==2
    % found time with temperature at initial temperature, so store output
    % and skip next simulations
    disp(['Found initial value with ' num2str(params.tau-te1,4) ' years from there until guessed steady-state-like conditions.']);
    y = y1;
    t = t1;
    ye3 = ye1;
    
else % ie1 == 1, negative emissions constraint is binding so we must
     % simulate the model where abatement is constrained
    
    disp(['Negative emission constraint binds at state '...
        mat2str(ye1,4) ', with ' num2str(params.tau-te1,4)...
        ' years from there until guessed steady-state-like conditions.']);
    
    %% Use the simulated results for after the constraint has stopped binding to find the
    % first time that negative emission constraint starts binding, or we hit initial temperature
    [t2,y2,te2,ye2,ie2] = ode_ghkt_constrained(ye1);
    
    % If we cannot find the first time, change the time horizon
    if isempty(ye2)
        params.tau = params.tau*4;
        [t2,y2,te2,ye2,ie2] = ode_ghkt_constrained(ye1);
        if isempty(ye2)
            params.tau = params.tau/16;
            [t2,y2,te2,ye2,ie2] = ode_ghkt_constrained(ye1);
        end
        if isempty(ye2)
            disp('Never found first time at which constraint started binding.');
            ceq = [-30000];
            return;
        end
    end
        
    if size(ie2,1)>1 % first instant should have had an event
        ie2 = ie2(end,:);
        ye2 = ye2(end,:);
        te2 = te2(end,:);
    end
    
    % New time
    t2=t2+(t1(end)-params.tau);
    
    if ie2 == 2
        % found time with initial temperature, so store output and skip next
        % simulations
        disp(['Found initial value with ' num2str(params.tau-te2,4) ' years from there until constraint stopped binding.']);
        y = vertcat(y1(:,1:6),y2(:,1:6));
        t = vertcat(t1,t2);
        ye3 = ye2;
        
    else % ie2==1, so found first time with full abatement
        
        disp(['Negative emission constraint first binds at state '...
            mat2str(ye2,4) ' and then binds for '...
            num2str(params.tau-te2,4) ' years.']);
        
        %% Now simulate backwards to the initial time
        [t3,y3,te3,ye3] = ode_ghkt_to_init(ye2);
        
        % If we cannot find the first time, change the time horizon
        if isempty(ye3)
            params.tau = params.tau*4;
            [t3,y3,te3,ye3] = ode_ghkt_to_init(ye2);
            if isempty(ye3)
                params.tau = params.tau/16;
                [t3,y3,te3,ye3] = ode_ghkt_to_init(ye2);
            end
            if isempty(ye3)
                disp('Never found initial value.');
                ceq = [-30000];
                return;
            end
        end
        
        disp(['Negative emission constraint first binds after '...
            mat2str(params.tau-te3,4) ' years.']);
        
        t3=t3+(t2(end)-params.tau);
        
        y = vertcat(y1(:,1:6),y2(:,1:6),y3(:,1:6));
        t = vertcat(t1,t2,t3);
                
    end
    
end

% Add emissions and time to matrix
if logic.hotelling
    y(:,5) = params.initial_emissions;
    y(:,6) = t;
    y(:,7) = y(:,2)+y(:,4);
    y(:,8) = y(:,7)-params.terminal_co2;
else
    y(:,7) = params.initial_emissions;
    y(:,8) = t;
    y(:,9) = y(:,4)+y(:,6);
    y(:,10) = y(:,9)-params.terminal_co2;
end

count = 0;

% Set equality constraint to be M1(t0) = initial M1, M2(t0) = initial M2,
% and for model of least-cost policy, T(t0) = initial temp
if logic.hotelling
    ceq = [ye3(2)-params.initial_co2_perm;
    ye3(4)-params.initial_co2_geo];
else
ceq = [ye3(3)-params.initial_temp
    ye3(4)-params.initial_co2_perm;
    ye3(6)-params.initial_co2_geo];
end

% Return error in initial CO2 and temperature
if ~isempty(ceq)
    if logic.hotelling
    disp(['True initial CO2 perm: ' num2str(params.initial_co2_perm,4), ' simulated initial CO2 perm: ' num2str(ye3(2))]);
    disp(['True initial CO2 geo: ' num2str(params.initial_co2_geo,4), ' simulated initial CO2 geo: ' num2str(ye3(4))]);
    else
        disp(['True initial temp: ' num2str(params.initial_temp,4), ' simulated initial temp: ' num2str(ye3(3))]);
        disp(['True initial CO2 perm: ' num2str(params.initial_co2_perm,4), ' simulated initial CO2 perm: ' num2str(ye3(4))]);
        disp(['True initial CO2 geo: ' num2str(params.initial_co2_geo,4), ' simulated initial CO2 geo: ' num2str(ye3(6))]);    
    end
else
    disp('Monotonically increasing CO2, readjusting feasiblity error.');
    ceq = [-30000];
end

it_count = it_count+1; % Advance counter
