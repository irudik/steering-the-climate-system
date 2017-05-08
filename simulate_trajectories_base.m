function [y] = simulate_trajectories_base(terminal_search)

%% Search for terminal conditions that solve the model

% Base model: search for lambda_T(tau) that results in T(t') and M(t')
% equal to the initial temperature and CO2 values at some time t'

global params timespace it_count logic

cineq = 0; % Automatically satisfy inequality constraints

disp(['Base least-cost run ' num2str(it_count) '. Terminal conditions: ' mat2str(terminal_search,6) '.']);

%% Find last time that negative emission constraint binds
[t1,y1,te1,ye1,ie1] = ode_base_to_constraint(terminal_search);

% If we cannot find the last time, change the time horizon
if isempty(ye1)
    params.tau = params.tau*2;
    [t1,y1,te1,ye1,ie1] = ode_base_to_constraint(terminal_search);
    if isempty(ye1)
        params.tau = params.tau/4;
        [t1,y1,te1,ye1,ie1] = ode_base_to_constraint(terminal_search);
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
    [t2,y2,te2,ye2,ie2] = ode_constrained_ghkt(ye1);
    
    % If we cannot find the first time, change the time horizon
    if isempty(ye2)
        params.tau = params.tau*4;
        [t2,y2,te2,ye2,ie2] = ode_base_constrained(ye1);
        if isempty(ye2)
            params.tau = params.tau/16;
            [t2,y2,te2,ye2,ie2] = ode_base_constrained(ye1);
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
        [t3,y3,te3,ye3] = ode_to_init_base(ye2);
        
        % If we cannot find the first time, change the time horizon
        if isempty(ye3)
            params.tau = params.tau*4;
            [t3,y3,te3,ye3] = ode_base_to_init(ye2);
            if isempty(ye3)
                params.tau = params.tau/16;
                [t3,y3,te3,ye3] = ode_base_to_init(ye2);
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



