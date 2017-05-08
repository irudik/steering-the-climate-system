function [y] = simulate_trajectories_ghkt(terminal_search)

%% Simulate GHKT trajectories using the model solution,
% it is the same as parameter_search functions but returns the trajectories
% instead of constraints

global params timespace it_count logic


%% Find last time that negative emission constraint binds
[t1,y1,te1,ye1,ie1] = ode_ghkt_to_constraint(terminal_search);

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
    % found time with initial temperature, so store output and skip next
    % simulations
    disp(['Found initial value with ' num2str(params.tau-te1) ' years from there until guessed steady-state-like conditions.']);
    y = y1;
    t = t1;
    ye3 = ye1;
    
else % ie1==1, so found time with full abatement
    
    disp(['Negative emission constraint binds at state ' mat2str(ye1) ', with ' num2str(params.tau-te1) ' years from there until guessed steady-state-like conditions.']);
    
    
    %% Use those results to find first time that negative emission constraint binds or hit initial temperature
    
    [t2,y2,te2,ye2,ie2] = ode_ghkt_constrained(ye1);
    
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
    
    t2=t2+(t1(end)-params.tau);
    
    if ie2==2
        % found time with initial temperature, so store output and skip next
        % simulations
        disp(['Found initial value with ' num2str(params.tau-te2) ' years from there until constraint stopped binding.']);
        y = vertcat(y1(:,1:6),y2(:,1:6));
        t = vertcat(t1,t2);
        ye3 = ye2;
        
    else % ie2==1, so found time with full abatement
        
        
        disp(['Negative emission constraint first binds at state ' mat2str(ye2) ' and then binds for ' num2str(params.tau-te2) ' years.']);
        
        
        %% Now simulate backwards to initial temperature
        
        [t3,y3,te3,ye3] = ode_ghkt_to_init(ye2);
        
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
        
        disp(['Negative emission constraint first binds after ' mat2str(params.tau-te3) ' years.']);
        
        t3=t3+(t2(end)-params.tau);
        
        y = vertcat(y1(:,1:6),y2(:,1:6),y3(:,1:6));
        t = vertcat(t1,t2,t3);
                
    end
    
end

% add emissions and time to matrix
if logic.hotelling
    y = [zeros(size(y(:,1))) y(:,1) zeros(size(y(:,1))) y(:,2:end)];
    y(:,7) = params.initial_emissions;
    y(:,8) = t;
    y(:,9) = y(:,4)+y(:,6);
    y(:,10) = y(:,9)-params.terminal_co2;
else
    y(:,7) = params.initial_emissions;
    y(:,8) = t;
    y(:,9) = y(:,4)+y(:,6);
    y(:,10) = y(:,9)-params.terminal_co2;
end