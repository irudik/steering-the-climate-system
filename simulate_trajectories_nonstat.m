function [y] = simulate_trajectories_nonstat(terminal_search)

%% Simulate non-stationary trajectories using the model solution,
% it is the same as parameter_search functions but returns the trajectories
% instead of constraints

global params timespace it_count logic

% Call the solver for the passed initial time for non-stationary processes
[t1,y1,te1,ye1,ie1] = ode_nonstat_to_constraint(terminal_search);

if logic.hotelling
    y = y1;
    
else
        
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
    
    
    
    if ie1 == 2
        
        % Hit initial time
        disp(['Found initial value with ' num2str(terminal_search(2)-te1) ' years from there until guessed steady-state-like conditions.']);
        y = y1;
        t = t1;
        ye3 = ye1;
        
    else % ie1==1, so found time with full abatement
        
        disp(['Negative emission constraint binds at state ' mat2str(ye1,4) ', with ' num2str(terminal_search(2)-te1,4) ' years from there until guessed steady-state-like conditions.']);
        
        
        %% Use those results to find first time that negative emission constraint binds or hit initial temperature
        ye1 = [ye1 te1];
        [t2,y2,te2,ye2,ie2] = ode_nonstat_constrained(ye1);
        
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
        
        t2=t2+(t1(end)-terminal_search(2));
        
        if ie2 == 2
            % found time with initial temperature, so store output and skip next
            % simulations
            disp(['Found initial value with ' num2str( terminal_search(2)-te2) ' years from there until constraint stopped binding.']);
            y = vertcat(y1,y2);
            t = vertcat(t1,t2);
            ye3 = ye2;
            
        else % ie2 == 1, so found time with full abatement
            
            
            disp(['Negative emission constraint first binds at state ' mat2str(ye2,4) ' and then binds for ' num2str( terminal_search(2)-te2,4) ' years.']);
            
            
            %% Now simulate backwards to initial temperature
            ye2(5) = te2;
            [t3,y3,te3,ye3,ie3] = ode_nonstat_to_init(ye2);
            
            if all(y3(end,:)==y3(end-1,:))
                y3=y3(1:end-1,:);
            end
            
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
            
            disp(['Negative emission constraint first binds after ' mat2str(terminal_search(2)-te3) ' years.']);
            
            t3=t3+(t2(end)-terminal_search(2));
            
            y = vertcat(y1(:,1:7),y2(:,1:7),y3(:,1:7));
            t = vertcat(t1,t2,t3);
            
        end
        
        
    end
   
    
    
end

