% Set directories for saving
root_dir = pwd;
if ~logic.stat_ems
    store_dir = [root_dir '\results\nonstat'];
else
    store_dir = [root_dir '\results\base'];
end



% If necessary, extend the simulation so that it reaches
% params.time_horizon
% Time vector up until the end of the simulated time horizon
temp_time = [(fliplr((inertia_out(1,7)+params.grid_res:params.grid_res:params.time_horizon)))'];

% How long we need to extend the simulated horizon to reach the full
% time horizon we selected
time_extend = length(temp_time);

if inertia_out(1,7) < params.time_horizon
    
    inertia_out = [zeros(size(temp_time,1),size(inertia_out,2)); inertia_out];
    
    inertia_out(1:size(temp_time,1),7) = temp_time;
    
    % extend temperature, co2, emissions, abatement for stationary emissions
    if logic.stat_ems
        
        inertia_out(1:time_extend,[4]) = inertia_out(time_extend+1,[4]); % stay at same CO2, emissions once we hit the target, approximating the path
        inertia_out(1:time_extend,[5]) = inertia_out(time_extend+1,[5]); % stay at same CO2, emissions once we hit the target, approximating the path
        inertia_out(1:time_extend,[6]) = inertia_out(1:time_extend,5) - params.delta*(inertia_out(1:time_extend,4)-params.mpre); % abatement to keep mdot=0
        
        if ~logic.hotelling
            for i=1:time_extend % Simulate resulting temperature
                inertia_out(time_extend-i+1,3) = (1-params.grid_res*params.phi)*inertia_out(time_extend-i+2,3) +...
                    params.grid_res*params.phi*params.s*params.alpha*log(inertia_out(time_extend-i+2,4)/params.mpre);
            end
        end
        
    % Adjust abatement/emissions beyond epsilon ball of the target if non-stationary emissions
    elseif ~logic.stat_ems
        
        inertia_out(1:time_extend,[4]) = inertia_out(time_extend+1,[4]); % stay at same CO2, emissions once we hit the target, approximating the path
        inertia_out(:,5) = params.initial_emissions*exp(.0068*(inertia_out(:,7))); % emissions trajectory
        inertia_out(1:time_extend,6) = inertia_out(1:time_extend,5) - params.delta*(inertia_out(1:time_extend,4)-params.mpre); % abatement to keep mdot=0
        
        inertia_out(inertia_out(:,6)>inertia_out(:,5),6) = inertia_out(inertia_out(:,6)>inertia_out(:,5),5);

        for i=1:time_extend % Simulate resulting temperature
            inertia_out(time_extend-i+1,3) = (1-params.grid_res*params.phi)*inertia_out(time_extend-i+2,3) +...
                params.grid_res*params.phi*params.s*params.alpha*log(inertia_out(time_extend-i+2,4)/params.mpre);
        end
    end
end

% Abatement cost to maintain Adot = 0
inertia_out(1:end,11) = params.initial_gdp*params.a_0*params.sigma_0*(inertia_out(1:end,6)./inertia_out(1:end,5)).^params.a_2/params.a_2;

% Extend temperature shadow cost trajectories
if ~logic.hotelling
    for i=1:time_extend
        inertia_out(time_extend-i+1,1) = inertia_out(time_extend-i+2,1) + params.grid_res*(params.r+params.phi)*inertia_out(time_extend-i+2,1);
    end
end

% Construct shadow cost decomposition
if ~logic.hotelling
    
    % Hotelling component
    inertia_out(:,8) = inertia_out(end,2).*exp((params.r+params.delta).*inertia_out(:,7));
    
    % inertia component
    time_vec = flipud(inertia_out(:,7)); % flip time vector for computing integral
    m_vec = flipud(inertia_out(:,4)); % flip M vector
    
    % Integral term
    clear integral
    for time = 2:size(inertia_out,1) % integrate data, skip initial time since cannot integrate over zero distance
        expo = exp(-(params.phi-params.delta)*(time_vec(time)-time_vec(1:time))); % vector of the exponential
        integral(time-1,:) = trapz(time_vec(1:time),expo(1:time)*params.phi*params.s*params.alpha./m_vec(1:time)); % approximate the integral
    end
    
    integral = [0; integral];
    integral = flipud(integral);
    
    % Inertia component
    inertia_out(:,9) = inertia_out(:,1).*integral;
    
    % Calculated co2 shadow cost
    inertia_out(:,10) = inertia_out(:,8)-inertia_out(:,9); 
end

% Calculate CO2 shadow cost through marginal abatement cost
if logic.stat_ems
    inertia_out(1:time_extend,2) = params.initial_gdp*params.a_0*params.sigma_0*inertia_out(1:time_extend,6).^(params.a_2-1)./(inertia_out(1:time_extend,5).^params.a_2);
else
    inertia_out(1:time_extend,2) = params.initial_gdp*params.a_0*params.sigma_0*inertia_out(1:time_extend,6).^(params.a_2-1)./(inertia_out(1:time_extend,5).^params.a_2);
end

%Calculate temperature under the hotelling policy
if logic.hotelling
    
    inertia_out(end,3) = params.initial_temp;
    
    for i = 1:size(inertia_out,1)-1
        
        % discretize analytic temperature transition equation
        inertia_out(size(inertia_out,1)-i,3) = params.grid_res*params.phi*params.s*params.alpha*...
            log(inertia_out(size(inertia_out,1)-i+1,4)/params.mpre) + (1-params.grid_res*params.phi)*inertia_out(size(inertia_out,1)-i+1,3);
        
    end
end

% Present value of cumulative abatement cost
inertia_out(:,12) = exp(-params.r*inertia_out(:,7)).*inertia_out(:,11);

% flip matrices to be in correct order
inertia_out = flipud(inertia_out);

% Store in structure

% Temperature shadow cost
output.temp_shadow = inertia_out(:,1);

% CO2 shadow cost
output.co2_shadow = inertia_out(:,2);

% Temperature
output.temp = inertia_out(:,3);

% CO2
output.co2 = inertia_out(:,4);

% Emissions
output.emissions = inertia_out(:,5);

% Abatement
output.abatement = inertia_out(:,6);

% Time
output.time = inertia_out(:,7);

% Abatement cost
output.abatement_cost = 1000*inertia_out(:,11);

% Abatement cost first 200 years
output.abatement_cost_200 = 1000*inertia_out(1:200/params.grid_res,11);

% Present value of abatement cost
output.abatement_cost_pv = 1000*inertia_out(:,12);

% Present value of abatement cost first 200 years
output.abatement_cost_pv_200 = 1000*exp(-params.r*inertia_out(1:200/params.grid_res,7)).*inertia_out(1:200/params.grid_res,11);

% Total abatement cost first 200 years
output.total_abate_cost_200 = trapz(output.time(1:200/params.grid_res),output.abatement_cost_pv_200);

% Present value of total abatement cost
output.total_abate_cost = trapz(output.time,output.abatement_cost_pv);

% Hotelling component of CO2 shadow cost
output.hotelling_comp = inertia_out(:,8);

% Inertia component of CO2 shadow cost
output.inertia_comp = inertia_out(:,9);

% Implied CO2 shadow cost from both components
output.implied_co2_shadow = inertia_out(:,10);

% CO2 in parts per million
output.co2_ppm = output.co2/2.16;

% CO2 shadow price in $/tCO2
output.co2_shadow_perton = (1000*12/44)*output.co2_shadow;

% Hotelling component of CO2 shadow price in $/tCO2
output.co2_shadow_hotelling_perton = (1000*12/44)*output.hotelling_comp;

% Inertia component of CO2 shadow price in $/tCO2
output.co2_shadow_inertia_perton = (1000*12/44)*output.inertia_comp;

% Efficient carbon tax in $/tCO2
output.tax_perton = (1000*12/44)*params.a_0*params.sigma_0*...
    params.initial_gdp*(output.abatement).^(params.a_2-1)./(output.emissions.^params.a_2);

% Total abatement first 200 years
output.total_abatement = trapz(output.time(1:200/params.grid_res),output.abatement(1:200/params.grid_res));

% Store workspace in run-specific directory
cd(store_dir)
if logic.vary_params
    save(['target_' num2str(params.terminal_temp)...
        '_hotelling_' num2str(logic.hotelling) '_r_' num2str(fields_r{rs}) '_inertia_' num2str(fields_phi{phis}) '.mat'],'output');
else
    save(['target_' num2str(params.terminal_temp)...
        '_hotelling_' num2str(logic.hotelling) '.mat'],'output');
end
cd(root_dir)
