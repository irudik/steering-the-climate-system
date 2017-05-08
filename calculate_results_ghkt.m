% Set directories for saving
root_dir = pwd;
store_dir = [root_dir '\results\ghkt\'];

% Abatement cost at each node
abatement = ((inertia_out(:,2)*params.psi_l+inertia_out(:,5)*params.psi_0*...
    (1-params.psi_l))*params.initial_emissions^params.a_2/...
    params.sigma_0/params.a_0/params.initial_gdp).^(1/(params.a_2-1));

abatement(abatement > params.initial_emissions,1) = params.initial_emissions;

% Abatement at time tau
abatement(1) = params.initial_emissions - ...
params.psi*(params.terminal_co2-inertia_out(1,4))/(params.psi_l+params.psi_0*(1-params.psi_l));

inertia_out(:,11) = abatement;

abatement_cost = params.sigma_0*params.a_0*params.initial_gdp*...
    (abatement./params.initial_emissions).^params.a_2;

% If necessary, extend the simulation so that it reaches
% params.time_horizon
if inertia_out(1,8) < params.time_horizon

    % Time vector up until the end of the simulated time horizon
    temp_time = [(fliplr((inertia_out(1,8)+params.grid_res:params.grid_res:params.time_horizon)))'];
    
    % How long we need to extend the simulated horizon to reach the full
    % time horizon we selected
    time_extend = length(temp_time);
    
    inertia_out = [zeros(size(temp_time,1),size(inertia_out,2)); inertia_out];
    
    inertia_out(1:size(temp_time,1),8) = temp_time;
    
    % Extend temperature, co2, emissions, abatement for stationary emissions
    inertia_out(1:time_extend,[9]) = inertia_out(time_extend+1,[9]); % stay at same CO2, emissions once we hit the target, approximating the path
    inertia_out(1:time_extend,[7]) = inertia_out(time_extend+1,[7]); % stay at same CO2, emissions once we hit the target, approximating the path
    
    % Abatement, M1, M2 to keep M1dot + M2dot = 0
    for t = time_extend:-1:1 
        inertia_out(t,11) = params.initial_emissions - ...
            params.psi*inertia_out(t+1,6)/(params.psi_l+params.psi_0*(1-params.psi_l));
        inertia_out(t,4) = inertia_out(t+1,4) + ...
            params.grid_res*params.psi_l*(params.initial_emissions-inertia_out(t,11));
        inertia_out(t,6) = inertia_out(t+1,6) + ...
            params.grid_res*(params.psi_0*(1-params.psi_l)*(params.initial_emissions-inertia_out(t,11))-params.psi*inertia_out(t+1,6));
    end
    
    if any(abs(inertia_out(:,4)+inertia_out(:,6) - inertia_out(:,9))>1e-6)
        error('Error in CO2 simulation');
    end
    
    for i=1:time_extend % Simulate resulting temperature
        inertia_out(time_extend-i+1,3) = (1-params.grid_res*params.phi)*inertia_out(time_extend-i+2,3) +...
            params.grid_res*params.phi*params.s*params.alpha*log((inertia_out(time_extend-i+2,4)+inertia_out(time_extend-i+2,6))/params.mpre);
    end
    
    % Extend temperature shadow costs
    for i=1:time_extend
        inertia_out(time_extend-i+1,1) = inertia_out(time_extend-i+2,1) + params.grid_res*(params.r+params.phi)*inertia_out(time_extend-i+2,1);
    end

        
end

% Carbon price (conditional on not binding) ($/tCO2)
co2_shadow_perton = 1000*12/44*params.a_0*params.sigma_0*params.initial_gdp*inertia_out(:,11).^(params.a_2-1)/params.initial_emissions^params.a_2;
inertia_out(:,12) = co2_shadow_perton;

% Abatement cost
abatement_cost = params.sigma_0*params.a_0*params.initial_gdp*(inertia_out(:,11)./params.initial_emissions).^params.a_2;
inertia_out(:,14) = abatement_cost;

% Construct shadow cost decomposition

% Hotelling component M1
inertia_out(:,15) = inertia_out(end,2).*exp((params.r).*inertia_out(:,8));

time_vec = flipud(inertia_out(:,8)); % flip time vector for computing integral
m_vec = flipud(inertia_out(:,4)+inertia_out(:,6)); % flip M vector

% Integral term M1
clear integral
for time = 2:size(inertia_out,1) % integrate data, skip initial time since cannot integrate over zero distance
    expo = exp(-(params.phi)*(time_vec(time)-time_vec(1:time))); % vector of the exponential
    integral(time-1,:) = trapz(time_vec(1:time),expo(1:time)*params.phi*params.s*params.alpha./m_vec(1:time)); % approximate the integral
end

integral = [0; integral];
integral = flipud(integral);

% inertia component
inertia_out(:,17) = inertia_out(:,1).*integral; 

% Calculated M1 shadow cost
inertia_out(:,18) = inertia_out(:,15)-inertia_out(:,17);

% Hotelling component M2
inertia_out(:,19) = inertia_out(end,5).*exp((params.r+params.psi).*inertia_out(:,8));

time_vec = flipud(inertia_out(:,8)); % flip time vector for computing integral
m_vec = flipud(inertia_out(:,4)+inertia_out(:,6)); % flip M vector

% Integral term M2
clear integral
for time = 2:size(inertia_out,1) % integrate data, skip initial time since cant integrate over zero distance
    expo = exp(-(params.phi-params.psi)*(time_vec(time)-time_vec(1:time))); % vector of the exponential
    integral(time-1,:) = trapz(time_vec(1:time),expo(1:time)*params.phi*params.s*params.alpha./m_vec(1:time)); % approximate the integral
end

integral = [0; integral];
integral = flipud(integral);

% Inertia component M2
inertia_out(:,20) = inertia_out(:,1).*integral;

% Present value of abatement cost trajectory
inertia_out(:,21) = exp(-params.r*inertia_out(:,8)).*inertia_out(:,14);

% flip matrices to be in correct order
inertia_out = flipud(inertia_out);

% Store in structure

% Temperature shadow cost
output.temp_shadow = inertia_out(:,1);

% M1 shadow cost
output.m1_shadow = inertia_out(:,2);

% Temeperature
output.temp = inertia_out(:,3);

% M1 CO2 stock
output.m1 = inertia_out(:,4);

% M2 shadow cost
output.m2_shadow = inertia_out(:,5);

% M2 CO2 stock
output.m2 = inertia_out(:,6);

% Emissions
output.emissions = inertia_out(:,7);

% Total CO2 stock
output.co2 = inertia_out(:,9);

% Abatement
output.abatement = inertia_out(:,11);

% Time
output.time = inertia_out(:,8);

% Abatement cost
output.abatement_cost = 1000*inertia_out(:,14);

% Abatement cost first 200 years
output.abatement_cost_200 = 1000*inertia_out(1:200/params.grid_res,11);

% Present value of abatement cost
output.abatement_cost_pv = 1000*inertia_out(:,21);

% Present value of abatement cost first 200 years
output.abatement_cost_pv_200 = 1000*exp(-params.r*inertia_out(1:200/params.grid_res,8)).*inertia_out(1:200/params.grid_res,14);

% Total abatement cost first 200 years
output.total_abate_cost_200 = trapz(output.time(1:200/params.grid_res),output.abatement_cost_pv_200);

% Present value of total abatement cost
output.total_abate_cost = trapz(output.time,output.abatement_cost_pv);

% Hotelling component of M1 shadow cost
output.hotelling_comp_m1 = inertia_out(:,15);

% Inertia component of M1 shadow cost
output.inertia_comp_m1 = inertia_out(:,17);

% Hotelling component of M2 shadow cost
output.hotelling_comp_m2 = inertia_out(:,19);

% Inertia component of M2 shadow cost
output.inertia_comp_m2 = inertia_out(:,20);

% CO2 in parts per million
output.co2_ppm = output.co2/2.16;

% CO2 in parts per million
output.m1_ppm = inertia_out(:,4)/2.16;

% CO2 in parts per million
output.m2_ppm = inertia_out(:,6)/2.16;

% CO2 shadow price in $/tCO2
output.co2_shadow_perton = inertia_out(:,12);

% Hotelling component of M1 shadow price in $/tCO2
output.co2_shadow_hotelling_m1_perton = (1000*12/44)*output.hotelling_comp_m1;

% Inertia component of M1 shadow price in $/tCO2
output.co2_shadow_inertia_m1_perton = (1000*12/44)*output.inertia_comp_m1;

% Hotelling component of M2 shadow price in $/tCO2
output.co2_shadow_hotelling_m2_perton = (1000*12/44)*output.hotelling_comp_m2;

% Inertia component of M2 shadow price in $/tCO2
output.co2_shadow_inertia_m2_perton = (1000*12/44)*output.inertia_comp_m2;

% Efficient carbon tax in $/tCO2
output.tax_perton = (1000*12/44)*params.a_0*params.sigma_0*...
    params.initial_gdp*(output.abatement).^(params.a_2-1)./(output.emissions.^params.a_2);

% M1 hotelling tax component in $/tCO2
output.tax_hotelling_m1_perton = (1000*12/44)*output.hotelling_comp_m1;

% M1 inertia tax component in $/tCO2
output.tax_inertia_m1_perton = (1000*12/44)*output.inertia_comp_m1;

% M2 hotelling tax component in $/tCO2
output.tax_hotelling_m2_perton = (1000*12/44)*output.hotelling_comp_m2;

% M2 inertia tax component in $/tCO2
output.tax_inertia_m2_perton = (1000*12/44)*output.inertia_comp_m2;

% Total abatement first 200 years
output.total_abatement = trapz(output.time(1:200/params.grid_res),output.abatement(1:200/params.grid_res));

% Store workspace in run-specific directory
cd(store_dir)
save(['target_' num2str(params.terminal_temp) '_hotelling_' num2str(logic.hotelling) '.mat'],'output');
cd(root_dir)
