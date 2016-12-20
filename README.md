# Replication code for: "Steering the Climate System: Using Inertia to Lower the Cost of Policy"

This is the main file to numerically solve the temperature-targeting framework in Lemoine and Rudik (2016). The first part of the file initializes parameters and switches. The second part solves the Hotelling / CO2 target model with stationary and non-stationary emissions and abatement cost. The next part solves the inertia / temperature target model with stationary and non-stationary emissions and abatement cost. The final part solves the model with the GHKT decay structure. The model is a system of ordinary differential equations that we opt to solve using a reverse-shooting style method which increases stability of the solvers. Depending on which setting we are analyzing, the ode solver finds optimal trajectories conditional on being passed a value of:
  1) The temperature shadow cost at the time the target is hit
  2) The 'start time' of the non-stationary processes which is not
     necessarily 0
  3) Both via a solver
  4) All of the above and initial time via a solver
By searching for the 'start time' of exogenous processes we are effectively searching for the time at which we hit the target. All other states/variables are pinned down by the efficient conditions to hit a target. Given a returned solution trajectory based off the passed shadow values and/or 'start time,' we find the numerically solved initial time of the problem by searching for when the CO2 trajectory is equal to its known, exogenously given initial value in 2005. The solver then maximizes an objective function that is a constant 0 subject to the constraint(s):
  1) Exogenous initial temperature - temperature at numerically solved initial time = 0
  2) Exogenous initial emissions/gdp - initial emissions/gdp at numerically solved initial time = 0 
  3) Both 
Once these conditions are satisfied, we have the temperature shadow cost at the time we arrive at the target and/or the correct time we arrive at the constraint and therefore have pinned down all state and co-state trajectories and have fully solved the problem. Note that at these very fine tolerances the ODE solver will throw warnings from not being able to take a small enough stepsize. These can be ignored. We opt for the warnings with a low tolerance instead of increasing the tolerance solely to avoid the warnings.
