# "Steering the Climate System: Using Inertia to Lower the Cost of Policy"
## By: Derek Lemoine and Ivan Rudik

This code is written in MATLAB and requires the KNITRO solver. As is, the model replicates results found in the paper but parameters can be changed to explore other settings. Initial guesses for the solver are set near the solution to reduce run time.

### Instructions:
The main file is `MAIN_SCRIPT.m`. At the top of the file, the user must indicate the location of the KNITRO options file using `options_file`. Next, the user selects the type of model to run using the `logic` structure. Setting `logic.stat_ems = 0` will run the model with non-stationary emissions following the DICE-2007 model. Setting `logic.ghkt = 1` will alter the carbon decay structure to match that of Golosov, Hassler, Krusell, and Tsyvinski (2014). `logic.hotelling` determines whether, for a given model, you will perform a run accounting for inertia in warming `logic.hotelling = 0` or not `logic.hotelling = 1`. When running the base model, setting `logic.vary_params = 1` will also solve the model for an inertia parameter that is 25% higher and lower, and then for a discount rate of 1.4% from Stern (2007). Note two things:

1. To perform a non-stationary emissions run, you must use the base decay structure, and to use the GHKT decay structure, you must have stationary emissions. If you attempt to run both, the code will default to a GHKT model with stationary emissions.

2. You can only use the code to automatically vary discounting and the level of inertia for the base model in the main text. However, alternative parameterizations for the GHKT or non-stationary emissions models can be explored by hard coding changes to the parameter values.

Finally, the user must determine which temperature targets to loop over and solve the model. `params.first_temp` and `params.final_temp` determine the lowest and highest temperature target to use, and `params.temp_incr` determines the size of the increment in the temperature target as we loop from the lowest to the highest target.

### How the algorithms work
We solve all the models using reverse-shooting style algorithms (see Judd (1998)). We use the necessary conditions to pin down as many terminal conditions as possible and then search over the free terminal conditions that yield the correct initial conditions at initial time `t0`. Here we take the terminal time to be the time at which the policymaker hits the temperature constraint, after the constraint is hit the policymaker simply abates to maintain the constraint. As with most optimal control problems, the problem is significantly more stable when run in reverse time so we solve the problem backwards.

#### Base model (main text)
For the base Hotelling model, we do not need to search over free terminal conditions. Here we have one CO2 state and one CO2 costate. The temperature target pins down the CO2 state at the terminal time. The target also implies that CO2 cannot be changing at the terminal time, so the velocity of CO2 must be zero, which implies the level of abatement at the terminal time must be the level that maintains CO2 equal to the target. This level of abatement then implies the CO2 costate from the maximality condition. All the terminal conditions are accounted for and the system can simply be simulated in reverse.

The base inertia model follows the same logic as the base Hotelling model. The target pins down terminal temperature and CO2, when then pins down terminal abatement and the terminal CO2 costate. however the temperature costate at the terminal time is free. Therefore, we search over the terminal temperature costate with the objective equating the initial temperature parameter to the simulated temperature at the time of the simulation where the simulated CO2 is within a given tolerance of the given initial CO2 stock. 

#### GHKT decay model (appendix)


#### Non-stationary emissions model (appendix)
