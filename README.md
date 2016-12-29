# "Steering the Climate System: Using Inertia to Lower the Cost of Policy"
## By: Derek Lemoine and Ivan Rudik

This code is written in MATLAB and requires the KNITRO solver. As is, the model replicates results found in the paper but parameters can be changed to explore other settings. Initial guesses for the solver are set near the solution to reduce run time.

### Instructions:
The main file is `MAIN_SCRIPT.m`. At the top of the file, the user must indicate the location of the KNITRO options file using `options_file`. Next, the user selects the type of model to run using the `logic` structure. Setting `logic.stat_ems = 0` will run the model with non-stationary emissions. Setting `logic.ghkt = 1` will alter the carbon decay structure to match that of Golosov, Hassler, Krusell, and Tsyvinski (2014). `logic.hotelling` determines whether, for a given model, you will perform a run accounting for inertia in warming `logic.hotelling = 0` or not `logic.hotelling = 1`. Note two things:

1. To perform a non-stationary emissions run, you must use the base decay structure, and to use the GHKT decay structure, you must have stationary emissions. If you attempt to run both, the code will default to a GHKT model with stationary emissions.

2. You can only use the code to automatically vary discounting and the level of inertia for the base model in the main text. However, alternative parameterizations for the GHKT or non-stationary emissions models can be explored by hard coding changes to the parameter values.

Finally, the user must determine which temperature targets to loop over and solve the model. `params.first_temp` and `params.final_temp` determine the lowest and highest temperature target to use, and `params.temp_incr` determines the size of the increment in the temperature target as we loop from the lowest to the highest target.
