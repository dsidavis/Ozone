# TODO

+ Update C code (src/deltaX.c, src/timeLoop.c) to reflect changes in R/model_functions.R

+ clean up the existing code in scripts

+ document the functions

+ turn repo into package

+ add code to compile the models on package install
 
+ write a vignette for Ed/others to get started

+ [done] update testing scripts to reflect changes in function args

+ [done] Check simultated data for similarity to observed data

+ Resolve issues with fit not resolving to unique solution
  - No easy fix. Need to do a grid search (Ed's model) or random inits and focus on most optimum values (Williams)

+ [done] make functions to read data into long format. We need each row to represent a single measurement, e.g.,

O3 | FEV | VE | time 
---|-----|----|-----
0 | 4.5 | 7.2 | 120

Where,

  - O3: the ozone 
  - FEV: the FEV measured at time t
  - VE: the VE at time t
  - time: the time of the measurement in minutes
  
To do this, we need to reference data on multiple sheets of the xls files. We also need to pull and cross reference the experiment to convert the measurement (e.g., FEV0_1) into the minutes since exposure.

+ [done] Simulate data and ensure the fitting algorithm can recover known parameters
 
