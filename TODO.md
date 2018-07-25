# TODO

+ make functions to read data into long format. We need each row to represent a single measurement, e.g.,

O3 | FEV | VE | time 
---|-----|----|-----
0 | 4.5 | 7.2 | 120

Where,

  - O3: the ozone 
  - FEV: the FEV measured at time t
  - VE: the VE at time t
  - time: the time of the measurement in minutes
  
To do this, we need to reference data on multiple sheets of the xls files. We also need to pull and cross reference the experiment to convert the measurement (e.g., FEV0_1) into the minutes since exposure.

+ Simulate data and ensure the fitting algorithm can recover known parameters
