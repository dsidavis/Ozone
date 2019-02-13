# Functions for simulating data according to each model

source("R/model_functions.R")
# source("R/sim_funs.R")

# Using z4901_HR7PRED.xls as an example

time = c(55, 115, 180, 270, 330, 395)

o3 = rep(c(0.001, 0.08, 0.1), 3)

ve = rnorm(n = length(time), mean = 40, sd = 5)

dos = 1050 # mean of m/f from paper


experimentFEV1(o3, ve, time, Dos = dos, K = 0.69, A = -0.013)

# Attempt to duplicate the Excel spreadsheet - "Pred modV5.xlsx"
# Cols AC - AH
#

time = c(0,50, 60, 110, 120, 170, 215, 265, 275, 325, 335, 385, 399)
         
o3 = c(0.0005, 0.0005, 0.0004, 0.0004, 0.0005, 0.0005, 0.0005, 0.0005,
       0.0006, 0.0006, 0.0005, 0.0005)

ve = rep(c(38.3, 15.0792453), 6)

experimentFEV1(o3, ve, time, Dos = 1400, A = -0.0246633, K = 35/0.693)
# Appears to give exact same numbers
# Now one with change in dFEV1

o3 = c(0.0793933333, 0.0793933333, 0.07981, 0.07981,
       0.08029, 0.08029,
       0.07958, 0.07958,
       0.0803966666666, 0.0803966666666,
       0.07977333333333, 0.07977333333333) 

ans = experimentFEV1(o3, ve, time, Dos = 1400, A = -0.0246633, K = log(2)/35)

expected_ans = c(0, 0, 0, 0, 0, 0, -1.72087088354959e-94, -1.90017077526133e-49, 
-4.43485912124668, -4.19394040921599, -6.22098178707323, -5.47539345520911
)

stopifnot(all.equal(ans,expected_ans))
