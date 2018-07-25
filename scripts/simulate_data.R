# Functions for simulating data according to each model

source("R/model_functions.R")
# source("R/sim_funs.R")

# Using z4901_HR7PRED.xls as an example

time = c(55, 115, 180, 270, 330, 395)

o3 = c(0, 0.08, 0.1)

ve = rnorm(n = length(time), mean = 40, sd = 5)

dos = 1050 # mean of m/f from paper

uos = UOS(o3[2], ve, time, Dos = dos)

FEV(UOS = uos, K = 0.2, A = -0.02, fev_base = 4.6)
