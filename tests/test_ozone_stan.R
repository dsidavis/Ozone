# Testing the Stan version of Ed's model

source("R/model_functions.R")
mod = stan_model("src/ozone.stan")

time = c(50, 60, 110, 120, 170, 215, 265, 275, 325, 335, 385, 399)
ve = rep(c(38.3, 15.0792453), 6)
o3 = c(0.0793933333, 0.0793933333, 0.07981, 0.07981,
       0.08029, 0.08029,
       0.07958, 0.07958,
       0.0803966666666, 0.0803966666666,
       0.07977333333333, 0.07977333333333) 

expected_ans = c(0, 0, 0, 0, 0, 0, -1.72087088354959e-94, -1.90017077526133e-49, 
-4.43485912124668, -4.19394040921599, -6.22098178707323, -5.47539345520911
)

data_list = list(n_blocks = length(o3),
                 n_fev1 = length(expected_ans),
                 o3 = o3,
                 ve = ve,
                 t_stop = time,
                 fev1_pts = time,
                 fev1 = expected_ans)

optimizing(mod, data = data_list)
