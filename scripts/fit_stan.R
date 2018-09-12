# Fit the data from a single sheet
library(rstan)
source("R/fit_funs.R")
source("R/model_functions.R")
source("R/readUCD.R")
f = "data/UCD_WCA2002.xls"

d = readExp(f)

d = lapply(d, deltaFEV1)

mod = stan_model("src/ozone.stan")



tt = d[[2]]

data_list = list(n_blocks = nrow(tt),
                 n_fev1 = sum(!is.na(tt$dFEV1)),
                 o3 = tt$O3,
                 ve = tt$VE,
                 t_stop = tt$endTime,
                 fev1_pts = tt$endTime[!is.na(tt$dFEV1)],
                 fev1 = tt$dFEV1[!is.na(tt$dFEV1)])

optimizing(mod, data = data_list)
