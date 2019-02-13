# Fit the data from a single sheet
library(rstan)
source("R/fit_funs.R")
source("R/model_functions.R")
source("R/readUCD.R")
f = "data/UCD_WCA2002.xls"

d = readExp(f)

d = lapply(d, deltaFEV1)

mod = stan_model("src/ozone.stan")

ind_data = list(O3 = rep(rep(c(0.123, 0), each = 2), 3),
                Ve = rep(c(30, 13), length.out = 12), 
                t_stop = c(50, 63, 113, 126, 176, 200, 260, 273, 323, 336, 386, 399))

optimizing(mod, data = list(n_blocks = 12,
                            n_fev1 = 12,
                            o3 = ind_data$O3,
                            ve = ind_data$Ve,
                            t_stop = ind_data$t_stop,
                            fev1_pts = ind_data$t_stop,
                            fev1

fits = lapply(seq(along = d), function(j){ 
try({
    tt = d[[j]]
    i = which(is.na(tt$VE))
    tt$VE[i] = tt$VE[i-1]
    data_list = list(n_blocks = nrow(tt),
                     n_fev1 = sum(!is.na(tt$dFEV1)),
                     o3 = tt$O3,
                     ve = tt$VE,
                     t_stop = tt$endTime,
                     fev1_pts = tt$endTime[!is.na(tt$dFEV1)],
                     fev1 = tt$dFEV1[!is.na(tt$dFEV1)] * 100)


# fit = sampling(mod, data = data_list, iter = 200)

    optimizing(mod, data = data_list, init = list(dos = 1000, k = .6, a = -0.001))
})
})
data_list
