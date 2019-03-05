# Testing the Williams model against "known" values
library(rstan)
mod = stan_model("src/williams.stan")

source("scripts/munge_data_williams.R")


init_list = list(B = c(1.475, -0.4847, 0.0068, 49.2633, 0.0042, 1.0, 0, 0.0001, 10),
                 U = rep(0.01, 540)) 
ans$sigma_U = 1
fit = optimizing(mod, data = ans, init = init_list, verbose = TRUE)

