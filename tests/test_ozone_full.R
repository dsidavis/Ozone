# Testing Ed's model fitted to the full dataset
library(rstan)
source("scripts/munge_data_williams.R")

mod = stan_model("src/ozone_all.stan")

inits = list(dos = 1100, a = -0.001, k = 0.02, sigma = 0.1)
fit = optimizing(mod, data = ans, init = inits, verbose = TRUE)

