# Testing Ed's model fitted to the full dataset
library(rstan)
source("scripts/munge_data_williams.R")

mod = stan_model("src/ozone_all.stan")


ans$dos = 900
ans$k = 0.002
ans$a = -0.002
ans$sigma = 5

inits = list(dos = 1100, a = -0.001, k = 0.02, sigma = 0.1)
optimizing(mod, data = ans, verbose = TRUE)


