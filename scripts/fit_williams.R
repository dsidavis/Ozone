# Fit the model for McDonnel et al

# SAS code for model:
library(rstan)
mod = stan_model("src/williams.stan")
ans$sigma_U = 1
optimizing(mod, data = ans)
