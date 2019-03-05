# Fit the model for McDonnel et al

# SAS code for model:
library(rstan)
mod = stan_model("src/williams.stan")


source("scripts/munge_data_williams.R")


init_list = list(B = c(1.475, -0.4847, 0.0068, 49.2633, 0.0042, 1.0, 0, 0.0001, 10),
                 U = rep(0.01, 540)) 
ans$sigma_U = 1
# b = replicate(100, system.time({fit = optimizing(mod, data = ans, init = init_list)}))
fit = optimizing(mod, data = ans, init = init_list, verbose = TRUE)

fit = sampling(mod, data = ans, cores = 4)

b = replicate(100, optimizing(mod, data = ans)$par[1:9])
summary(b)
b = replicate(100, optimizing(mod, data = ans)$value)

fit = optimizing(mod, data = ans, algorithm = "Fixed_param", init = init_list)


mod2 = 
