# Fit the model for McDonnel et al

# SAS code for model can be found at the end of
# William F. McDonnell, Paul
# W. Stewart & Marjo V. Smith (2013)
# Ozone exposure-response model for
# lung function changes: an alternate variability structure,
# Inhalation Toxicology, 25:6, 348-353, DOI:
# 10.3109/08958378.2013.790523

library(rstan)
mod = stan_model("src/williams.stan")


source("scripts/munge_data_williams.R")

# List from original code
init_list = list(B = c(1.475, -0.4847, 0.0068, 49.2633, 0.0042, 1.0, 0, 0.0001, 10),
                 U = rep(0.01, 540))

# Based on optim
init_list = list(B = c(5.24510075177552, -0.215188297213872, 2.23839744913426, 
                       20.9084495320374, -1.26311776378368, 1.36166512931033, 
                       0, -0.128347714260874, 13.9865428397586),
                 U = rep(0.01, 540))

ans$sigma_U = 1
# b = replicate(100, system.time({fit = optimizing(mod, data = ans, init = init_list)}))
fit = optimizing(mod, data = ans, init = init_list, verbose = TRUE)

b = replicate(100, try({optimizing(mod, data = ans,
                             init = sapply(init_list, function(x)
                                 rnorm(length(x), x, abs(x/8)),
                                 simplify = FALSE))}))

b = replicate(100, try({optimizing(mod, data = ans)}))

ww = do.call(rbind, b["par",])[,1:9]
ww = cbind(ww, ll = unlist(b["value",]))
i = ww[,"ll"] > -7500
summary(ww[i,])
pairs(ww[i,])
ww[which.max(ww[,"ll"]),]

