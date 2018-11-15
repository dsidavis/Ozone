library(pso)

opt_fun = optim
# Fit the data from a single sheet
source("R/fit_funs.R")
source("R/model_functions.R")
source("R/readUCD.R")
f = "data/UCD_WCA2002.xls"

d = readExp(f)

d = lapply(d, deltaFEV1)

fit = optim(list(Dos = 1, K = 0.040, A = -0.040), fit_FEV1b,
            d = d, method = "L-BFGS-B",
            lower = c(5,log(2)/1000, -0.15),
            upper = c(2500, log(2)/1, 0))

i = sapply(d, function(x) unique(x$person))

fit2 = lapply(unique(i), function(j){
    tmp = d[i == j]
    optim(c(Dos = 900, K = 0.040, A = -0.040), fit_FEV1b,
          d = tmp, #method = "L-BFGS-B",
          lower = c(5,log(2)/1000, -0.15),
          upper = c(2500, log(2)/1, 0))
})

summary(do.call(rbind, lapply(fit2, '[[', "par")))
