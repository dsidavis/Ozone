# which have a >0 dFEV1
source("R/model_functions.R")
source("R/fit_funs.R")

i = sapply(sas[[3]], function(x) any(x$dFEV1 > 0, na.rm = TRUE))
sas[[1]][i]
fit_FEV1b(c(1100, 0.69, -.02), sas[[1]])

# Try the optimization
fit = optim(par = c(Dos = 900, K = 0.65, A = -0.020), fn = fit_FEV1b,
            d = sas[[3]],
            lower = c(5,log(2)/1000, -0.15),
            upper = c(2500, log(2)/1, 0),
            control = list(parscale = c(1000, 0.1, 0.01),
                           trace = 6))
