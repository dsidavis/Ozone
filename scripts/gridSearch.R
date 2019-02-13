# Fit the data from a single sheet
source("R/fit_funs.R")
source("R/model_functions.R")
source("R/readUCD.R")
f = "data/UCD_WCA2002.xls"

d = readExp(f)
d = lapply(d, deltaFEV1)


numIntervals = 50
dos = seq(5, 2500, length.out = numIntervals)
#dos = exp(seq(1,9, length.out = numIntervals))
K = seq(log(2)/1000, log(2)/1, length.out = numIntervals)
A = seq(-0.15, 0, length.out = numIntervals)
#A = seq(0, 1, length.out = numIntervals) 
#pars = c(0.69, -0.0135)
pars = expand.grid(DOS = dos, K = K, A = A)

gr = apply(pars, 1, function(x) fit_FEV1b(unlist(x), d = d))

library(parallel)
cl = makeCluster(60, "FORK")
grp = parallel::parApply(cl, pars, 1, function(x) fit_FEV1b(pars = unlist(x), d = d))


