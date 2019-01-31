source("R/fit_funs.R")
source("R/model_functions.R")
source("R/readUCD.R")
f = "data/UCD_WCA2002.xls"

d = readExp(f)
d = lapply(d, deltaFEV1)


a = fit_FEV1b(c(1, log(2)/1000, -0.15), d = d)

source("R/deltaX2.R")
dyn.load("src/deltaX.so")
b = fit_FEV1b(c(1, log(2)/1000, -0.15), d = d)
