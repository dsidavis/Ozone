# Fit the data from a single sheet
source("R/fit_funs.R")
source("R/model_functions.R")

f = "data/UCD_WCA2002.xls"

d = readExp(f)

d = lapply(d, deltaFEV1)

fit = optim(list(Dos = 200, K = 0.040, A = -0.040), fit_FEV1b, d = d, method = "BFGS")
