source("R/fit_funs.R")
source("tests/test_model_funs.R")

ind_data$FEV1 = aa[-13]
names(ind_data)[3] = "t_end"
ind_data$t_end = ind_data$t_end[-1]
ind_data = as.data.frame(ind_data)
ind_data$subject = 1
dd = lapply(1:100, function(i) {
    x = ind_data
    x$FEV1 = rnorm(nrow(x), mean = x$FEV1, sd = 0.1)
    x$subject = i
    x})

dd = do.call(rbind, dd)

fit_FEV1(df = ind_data)

optim(list(Dos = 1000, K = 0.15, A = -0.15), fit_FEV1, df = dd, method = "L-BFGS-B")
