source("R/fit_funs.R")
source("tests/test_model_funs.R")

dd = lapply(1:100, function(i) {
    ans = sim_experiment()
    ans$subject = i
    ans
})

dd = do.call(rbind, dd)

# fit_FEV1(df = ind_data)
fit_FEV1(df = d)
optim(list(Dos = 1000, K = .0400, A = -.0400), fit_FEV1, df = dd)

optim(list(Dos = 200, K = 10.040, A = -.040), f, d = d3, method = "BFGS")
