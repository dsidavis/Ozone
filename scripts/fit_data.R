library(pso)

opt_fun = optim
# Fit the data from a single sheet
source("R/fit_funs.R")
source("R/model_functions.R")
source("R/readUCD.R")
f = "data/UCD_WCA2002.xls"

d = readExp(f)

d = lapply(d, deltaFEV1)

fit = optim(list(Dos = 3500, K = 0.020, A = -0.020), fit_FEV1b,
            d = d,
            method = "L-BFGS-B",
            lower = c(5,log(2)/1000, -0.15),
            upper = c(2500, log(2)/1, 0))


fit2 = psoptim(list(Dos = 1100, K = 0.040, A = -0.040, sigma = 0.25), fit_FEV1b,
             d = d, cost = "loglik")
           #  method = "L-BFGS-B",
           # lower = c(5,log(2)/1000, -0.15),
           # upper = c(2500, log(2)/1, 0))
# Confirmed - get similar answers regardless of cost function (good sanity check)

# Grid search over values of Dos
dos = exp(seq(1,9, length.out = 20))
pars = c(0.69, -0.0135)

ans = sapply(dos, function(x) fit_FEV1b(pars = c(x, pars), d = d))

pdf("dos_sse.pdf")
plot(ans ~ dos, type = "b",
     xlab = "'Dos' value",
     ylab = "Cost",
     main = "SSE cost")
text(y = pretty(ans)[2], x = 800,
     paste(c("K =", "A ="), pars, collapse = "\n"), adj = c(0,0))
dev.off()
# Sanity check - should be similar
ans2 = sapply(dos, function(x) fit_FEV1b(pars = c(x, pars, 0.5), d = d, cost = "loglik"))

pdf("dos_loglik.pdf")
plot(ans2 ~ dos, type = "b",
     xlab = "'Dos' value",
     ylab = "Cost",
     main = "Log-likelihood cost")
text(y = pretty(ans2)[2], x = 800,
     paste(c("K =", "A ="), pars, collapse = "\n"), adj = c(0,0))
dev.off()                                                                                                                                                                                                                                
# Similar responses - confirmed

i = sapply(d, function(x) unique(x$person))

fit2 = lapply(unique(i), function(j){
    tmp = d[i == j]
    optim(c(Dos = 500, K = 0.020, A = -0.020), fit_FEV1b,
          d = tmp, method = "L-BFGS-B",
          lower = c(5,log(2)/1000, -0.15),
          upper = c(2500, log(2)/1, 0))
})

summary(do.call(rbind, lapply(fit2, '[[', "par")))
