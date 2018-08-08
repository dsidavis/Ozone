source("R/model_functions.R")

# Using typical values from UCD_WCA2000b.xls and Schelegle et al.
O3 = 0.1
VE = 30
t = 1:50
dos = 1100
uos = UOS(O3, VE, t, Dos = dos)

stopifnot(length(uos) == length(t))
stopifnot(identical(uos[length(uos)], 5.88))

k = 0.02
a = -0.02
fev0 = 0
fev = FEV1(deltaX(uos, fev_base = fev0, K = k), A = a)

stopifnot(length(uos) == length(fev))
stopifnot(identical(fev[1], 0))

# simplified version of data.
ind_data = list(O3 = rep(rep(c(0.123, 0), each = 2), 3),
                Ve = rep(c(30, 13), length.out = 12), 
                t_stop = c(0, 50, 63, 113, 126, 176, 200, 260, 273, 323, 336, 386, 399))

aa = experimentFEV1(ind_data$O3, ind_data$Ve, ind_data$t_stop, Dos = dos, K = k, A = a)
