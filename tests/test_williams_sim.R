source("R/sim_williams.R")
source("scripts/munge_data_williams.R")

n_sim = 20
sim_idx = seq(n_sim)
XB5s = lapply(sim_idx, function(i, B5, B6){
    
    with(ans, get_XB5(Cm[i,], Cs[i,], Vs = 0, Ve[i,],
                            BSA[i], Time[i,],
                            B5, B6))
}, B5 = 0.003224, B6 = 0.886759)

 m = lapply(sim_idx, function(i, age_c, BMI_c, B1, B2, B3, B4, B8, B9){
    get_pop_median(XB5s[[i]], ans$age[i] - 23.8, ans$BMI[i] - 23.1, 
                          B1, B2, B3,
                          B4, B8, B9)

}, B1 = 11.091059, B2 = 0.287275, B3 = 0.014862,
B4 = 13.449713, B8 = 0.546653, B9 = 59.949843)

