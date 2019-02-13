# Munge the data from SAS to fit in the Williams model
# M. Espe Feb 2019

### Data needed ####
# int max_timepts; // max number of time points
# int max_n_dFEV1;
# int n_obs; // individual x study x exposure = obs
# int n_ind; // individuals
# int n_dFEV1[n_obs]; // n dFEV1 measurements per obs
# int n_timepts[n_obs]; // n timepts per obs
# int ind[n_obs]; 
# vector[n_ind] age;
# vector[n_ind] BMI;
# vector[n_obs] BSA;
# // These are padded with zeros
# vector[max_timepts] Ve[n_obs];
# vector[max_timepts] Vs[n_obs];
# vector[max_timepts] Cm[n_obs];
# vector[max_timepts] Cs[n_obs];
# vector[max_timepts] Time[n_obs];
# int dFEV1_measure_idx[n_obs, max_n_dFEV1];
# vector[max_n_dFEV1] dFEV1[n_obs];
# real<lower = 0> sigma_U; // prior sd of random effects

tars = list.files("~/Downloads", pattern = "*bdat$", full.names = TRUE)

sas = lapply(tars, haven::read_sas)


