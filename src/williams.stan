functions{
  vector get_XB5(vector Cm, vector Cs, vector Vs, vector Ve,
				 vector BSA, vector Time,
				 real B5, real B6){
	int N = size(Cm);
	vector[N] XB5 = rep_vector(0.0, N);
	real Vm;
	real Ta;
	real Tb;
	real TD;
	real XB5_previous;
	
	for(j in 1:N){
	  Vm = (Ve[j] ^ B6) / (BSA[j] ^ B6);
	  Ta = j > 1 ? Time[j-1] : 0.0;
	  Tb = Time[j];
	  TD = (Tb - Ta);
	  XB5_previous = j > 1 ? XB5[j-1] : 0;
	  
	  XB5[j] = XB5_previous * (exp(-B5 * TD)) +
		(Cm * Vm * (B5^-1)) * (1 - exp(-B5 * TD)) +
		(Cm * Vs * (B5^-2)) * (((1 - B5 * Ta) * exp(-B5 * TD)) - (1-B5*Tb)) +
		(Cs * Vm * (B5^-2)) * (((1 - B5 * Ta) * exp(-B5 * TD)) - (1-B5*Tb)) +
		(Cs * Vs * (B5^-3)) * (((-2 + (2 * B5 * Ta) -
								 (B5^2 * Ta^2)) * exp(-B5 * TD)) -
							   (-2 + (2 * B5 * Tb) - (B5^2 * Tb^2)));
    }
	return XB5;
  }

  vector get_pop_median (vector XB5, real age_c, real BMI_c, 
						 real B1, real B2, real B3,
						 real B4, real B8, real B9) {
	int N = size(XB5);
	vector[N] XB5G;
	real F1 = B1 + B2 * age_c + B8 * BMI_c;
	vector[N] T1;
	real T2 = 1 + B4;
	vector[N] Median;
	
	for(n in 1:N)
	  XB5G[n] = XB5[n] <= B9 ? XB5[n] - B9 : 0;

	T1 = 1 + B4 * exp(-B3 .*XB5G);
	Median = F1 * (1/T1 - 1/T2);
	
	return Median;
  }
}

data{
  int max_timepts; // max number of time points
  int max_n_dFEV1;
  int n_obs; // individual x study x exposure = obs
  int n_ind; // individuals
  int n_dFEV1[n_obs]; // n dFEV1 measurements per obs
  int n_timepts[n_obs]; // n timepts per obs
  int ind[n_obs]; 
  vector[n_ind] age;
  vector[n_ind] BMI;
  vector[n_obs] BSA;
  // These are padded with zeros
  vector[max_timepts] Ve[n_obs];
  vector[max_timepts] Vs[n_obs];
  vector[max_timepts] Cm[n_obs];
  vector[max_timepts] Cs[n_obs];
  vector[max_timepts] BSA[n_obs];
  vector[max_timepts] Time[n_obs];
  vector[max_n_dFEV1] dFEV1_measure_idx[n_obs];
  vector[max_n_dFEV1] dFEV1[n_obs];
}

transformed data{
  vector[n_ind] age_c = age - 23.8;
  vector[n_ind] BMI_c = BMI - 23.1;
}

parameters{
  real B1;
  real B2;
  real B3;
  real B4;
  real B5;
  real B6;
  real B8;
  real B9;
  vector[n_ind] U; // random effects by ind.
}

model{
  
  for(n in 1:n_obs){
	int idx = n_timepts[n];
	vector[idx] XB5 = get_XB5(Cm[n][:idx], Cs[n][:idx],
									   Vs[n][:idx], Ve[n][:idx],
									   BSA[ind[n]], Time[n][:idx],
									   B5, B6);
	vector[idx] med = get_pop_median(XB5, age_c[ind[n]], BMI_c[ind[n]],
											  B1, B2, B3, B4, B8, B9);
	

  }
}
generated quantities{

}
