functions{
  vector get_XB5(vector Cm, vector Cs, vector Vs, vector Ve,
				 vector BSA, int[] Time,
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
	  Ta = j > 1 ? Time[j-1] : 0;
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
  int n_ind;
  vector[n_ind] age;
  vector[n_ind] BMI;
  
	
}

parameters{

}
model{

}
generated quantities{

}
