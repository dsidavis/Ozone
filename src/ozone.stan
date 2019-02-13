/*
  Model for change in %FEV1 after exposure to O3

  Matt Espe
  Aug 2018
*/

functions{
  
  vector seq(int start, int end){
	int N = (1 + end) - start;
	vector[N] x;
	for(i in 1:N)
	  x[i] = start + i - 1;
	return x;	
  }
  
  vector UOS(real DR, vector CumFrDos, real Dos, vector t) {
	int n = dims(t)[1];
	vector[n] UOS;
	for(i in 1:n)
	  UOS[i] = DR / (1 + exp(-20 * (t[i] - (t[i] / CumFrDos[i]))));
	
	return UOS;
  }
  
  vector deltaX(vector uos, real k, real x_prev){
	int n = dims(uos)[1];
	real r = (1 - exp(-k));
	vector[n] ceoss = uos / k;
	vector[n] x = rep_vector(0.0, n);

	x[1] = x_prev;

	if(n > 1)
	  for(i in 2:n)
		x[i] = x[i-1] + ((ceoss[i-1] - x[i-1]) * r);
	
	return x;
  }

  vector experimentFEV1(vector o3, vector ve, int[] t_stop,
						real dos, real k, real a){
	int m = dims(o3)[1];
	int n = max(t_stop);
	vector[n] dFEV1 = rep_vector(0, n);
	real x_previous = 0;
	real FrDos;
	real FrDos_previous = 0;
	real uos = 0;
	real dr = 0;
	
	real cd = 0;

	for(i in 1:(n - 1)){
	  int n_t = t_stop[i+1] - t_stop[i]+1;
	  vector[n_t] t = seq(t_stop[i]+1, t_stop[i+1]);
	  vector[n_t] CumFrDos;
	  CumFrDos[1] = FrDos_previous;
	  
	  dr = o3[i] * ve[i] * 1.96;
	  FrDos = dr / dos;
	  for(j in 1:n_t)
		CumFrDos[j] = FrDos_previous + FrDos * j;
	  
		
		
	}
  }
}

data{
  int n_blocks;
  int n_fev1;
  vector[n_blocks] o3;
  vector[n_blocks] ve;
  int t_stop[n_blocks];
  int fev1_pts[n_fev1];
  vector[n_fev1] fev1;
}

parameters{
  real dos;
  real k;
  real a;
  real<lower=0> sigma;
}

model{
  vector[max(t_stop)] pred_fev1;

  dos ~ normal(1000, 1000);
  k ~ normal(0,1);
  a ~ normal(0,1);

  pred_fev1 = experimentFEV1(o3, ve, t_stop, dos, k, a);
  target += normal_lpdf(fev1 | pred_fev1[fev1_pts], sigma);
}
generated quantities{

}
