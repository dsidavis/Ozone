/*
  Model for change in %FEV1 after exposure to O3

  Matt Espe
  Aug 2018
*/

functions{

  real UOS(real DR, real CD, real Dos, real t) {
	return (DR / (1 + exp(-20 * (t - (Dos * (t / CD))))));
	  }
  
  real deltaX(real uos, real k, real x_prev){
	return (x_prev + ((uos/k) - x_prev) * (1 - exp(-k)));
	  }

  vector experimentFEV1(vector o3, vector ve, int[] t_stop,
						real dos, real k, real a){
	int m = dims(o3)[1];
	int n = max(t_stop);
	vector[n] ans = rep_vector(0, n);
	real x_last = 0;
	int t = 1;
	int cur_block = 1;
	real uos = 0;
	real dr = 0;
	real cd = 0;
	
	while(cur_block <= m){
	  while(t <= t_stop[cur_block]){
		dr = o3[cur_block] * ve[cur_block] * 1.96;
		cd += dr;
		uos = UOS(dr, cd, dos, t);
		x_last = deltaX(uos, k, x_last);
		ans[t] = x_last;
		t += 1;
	  }
	  cur_block += 1;
	}
	
	return(ans * a);
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
