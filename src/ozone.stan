/*
  Model for change in %FEV1 after exposure to O3

  Matt Espe
  Aug 2018
*/

functions{

  vector UOS(real DR, vector CumFrDos, real Dos, vector t) {
	int n = dims(t)[1];
	vector[n] UOS = DR / (1 + exp(-20 * (t - (t / CumFrDos))));
	
	return UOS;
  }
  
  vector deltaX(vector uos, real k, real x_prev){
	int n = dims(uos)[1];
	real r = (1 - exp(-k));
	vector coess = uos / k;
	vector x = rep_vector(0.0, n);

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
	real uos = 0;
	real dr = 0;
	real cd = 0;

	if(t_stop[1] != 0)
	  t_stop = append_array(0.0, t_stop);

	for(i in 1:(n - 1){
		int n_t = t_stop[i+1] - t_stop[i]+1;
		

		
	  }
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
