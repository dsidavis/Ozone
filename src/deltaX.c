#include <Rdefines.h>


/* 
  fev_base  numeric(1) - scalar real
  ceoss  numeric(n_t)
  r_r   scalar containing named value for K

  r_K - scalar
 */
SEXP
R_deltaXLoop(SEXP r_ceoss,  SEXP fev_base, SEXP r_K)
{
    int i, n_t = Rf_length(r_ceoss);
    SEXP r_x;
    PROTECT(r_x = NEW_NUMERIC(n_t));
    double *x = REAL(r_x);
    double *ceoss = REAL(r_ceoss);
    x[0] = REAL(fev_base)[0];
//    double r = REAL(r_r)[0];
    double r = 1.0 - exp( - REAL(r_K)[0]);

/*
        for(i in 2:n_t)
            # sum as we go
            x[i] = x[i-1] + ((ceoss[i - 1] - x[i - 1]) * r)    
*/
    for(i = 1; i < n_t; i++) {
	x[i] = x[i-1] + ((ceoss[i - 1] - x[i - 1]) * r);
    }
    UNPROTECT(1);
    return(r_x);
}
