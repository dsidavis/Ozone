
void
timeLoop(double *x, int *len, double *K, double *UOS, double *r)
{
    for(int i = 1; i < *len; i++)
	x[i] = x[i-1] + (UOS[i-1]/(*K) - x[i-1]) * r[0];
}
