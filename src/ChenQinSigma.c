#include<R.h>
#include<Rmath.h>

void tracehat(	double *X,
		double *Xbar,
		int *J,
		int *K,
		int *n,
		int *p,
		double *result,
		double *xbar_jk)
{

	unsigned int l,m,lj,lk;
	double /* xbar_jk[*p], */ accum = 0, accum_a = 0, accum_b = 0, accum_c = 0;

	for( l = 0 ; l < (*n) * (*n - 1 ) / 2 ; l++ )
	{
		for( m = 0 ; m < *p ; m++ )
		{
			xbar_jk[m] = ( *n * Xbar[m] - X[J[l] * (*p) + m] - X[K[l] * (*p) + m] ) / ( *n - 2 );
		}

		accum_a = 0;
		accum_b = 0;
		accum_c = 0;		

		for( m = 0 ; m < *p ; m++ )
		{
			lj = J[l] * (*p) + m;
			lk = K[l] * (*p) + m;

			accum_a += X[lj] * X[lk];
			accum_b += X[lj] * xbar_jk[m];
			accum_c += X[lk] * xbar_jk[m];
		}
	
		accum += (accum_a - accum_b) * (accum_a - accum_c ) /  (*n * (*n - 1)) ;		
	} 
	
	*result = 2 * accum;
}

void trace12hat()
{

}
