#include "lib/lu.h"
#include <stdio.h>
#include <stdlib.h>


bool equal( double a, double b )
{
  if ( Abs(a-b) < 0.0001)
  {
    return true;
  }
  else
  {
    return false;
  }
}
int main( int argc, char *argv[] )
{
    FILE *fp;
    double *GP_x, *GP_x_v2, *GP_x_v3, *x, *test, *GP_x_amd; 
    double **A;
    double start_lu, finish_lu, start_GP, start_GP_v2, start_dense, finish_GP, finish_GP_v2, start_GP_v3, finish_GP_v3, start_GP_amd, 
    finish_GP_amd;
    double time_GP, time_GP_v2, pivot_ratio;
    double   *a;
    int      *asub, *xa;
    int m, n, nnz;
    int i, j;
    fp = fopen(argv[1], "r");
    dreadMM(fp, &m, &n, &nnz, &a, &asub, &xa);
//     printf("-----------------------super_LU----------------------\n");
//     start_lu = microtime();
//     x = ( double *)malloc(sizeof(double) * n );
//     x = super_lu( argv[1], x );
//     finish_lu = microtime() - start_lu;
//     for ( i = 0; i < n; i++ )
//     {
//       printf(" x[%d] = %lf\n", i, x[i]);
//     }
//     test = x + n;
//     for ( i = 0; i < 3*n; i++ )
//     {
//       printf(" test[%d] = %lf\n ", i, test[i]);
//     }
//     printf("The average cost time of super_LU function is: %f seconds\n", finish_lu);
    printf("-----------------------LU_GP----------------------\n");
    start_GP = microtime();
    GP_x = lu_gp( a, asub, xa, n );
    finish_GP = microtime() - start_GP;
    pivot_ratio = ( GP_x[n] / finish_GP ) * 100;
    printf(" the pivot time ration is: %lf\n ", pivot_ratio);
    printf("The average cost time of LU_GP function is: %f seconds\n", finish_GP);
//     printf("-----------------------LU_GP_v2----------------------\n");
//     start_GP_v2 = microtime();
//     GP_x_v2 = lu_gp_v2( a, asub, xa, n );
//     finish_GP_v2 = microtime() - start_GP_v2;
//     printf("The average cost time of LU_GP_v2 function is: %f seconds\n", finish_GP_v2);
//     printf("-----------------------LU_GP_v3----------------------\n");
//     start_GP_v3 = microtime();
//     GP_x_v3 = lu_gp_v3( a, asub, xa, n );
//     finish_GP_v3 = microtime() - start_GP_v3;
//     for ( i = 0; i < n; i++ )
//     {
//       printf(" x[%d] = %lf\n", i, GP_x_v3[i]);
//     }
//     printf("The average cost time of LU_GP_v23function is: %f seconds\n", finish_GP_v3);
    printf("-----------------------LU_SPARSE----------------------\n");
    start_GP_amd = microtime();
    GP_x_amd = lu_sparse( a, asub, xa, n );
    finish_GP_amd = microtime() - start_GP_amd;
    //pivot_ratio = ( GP_x[n] / finish_GP ) * 100;
    //printf(" the pivot time ration is: %lf\n ", pivot_ratio);
    printf("The average cost time of LU_SPARSE function is: %f seconds\n", finish_GP_amd);
    printf("Normal end of execution");
    return 0;
}
