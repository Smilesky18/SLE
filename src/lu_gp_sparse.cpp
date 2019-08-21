# include <stdio.h>
# include <stdlib.h>
# include "lu.h"

bool equal( double a, double b )
{
  if ( Abs(a-b) < 0.001)
  {
    return true;
  }
  else
  {
    return false;
  }
}

double Abs(double x)
{
  return x < 0 ? -x : x;
}

void* lu_gp_sparse(double *a, int *asub, int *xa, int n, int nzl, int nzu, int *perm_c, int *perm_r, int *asub_L, int *xa_L, int *asub_U, int *xa_U, double *l_data, double *u_data)
{
  int sum_l = 0, sum_u = 0;
  double *L, *U, *xx;
  double U_diag;
  int i, j, r, k, current_column;
  int *row_num;
  double *l, *u, l_value;
  double temp = 0.0;
  L = ( double * )malloc(sizeof(double *) * nzl);
  U = ( double * )malloc(sizeof(double *) * nzu);
  int sum_nonz_L = 0, sum_nonz_U = 0, sum_pviot_num = 0;
  double start, finish, start_xx, finish_xx, sum_pivot = 0.0, start_cal, finish_cal, sum_cal = 0.0, start_L, start_U, end_L, end_U, sum_xx = 0.0, sum_L = 0.0, sum_U = 0.0, start_compare, end_compare, start_xx_0, end_xx_0, sum_xx_0 = 0.0;
  int m = n + 1;
  xx = ( double *)malloc(sizeof(double) * n );
  row_num = ( int *)malloc(sizeof(int) * n );
  
  for ( i = 0; i < n; i++ )
  {
    row_num[i] = i;
    xx[i] = 0;
    L[xa_L[i]] = 0.0;
  }
//   L[0] = 0;
  // column-oriented G/P algorithm without partial pivoting
  start = microtime();
  for ( k = 0; k < n; k++ )
  {
    current_column = perm_c[k];
    start_xx = microtime();
    for ( j = xa[current_column]; j < xa[current_column+1]; j++ )
    {
      xx[perm_r[asub[j]]] = a[j];
    }
    finish_xx = microtime() - start_xx;
    sum_xx += finish_xx;
    start_cal = microtime();
    for ( j = 0; j < k; j++ )
    {
        if( xx[j] != 0 )
        {
            for ( i = xa_L[j]+1; i < xa_L[j+1]; i++ )
            {
                xx[asub_L[i]] = xx[asub_L[i]] - xx[j]*L[i];
            }
        }
    }
    finish_cal = microtime() - start_cal;
    sum_cal += finish_cal;
//     printf("every_time_cal = %lf\n", finish_cal);
    U[xa_U[k]] = xx[k];
    start_U = microtime();
    for ( i = xa_U[k]; i < xa_U[k+1]; i++ )
    {
      U[i] = xx[asub_U[i]];
    }
    end_U = microtime() - start_U;
    sum_U += end_U;
    U_diag = U[i-1];
    start_L = microtime();
    for ( i = xa_L[k]; i < xa_L[k+1]; i++ )
    {
      L[i] = xx[asub_L[i]] / U_diag;
    }
    end_L = microtime() - start_L;
    sum_L += end_L;
    start_xx_0 = microtime();
    memset(xx, 0, sizeof(double)*n);
//     for ( i = xa_L[k]; i < xa_L[k+1]; i++ )
//     {
//         xx[asub_L[i]] = 0.0;
//     }
//     for ( i = xa_U[k]; i < xa_U[k+1]; i++ )
//     {
//         xx[asub_U[i]] = 0.0;
//     }
//     for ( i = 0; i < n; i++ )
//     {
//       xx[i] = 0;
//     }
    end_xx_0 = microtime() - start_xx_0;
    sum_xx_0 += end_xx_0;
  }
  finish = microtime() - start;
    printf("sum_cal = %lf\n", sum_cal);
    printf("sum_xx = %lf\n", sum_xx);
    printf("sum_L = %lf\n", sum_L);
    printf("sum_U = %lf\n", sum_U);
    printf("sum_xx_0 = %lf\n", sum_xx_0);
    printf("sum_gp = %lf\n", finish);
  start_compare = microtime();
    for ( i = 0; i < nzl; i++ )
    {
        if ( equal(L[i], l_data[i]) ) continue;
        else sum_l++;
    }
    for ( i = 0; i < nzu; i++ )
    {
        if ( equal(u_data[i], U[i]) ) continue;
        else 
        {
            sum_u++;
            printf(" %lf  %lf \n", u_data[i], U[i]);
        }
    }
    end_compare = microtime() - start_compare;
    printf("sum_l = %d\n", sum_l);
    printf("sum_u = %d\n", sum_u);
    printf("the time of comparing is: %lf\n", end_compare);
}

