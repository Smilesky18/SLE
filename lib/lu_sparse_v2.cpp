# include <stdio.h>
# include <stdlib.h>
# include "lu.h"

double* lu_sparse_column(double *a, int *asub, int *xa, int n)
{
  double **L, **U, *xx;
  double sum_y = 0.0, sum_x = 0.0;
  double *y, *x;
  double *check_sum, sum = 0.0;
  double check_max; 
  int i, j, r, k;
  int *p, *row_num;
  double *l, *u, l_value;
  double temp = 0.0;
  int record_order_temp;
  FILE *fp_L = fopen("result/Sparse_GP_L.txt", "w");
  FILE *fp_U = fopen("result/Sparse_GP_U.txt", "w");
  FILE *result = fopen("result/Sparse_GP_solution.txt", "w");
  FILE *pviot = fopen("result/Sparse_GP_pviot.txt", "w");
  L = ( double ** )malloc(sizeof(double *) * n);
  U = ( double ** )malloc(sizeof(double *) * n);
  int sum_nonz_L = 0, sum_nonz_U = 0, sum_pviot_num = 0;
  double sum_element = n * n;
  double L_ratio, U_ratio;
  double start, finish, start_l, finish_l, time_l, start_checksum, finish_checksum, time_checksum, start_max, finish_max, 
  time_max, start_exchange, finish_exchange, time_exchange, start_u, finish_u, time_u, start_LU, finish_LU, time_LU;
  for ( i = 0; i < n; i++ )
  {
    L[i] = ( double * )malloc(sizeof(double) * n);
    U[i] = ( double * )malloc(sizeof(double) * n);
  }
  for ( i = 0; i < n; i++ )
  {
    L[i][i] = 1.0;
  }
  xx = ( double *)malloc(sizeof(double) * n );
  x = ( double *)malloc(sizeof(double) * n );
  row_num = ( int *)malloc(sizeof(int) * n );
  for ( i = 0; i < n; i++ )
  {
    row_num[i] = i;
  }
  // column-oriented G/P algorithm with partial pivoting
  start = microtime();
  for ( k = 0; k < n; k++ )
  {
    record_order_temp = k;
    // xx = A[;k]
    for ( j = xa[k]; j < xa[k+1]; j++ )
    {
      xx[row_num[asub[j]]] = a[j];
    }
    start_checksum = microtime();
    for ( j = 0; j < k; j++ )
    {
      for ( i = j+1; i < n; i++ )
      {
	xx[i] -= xx[j]*L[i][j];
      }
    }
    finish_checksum = microtime() - start_checksum;
    time_checksum += finish_checksum;
    // find the max value in u[k:n,k]
    temp = xx[k];
    for ( i = k+1; i < n; i++ )
    {
      if ( Abs(xx[i]) > Abs(temp) )
      {
	temp = xx[i];
	record_order_temp = i;
      }
    }
    // change the k'th and max value
    if ( record_order_temp != k )
    {
      //row_num[k] = record_order_temp;
      if ( row_num[k] != k )
      {
	row_num[record_order_temp] = k;
      }
      else
      {
	row_num[record_order_temp] = k;
        row_num[k] = record_order_temp;
      } 
      // change the k'th and max value in xx
      temp = xx[k];
      xx[k] = xx[record_order_temp];
      xx[record_order_temp] = temp;
      // change the k'th and max value in L
      for ( i = 0; i < k; i++ )
      {
	temp = L[k][i];
	L[k][i] = L[record_order_temp][i];
	L[record_order_temp][i] = temp;
      }
    }
    U[k][k] = xx[k];
    for ( i = 0; i < k; i++ )
    {
      U[i][k] = xx[i];
    }
    for ( i = k+1; i < n; i++ )
    {
      L[i][k] = xx[i] / xx[k];
    }
    for ( i = 0; i < n; i++ )
    {
      xx[i] = 0;
    }
  }
  finish = microtime() - start;
  //printf(" before the LU decomposition ");
  printf("The time of LU_column decomposition is %lf\n", finish);
  printf("The time of check_sum is %lf\n", time_checksum);
  y = ( double * )malloc(sizeof(double) * n);
  for ( i = 0; i < n; i++ )
  {
    for( j = 0; j < n; j++ )
    {
      //printf("L[%d][%d]=%lf\n", i, j, L[i][j]);
      fprintf(fp_L, "L[%d][%d]=%lf\n", i, j, L[i][j]);
    }
  }
  for ( i = 0; i < n; i++ )
  {
    for( j = 0; j < n; j++ )
    {
      //printf("U[%d][%d]=%lf\n", i, j, U[i][j]);
      fprintf(fp_U, "U[%d][%d]=%f\n", i, j, U[i][j]);
    }
  }
  y[0] = 1.0;  
  for ( i = 1; i < n; i++ )
  {
    for ( k = 0; k < i; k++ )
    {
      sum_y += L[i][k]*y[k]; 
    }
    y[i] = 1.0 - sum_y;
    sum_y = 0.0;
  }
  x[n-1] = y[n-1] / U[n-1][n-1];
  for ( i = n-2; i >= 0; i-- )
  {
    for ( k = i + 1; k < n; k++ )
    {
      sum_x += U[i][k]*x[k];
    }
    x[i] = ( y[i] - sum_x ) / U[i][i];
    sum_x = 0.0;
  }
  for ( i = 0; i < n; i++ )
  {
    //printf("x[%d] = %f\n", i, x[i]);
    fprintf(result, "x[%d]=%lf\n", i, x[i]);
  }
  for ( i = 0; i < n; i++ )
  {
    //free(A[i]);
    free(L[i]);
    free(U[i]);
    //free(col_list[i]);
  }
  //free(x);
  free(y);
  free(check_sum);
  //free(change_order);
  return x;
}

