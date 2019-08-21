# include <stdio.h>
# include <stdlib.h>
# include "lu.h"

double* lu_sparse(double *a, int *asub, int *xa, int n, int nzl, int nzu, int *perm_c, int *perm_r)
{
  double *L, *U;
  double sum_y = 0.0, sum_x = 0.0;
  double *y, *x;
  double *check_sum, sum = 0.0;
  double check_max; 
  int i, j, r, k, current_column;
  int *p, *row_num;
  double *l, *u, l_value;
  double temp = 0.0;
  int record_order_temp;
  FILE *fp_L = fopen("result/Sparse_L.txt", "w");
  FILE *fp_U = fopen("result/Sparse_U.txt", "w");
  FILE *result = fopen("result/Sparse_solution.txt", "w");
  FILE *pviot = fopen("result/Sparse_pviot.txt", "w");
  L = ( double * )malloc(sizeof(double *) * n);
  U = ( double * )malloc(sizeof(double *) * n);
  int sum_nonz_L = 0, sum_nonz_U = 0, sum_pviot_num = 0;
  double sum_element = n * n;
  double L_ratio, U_ratio;
  double start, finish, start_l, finish_l, time_l, start_checksum, finish_checksum, time_checksum, start_max, finish_max, 
  time_max, start_exchange, finish_exchange, time_exchange, start_u, finish_u, time_u, start_LU, finish_LU, time_LU;
  
  p = ( int * )malloc(sizeof(int) * n );
  row_num = ( int *)malloc(sizeof(int) * n );
  l = ( double *)malloc(sizeof(double) * n );
  u = ( double *)malloc(sizeof(double) * n );
  x = ( double *)malloc(sizeof(double) * n );
  check_sum = ( double *)malloc(sizeof(double) * n );
  for ( i = 0; i < n; i++ )
  {
    p[i] = xa[i];
    row_num[i] = xa[i+1] - xa[i];
//     printf(" l[%d] = %lf\n ", i, l[i]);
//     printf(" u[%d] = %lf\n ", i, u[i]);
  }
  //printf(" before the LU decomposition ");
  start = microtime();
  for ( r = 0; r < n; r++ )
  {
    start_l = microtime();
    // current_column存储的是要计算的当前列
    current_column = perm_c[r];
    for ( i = r; i < n; i++ )
    {
      if ( asub[p[i]] == r && row_num[i] > 0 )
      {
	  l[i] = a[p[i]];
	  p[i]++;
	  row_num[i]--;
      }
    }
    finish_l = microtime() - start_l;
    time_l += finish_l;
    
    start_checksum = microtime();
    for ( i = r; i < n; i++ )
    {
      for ( k = 0; k < r; k++ )
      {
	sum += L[i][k]*U[k][r];
      }
      check_sum[i] = l[i] - sum;
      l[i] = l[i] - sum;
      sum = 0.0;
    }         
    finish_checksum = microtime() - start_checksum;
    time_checksum += finish_checksum;
   
    check_max = check_sum[r];
    record_order_temp = r;
    start_max = microtime();
    for ( k = r+1; k < n; k++ )
    {
      if ( Abs(check_sum[k]) > Abs(check_max) )
      {
	check_max = check_sum[k];
	record_order_temp = k;
      }
    }
    finish_max = microtime() - start_max;
    time_max += finish_max;
    
    fprintf(pviot, "perm_r[%d] = %d\n", r, record_order_temp);
    start_exchange = microtime();
    U[r][r] = check_max;
    i = p[r];
    j = row_num[r];
    l_value = l[r];
    p[r] = p[record_order_temp];
    row_num[r] = row_num[record_order_temp];
    l[r] = l[record_order_temp];
    p[record_order_temp] = i;
    row_num[record_order_temp] = j;
    l[record_order_temp] = l_value;
    // exchange the row of r and record_order_temp 
    if ( record_order_temp != r )
    {
      sum_pviot_num ++;
      for ( i = 0; i < r; i++ )
      {
	temp = L[r][i];
	L[r][i] = L[record_order_temp][i];
	L[record_order_temp][i] = temp;
      }
    }
    finish_exchange = microtime() - start_exchange;
    time_exchange += finish_exchange;
    
    start_u = microtime();
    for ( j = 0; j < row_num[r]; j++ )
    {
      u[asub[p[r]+j]] = a[p[r]+j];
    }
    finish_u = microtime() - start_u;
    time_u += finish_u;
    
    start_LU = microtime();
    for ( i = r+1; i < n; i++ ) 
    {
      L[i][r] = l[i] / U[r][r];
      for ( k = 0; k < r; k++ )
      { 
	sum += L[r][k]*U[k][i];
      }
      U[r][i] = u[i] - sum;
      sum = 0.0;
    }
    finish_LU = microtime() - start_LU;
    time_LU += finish_LU;
    
    for ( i = 0; i < n; i++ )
    {
      l[i] = 0.0;
      u[i] = 0.0;
    }
  }
  finish = microtime() - start;
  printf("sum_pviot_num = %d\n", sum_pviot_num);
  printf("The time of LU decomposition is %lf\n", finish);
  printf("The time of generating l array is %lf\n", time_l);
  printf("The time of generating check_sum array is %lf\n", time_checksum);
  printf("The time of generating max value is %lf\n", time_max);
  printf("The time of exchanging the value is %lf\n", time_exchange);
  printf("The time of generating u array is %lf\n", time_u);
  printf("The time of generating LU value is %lf\n", time_LU);
  printf("The sum short time if %lf\n", time_l+time_max+time_exchange+time_u);
  y = ( double * )malloc(sizeof(double) * n);
  for ( i = 0; i < n; i++ )
  {
    for( j = 0; j < n; j++ )
    {
      if ( L[i][j] != 0 )
      {
	fprintf(fp_L, "L[%d][%d]=%lf\n", i, j, L[i][j]);
      }
    }
  }
  for ( i = 0; i < n; i++ )
  {
    for( j = 0; j < n; j++ )
    {
      //printf("U[%d][%d]=%lf\n", i, j, U[i][j]);
      if ( U[i][j] != 0 )
      {
	fprintf(fp_U, "U[%d][%d]=%f\n", i, j, U[i][j]);
      }
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
