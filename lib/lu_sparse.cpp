# include <stdio.h>
# include <stdlib.h>
# include "lu.h"

void lu_sparse(double *a, int *asub, int *xa, double *x, int n)
{
  double **L, **U;
  double sum_y = 0.0, sum_x = 0.0;
  double *y;
  double *check_sum, sum = 0.0;
  double check_max; 
  int i, j, r, k;
  int *p, *row_num;
  double *l, *u, l_value;
  double temp = 0.0;
  int record_order_temp;
  FILE *fp_L = fopen("result/Sparse_L.txt", "w");
  FILE *fp_U = fopen("result/Sparse_U.txt", "w");
  FILE *result = fopen("result/Sparse_solution.txt", "w");
  L = ( double ** )malloc(sizeof(double *) * n);
  U = ( double ** )malloc(sizeof(double *) * n);
  int sum_nonz_L = 0, sum_nonz_U = 0;
  double sum_element = n * n;
  double L_ratio, U_ratio;
  for ( i = 0; i < n; i++ )
  {
    L[i] = ( double * )malloc(sizeof(double) * n);
    U[i] = ( double * )malloc(sizeof(double) * n);
  }
  for ( i = 0; i < n; i++ )
  {
    L[i][i] = 1.0;
  }
  p = ( int * )malloc(sizeof(int) * n );
  row_num = ( int *)malloc(sizeof(int) * n );
  l = ( double *)malloc(sizeof(double) * n );
  u = ( double *)malloc(sizeof(double) * n );
  check_sum = ( double *)malloc(sizeof(double) * n );
  for ( i = 0; i < n; i++ )
  {
    p[i] = xa[i];
    row_num[i] = xa[i+1] - xa[i];
  }
  //printf(" before the LU decomposition ");
  for ( r = 0; r < n; r++ )
  {
    for ( i = r; i < n; i++ )
    {
      if ( asub[p[i]] == r && row_num[i] > 0 )
      {
	  l[i] = a[p[i]];
	  p[i]++;
	  row_num[i]--;
      }
    }
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
    check_max = check_sum[r];
    record_order_temp = r;
    for ( k = r+1; k < n; k++ )
    {
      if ( Abs(check_sum[k]) > Abs(check_max) )
      {
	check_max = check_sum[k];
	record_order_temp = k;
      }
    }
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
      for ( i = 0; i < r; i++ )
      {
	temp = L[r][i];
	L[r][i] = L[record_order_temp][i];
	L[record_order_temp][i] = temp;
      }
    }
    for ( j = 0; j < row_num[r]; j++ )
    {
      u[asub[p[r]+j]] = a[p[r]+j];
    }
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
    for ( i = 0; i < n; i++ )
    {
      l[i] = 0.0;
      u[i] = 0.0;
    }
  }
  y = ( double * )malloc(sizeof(double) * n);
  for ( i = 0; i < n; i++ )
  {
    for( j = 0; j < n; j++ )
    {
      //printf("L[%d][%d]=%lf\n", i, j, L[i][j]);
      fprintf(fp_L, "L[%d][%d]=%lf\n", i, j, L[i][j]);
      if ( L[i][j] != 0 )
	sum_nonz_L++;
    }
  }
  L_ratio = (1 - sum_nonz_L/sum_element)*100;
  //printf("sum_element = %f\n", sum_element);
  //printf("sum_nonz_L = %d\n", sum_nonz_L);
  //printf("The sparsity of L is: %lf%%\n", L_ratio);
  //printf("The element of U: ");
  //printf(" the value of sum_nonz_U is: %d\n", sum_nonz_U);
  for ( i = 0; i < n; i++ )
  {
    for( j = 0; j < n; j++ )
    {
      //printf("U[%d][%d]=%lf\n", i, j, U[i][j]);
      fprintf(fp_U, "U[%d][%d]=%f\n", i, j, U[i][j]);
      if ( U[i][j] != 0 )
	sum_nonz_U++;
    }
  }
  U_ratio = (1 - sum_nonz_U/sum_element)*100;
  //printf("sum_nonz_U = %d\n", sum_nonz_U);
  //printf("The sparsity of U is: %lf%%\n", U_ratio);
  //solve Ly=b && Ux=y
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
}

