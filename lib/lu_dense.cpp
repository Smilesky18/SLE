#include <stdio.h>
#include <stdlib.h>
#include "lu.h"

double Abs(double x)
{
  return x < 0 ? -x : x;
}
void lu_dense(double **A, double *x, int n)
{
  double **L, **U;
  double sum_U = 0.0, sum_y = 0.0, sum_x = 0.0;
  double *y;
  double *check_sum, sum = 0.0, check_max;
  int i, j, r, k;
  int *record_order;
  double temp = 0.0;
  int record_order_temp;
  FILE *fp_L = fopen("result/Dense_L.txt", "w");
  FILE *fp_U = fopen("result/Dense_U.txt", "w");
  FILE *result = fopen("result/Dense_solution.txt", "w");
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
      for ( j = 0; j < n; j++ )
      {
	L[i][j] = 0.0;
	U[i][j] = 0.0;
      }
    }
  y = ( double * )malloc(sizeof(double) * n);
  check_sum = ( double * )malloc(sizeof(double) * n);
  record_order = ( int * )malloc(sizeof(int) * n);
  for ( i = 0; i < n; i++ )
  {
    record_order[i] = i;
  }
  temp = A[0][0];
  for ( i = 1; i < n; i++ )
  {
    if ( Abs(temp) < Abs(A[i][0]) )
    {
      record_order[0] = i;
      temp = A[i][0];
    }
  }
  if ( record_order[0] != 0 )
  {
    record_order_temp = record_order[0];
    for ( i = 0; i < n; i++ )
    {
      temp = A[0][i];
      A[0][i] = A[record_order_temp][i];
      A[record_order_temp][i] = temp;
    }
    temp = 0.0;
  }
  //record_order[0] = record_order_temp;
  //printf( " A[1][0] = %f\n ", A[1][0]);
  //the first row of U
  for ( i = 0; i < n; i++ )
  {
    U[0][i] = A[0][i];
  }
  //the first column of L
  for ( i = 1; i < n; i++ )
  {
    L[i][0] = A[i][0]/U[0][0];
  }
  //printf( " L[1][0] = %f\n", L[1][0] );
  //The diagonal element of L is 1
  for ( i = 0; i < n; i++ )
  {
    L[i][i] = 1;
  }
  // start the LU from the second row/column
  for ( r = 1; r < n; r++ )
  {
    for ( i = r; i < n; i++ )
    {
      for ( k = 0; k < r; k++ )
      {
	sum += L[i][k]*U[k][r];
      }
      check_sum[i] = A[i][r] - sum;
      //printf("the first loop: sum[%d] = %f\n", i, sum);
      A[i][r] = A[i][r] - sum; 
      sum = 0.0;
    }
    /*for ( i = 0; i < n; i++ )
    {
      printf(" check_sum[%d] = %f\n", i, check_sum[i]);
    }*/
    check_max = check_sum[r];
    for ( k = r+1; k < n; k++ )
    {
      if ( Abs(check_sum[k]) > Abs(check_max) )
      {
	check_max = check_sum[k];
	record_order[r] = k;
      }
    }
    // exchange the row of r and record_order[r]
    if ( record_order[r] != r )
    {
      //printf(" row interchanges!! ");
      record_order_temp = record_order[r];
      for ( i = 0; i < n; i++ )
      {
	temp = A[r][i];
	A[r][i] = A[record_order_temp][i];
	A[record_order_temp][i] = temp;
      }
      for ( i = 0; i < r; i++ )
      {
	temp = L[r][i];
	L[r][i] = L[record_order_temp][i];
	L[record_order_temp][i] = temp;
      }
    }
    // get the r row of U and the r column of L
    U[r][r] = A[r][r];
    //printf("check_max = %f\n", U[r][r]);
    for ( i = r+1; i < n; i++ ) 
    {
      L[i][r] = A[i][r] / U[r][r];
      //printf(" L[%d][%d] = %f\n", i, r, L[i][r]);
      for ( k = 0; k < r; k++ )
      { 
	sum_U += L[r][k]*U[k][i];
	//printf(" sum_U = %f\n", sum_U);
      }
      U[r][i] = A[r][i] - sum_U;
      //printf(" U[%d][%d] = %f\n", r, i, U[r][i]);
      sum_U = 0.0;
    }    
  }
  //printf("The element of L: ");
  //printf("the value of sum_nonz_L is: %d\n", sum_nonz_L);
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
      fprintf(fp_U, "U[%d][%d]=%lf\n", i, j, U[i][j]);
      if ( U[i][j] != 0 )
	sum_nonz_U++;
    }
  }
  U_ratio = (1 - sum_nonz_U/sum_element)*100;
  //printf("sum_nonz_U = %d\n", sum_nonz_U);
  //printf("The sparsity of U is: %lf%%\n", U_ratio);
  //Ly = b Ux = y
  y[0] = A[0][n];  
  for ( i = 1; i < n; i++ )
  {
    for ( k = 0; k < i; k++ )
    {
      sum_y += L[i][k]*y[k]; 
    }
    y[i] = A[i][n] - sum_y;
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
    //printf("x[%d]=%lf ", i, x[i]);
    fprintf(result, "x[%d]=%lf\n", i, x[i]);
  }
  for ( i = 0; i < n; i++ )
  {
    //free(A[i]);
    free(L[i]);
    free(U[i]);
  }
  //free(x);
  free(y);
  free(check_sum);
  free(record_order);
}

