# include <stdio.h>
# include <stdlib.h>
# include "lu.h"

double* lu_gp_v2(double *a, int *asub, int *xa, int n)
{
  double **L, **U, *xx;
  double sum_y = 0.0, sum_x = 0.0;
  double *y, *x, *x_real;
  double *check_sum, sum = 0.0;
  double check_max; 
  int i, j, r, k;
  int change_k, change_record_order;
  int *p, *row_num;
  double *l, *u, l_value;
  double temp = 0.0;
  int record_order_temp;
  FILE *fp_L = fopen("result/Sparse_GP_L_v2.txt", "w");
  FILE *fp_U = fopen("result/Sparse_GP_U_v2.txt", "w");
  FILE *result = fopen("result/Sparse_GP_solution_v2.txt", "w");
  FILE *pviot = fopen("result/Sparse_GP_pviot_v2.txt", "w");
  L = ( double ** )malloc(sizeof(double *) * n);
  U = ( double ** )malloc(sizeof(double *) * n);
  int sum_nonz_L = 0, sum_nonz_U = 0, sum_pviot_num = 0;
  double sum_element = n * n;
  double L_ratio, U_ratio;
  double start, finish;
  int nnz = xa[n], *perm_c, *etree;
  SuperMatrix AA, AC;
  if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
  dCreate_CompCol_Matrix(&AA, n, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
  int nrow = AA.nrow;
  int ncol = AA.ncol;
  //printf(" nrow = %d\n ncol = %d\n", nrow, ncol);
  NCformat *Astore;
  NCPformat *ACstore;
  Astore = AA.Store;
  //double *value = (double *)Astore->nzval;
  //double *value = (double*) ((NCformat*) AA.Store)->nzval;
  double *value = ( double *) Astore->nzval;
  int *A_rowidx = Astore->rowind;
  int *A_colptr = Astore->colptr;
  superlu_options_t options;
  set_default_options(&options);
  options.Fact = DOFACT;
  /*for ( i = 0; i < nnz; i++ )
  {
    //printf("value[%d] = %lf\n", i, value[i]);
    printf(" row_index[%d] = %d\n", i, A_rowidx[i]);
  }
  for ( i = 0; i <= n; i++ )
  {
    printf(" col_ptr[%d] = %d\n", i, A_colptr[i]);
  }*/
  get_perm_c(3, &AA, perm_c);
  //sp_preorder(&options, &AA, perm_c, etree, &AC);
  /*for ( i = 0; i < n; i++ )
  {
    printf(" perm_c[%d] = %d\n", i, perm_c[i]);
  }*/
  etree = intMalloc(AA.ncol);
  /*if ( options.Fact == DOFACT)
  {
    printf(" It is DOFACT!! ");
  }*/
  sp_preorder(&options, &AA, perm_c, etree, &AC);
  //printf(" AC->nrow = %d\n ", AC.nrow);
  ACstore = AC.Store;
  int *AC_rowidx = ACstore->rowind;
  int *AC_colbeg = ACstore->colbeg;
  int *AC_colend = ACstore->colend;
  double *perm_value = ( double * ) ACstore->nzval;
  /*for ( i = 0; i < nnz; i++ )
  {
    //printf("perm_value[%d] = %lf\n", i, perm_value[i]);
    printf("AC_rowidx[%d] = %d\n", i, AC_rowidx[i]);
  }*/
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
  for ( i = 0; i < n; i++ )
  {
    L[i][i] = 1.0;
  }
  xx = ( double *)malloc(sizeof(double) * n );
  x = ( double *)malloc(sizeof(double) * n );
  x_real= ( double *)malloc(sizeof(double) * n );
  row_num = ( int *)malloc(sizeof(int) * n );
  for ( i = 0; i < n; i++ )
  {
    row_num[i] = i;
  }
  // column-oriented G/P algorithm with partial pivoting
  int col_num;
  start = microtime();
  for ( k = 0; k < n; k++ )
  {
    record_order_temp = k;
    // xx = A[;k]
    //col_num = ACstore->colend[k] - ACstore->colbeg[k];
    for ( j = ACstore->colbeg[k]; j < ACstore->colend[k]; j++ )
    {
      xx[row_num[asub[j]]] = a[j];
    }
    for ( j = 0; j < k; j++ )
    {
      for ( i = j+1; i < n; i++ )
      {
	xx[i] -= xx[j]*L[i][j];
      }
    }
    // find the max value in a[k:n,k]
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
      sum_pviot_num ++;
      for ( i = 0; i < n; i++ )
      {
	if ( row_num[i] == k )
	{
	  change_k = i;
	}
	if ( row_num[i] == record_order_temp )
	{
	  change_record_order = i;
	}
      }
      // exchange row_num[change_k] and row_num[change_record_order]
      row_num[change_k] = record_order_temp;
      row_num[change_record_order] = k;
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
      //printf(" k = %d\n ", k);
      //printf(" row_num[%d] = %d ", i, row_num[i]);
    }
    //printf(" \n ");
  }
  finish = microtime() - start;
  //printf(" before the LU decomposition ");
  printf("The time of LU_column decomposition is %lf\n", finish);
  printf(" The number of pvioting is %d\n", sum_pviot_num);
  y = ( double * )malloc(sizeof(double) * n);
  for ( i = 0; i < n; i++ )
  {
    for( j = 0; j < n; j++ )
    {
      //printf("L[%d][%d]=%lf\n", i, j, L[i][j]);
      if ( L[i][j] != 0 )
      {
	sum_nonz_L ++;
      }
      fprintf(fp_L, "L[%d][%d]=%lf\n", i, j, L[i][j]);
    }
  }
  for ( i = 0; i < n; i++ )
  {
    for( j = 0; j < n; j++ )
    {
      //printf("U[%d][%d]=%lf\n", i, j, U[i][j]);
      if ( U[i][j] != 0 )
      {
	sum_nonz_U ++;
      }
      fprintf(fp_U, "U[%d][%d]=%f\n", i, j, U[i][j]);
    }
  }
  printf(" #nnz of L is: %d   #nnz of U is: %d\n", sum_nonz_L, sum_nonz_U);
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
    x_real[i] = x[perm_c[i]];
    //printf("x[%d] = %f\n", i, x_real[i]);
    fprintf(result, "x[%d]=%lf\n", i, x_real[i]);
  }
  return x_real;
}

