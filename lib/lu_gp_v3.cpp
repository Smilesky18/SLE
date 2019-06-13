# include <stdio.h>
# include <stdlib.h>
# include "lu.h"

double* lu_gp_v3(double *a, int *asub, int *xa, int n)
{
  double *L_a,  *U_a, *xx, *U;
  int *L_asub, *L_xa, *U_asub, *U_xa;
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
  int lu_nnz = n*(n-1)/2;
  FILE *fp_L = fopen("result/Sparse_GP_L_v3.txt", "w");
  FILE *fp_U = fopen("result/Sparse_GP_U_v3.txt", "w");
  FILE *result = fopen("result/Sparse_GP_solution_v3.txt", "w");
  FILE *pviot = fopen("result/Sparse_GP_pviot_v3.txt", "w");
  L_a = (double *)malloc(sizeof(double) * lu_nnz);
  L_asub = (int *)malloc(sizeof(int) * lu_nnz);
  L_xa = (int *)malloc(sizeof(int) * n);
  U_a = (double *)malloc(sizeof(double) * lu_nnz);
  U_asub = (int *)malloc(sizeof(int) * lu_nnz);
  U_xa = (int *)malloc(sizeof(int) * n);
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
  NCformat *Astore;
  NCPformat *ACstore;
  Astore = AA.Store;
  double *value = ( double *) Astore->nzval;
  int *A_rowidx = Astore->rowind;
  int *A_colptr = Astore->colptr;
  superlu_options_t options;
  set_default_options(&options);
  options.Fact = DOFACT;
  get_perm_c(3, &AA, perm_c);
  etree = intMalloc(AA.ncol);
  sp_preorder(&options, &AA, perm_c, etree, &AC);
  printf(" etree[4] = %d\n ", etree[4]);
  ACstore = AC.Store;
  int *AC_rowidx = ACstore->rowind;
  int *AC_colbeg = ACstore->colbeg;
  int *AC_colend = ACstore->colend;
  double *perm_value = ( double * ) ACstore->nzval;
  xx = ( double *)malloc(sizeof(double) * n );
  x = ( double *)malloc(sizeof(double) * n );
  U = ( double *)malloc(sizeof(double) * n );
  x_real= ( double *)malloc(sizeof(double) * n );
  row_num = ( int *)malloc(sizeof(int) * n );
  
  for ( i = 0; i < n; i++ )
  {
    row_num[i] = i;
    xx[i] = 0;
    L_xa[i] = 0;
    U_xa[i] = 0;
    printf("perm_c[%d] = %d\n ", i, perm_c[i]);
  } 
  for ( i = 0; i < lu_nnz; i++ )
  {
    L_asub[i] = 0;
    U_asub[i] = 0;
  }
  // column-oriented G/P algorithm with partial pivoting
  //int col_num;
  int L_col_num, L_col, num = 0, ii;
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
      L_col_num = L_xa[j+1] - L_xa[j];
      L_col = L_xa[j];
      if ( xx[j] != 0 )
      {
	for ( i = 0; i < L_col_num; i++ )
	{
	//xx[i] -= xx[j]*L[i][j];
	  ii = L_col + i;
	  xx[row_num[L_asub[ii]]] = xx[row_num[L_asub[ii]]] - xx[j] * L_a[ii];
	}
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
      // L的换主元是个问题，因为只换一部分，不是换全行，所以没有办法用指针数组去存
    }
    U[k] = xx[k];
    num = U_xa[k];
    for ( i = 0; i <= k; i++ )
    {
      if ( xx[i] != 0 )
      {
	U_asub[num] = i; 
	U_a[num] = xx[i];
	//U_xa[k+1] ++;
	num ++;
      }
    }
    U_xa[k+1] = num;
    num = L_xa[k];
    for ( i = k+1; i < n; i++ )
    {
      if ( xx[i] != 0 )
      {
	L_asub[num] = i;
	L_a[num] = xx[i];
	num ++;
      }
    }
    L_xa[k+1] = num;
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
  y = ( double * )malloc(sizeof(double) * n);
  printf(" #nnz of L is: %d   #nnz of U is: %d\n", L_xa[n], U_xa[n]);
  //y[0] = 1.0;  
  int U_col_num, U_col;
  for ( i = 0; i < n; i++ )
  {
    y[i] = 1.0;
  }
  for ( i = 1; i < n; i++ )
  {
    L_col_num = L_xa[i] - L_xa[i-1];
    L_col = L_xa[i-1];
    for ( j = 0; j < L_col_num; j++ )
    {
      y[row_num[L_asub[L_col+j]]] = y[row_num[L_asub[L_col+j]]] - y[i-1] * L_a[L_col+j];
    }
  }
  x[n-1] = y[n-1] / U_a[U_xa[n] - 1];
  for ( i = n; i >= 0; i-- )
  {
    U_col_num = U_xa[i] - U_xa[i - 1];
    U_col = U_xa[i - 1];
    for ( j = 0; j < U_col_num; j++ )
    {
      y[U_asub[U_col+j]] =  y[U_asub[U_col+j]] - y[i-1] * U_a[U_col+j];
    }
  }
  for ( i = 0; i < n; i++ )
  {
    x[i] = y[i] / U[i];
  }
  for ( i = 0; i < n; i++ )
  {
    x_real[i] = x[perm_c[i]];
    //printf("x[%d] = %f\n", i, x_real[i]);
    fprintf(result, "x[%d]=%lf\n", i, x_real[i]);
  }
  return x_real;
}

