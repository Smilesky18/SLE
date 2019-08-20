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


double* lu_gp_sparse(double *a, int *asub, int *xa, int n, int nzl, int nzu, int *perm_c, int *perm_r, int *asub_L, int *xa_L, int *asub_U, int *xa_U, double *l_data, double *u_data)
{
  lu_data *lu_array_data;
  int sum_l = 0, sum_u = 0;
  double *L, *U, *xx;
  double sum_y = 0.0, sum_x = 0.0, U_diag;
  double *y, *x;
  double *check_sum, sum = 0.0;
  double check_max; 
  int i, j, r, k, current_column;
  int change_k, change_record_order;
  int *p, *row_num;
  double *l, *u, l_value;
  double temp = 0.0;
  int record_order_temp;
  FILE *fp_L = fopen("result/Sparse_GP_L.txt", "w");
  FILE *fp_U = fopen("result/Sparse_GP_U.txt", "w");
  FILE *result = fopen("result/Sparse_GP_solution.txt", "w");
  FILE *pivot = fopen("result/Sparse_GP_pivot.txt", "w");
  L = ( double * )malloc(sizeof(double *) * nzl);
  U = ( double * )malloc(sizeof(double *) * nzu);
  int sum_nonz_L = 0, sum_nonz_U = 0, sum_pviot_num = 0;
  double sum_element = n * n;
  double L_ratio, U_ratio;
  double start, finish, start_pivot, finish_pivot, sum_pivot = 0.0, start_cal, finish_cal, sum_cal = 0.0;
  int m = n + 1;
  xx = ( double *)malloc(sizeof(double) * n );
  x = ( double *)malloc(sizeof(double) * n );
  row_num = ( int *)malloc(sizeof(int) * n );
  
  for ( i = 0; i < n; i++ )
  {
    row_num[i] = i;
    xx[i] = 0;
    L[xa_L[i]] = 0.0;
  }
//   L[0] = 0;
  // column-oriented G/P algorithm with partial pivoting
  start = microtime();
  for ( k = 0; k < n; k++ )
  {
    current_column = perm_c[k];
//     record_order_temp = k;
    // xx = A[;k],xx中存储的是Ａ中当前列的值，按照其（置换后的）行号进行的存储
    for ( j = xa[current_column]; j < xa[current_column+1]; j++ )
    {
      xx[perm_r[asub[j]]] = a[j];
    }
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
//     U[xa_U[k]] = xx[k];
    for ( i = xa_U[k]; i < xa_U[k+1]; i++ )
    {
      U[i] = xx[asub_U[i]];
    }
    U_diag = U[i-1];
    for ( i = xa_L[k]; i < xa_L[k+1]; i++ )
    {
      L[i] = xx[asub_L[i]] / U_diag;
    }
    for ( i = 0; i < n; i++ )
    {
      xx[i] = 0;
    }
  }
//   for ( i = 0; i < n; i++ )
//   {
//     L[xa_L[i]] = 1.0;
//   }
  finish = microtime() - start;
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
    
    printf("sum_l = %d\n", sum_l);
    printf("sum_u = %d\n", sum_u);
//   printf("\nL's values are: ");
//   for ( i = 0; i < nzl; i++ ) printf("%lf ", L[i]);
//   printf("\nU's values are: ");
//   for ( i = 0; i < nzu; i++ ) printf("%lf ", U[i]);
//   printf("The time of LU_column decomposition is %lf\n", finish);
//   //printf(" The number of pvioting is %d\n", sum_pviot_num);
//   printf(" sum pivot time is: %lf\n ", sum_pivot);
//   printf(" sum pivot number time is: %d\n ", sum_pviot_num);
// //   printf(" sum calculate time is: %lf\n ", sum_cal);	
//   y = ( double * )malloc(sizeof(double) * n);
//   for ( i = 0; i < n; i++ )
//   {
//     for( j = 0; j < n; j++ )
//     {
//       //printf("L[%d][%d]=%lf\n", i, j, L[i][j]);
//       if ( L[j][i] != 0 )
//       {
//         sum_nonz_L ++;
// //         printf("L[%d][%d]=%lf\n", j, i, L[j][i]);
// // 	fprintf(fp_L, "L[%d][%d]=%lf\n", i, j, L[i][j]);
//       }
//       
//     }
//   }
//   for ( i = 0; i < n; i++ )
//   {
//     for( j = 0; j < n; j++ )
//     {
//       if ( U[j][i] != 0 )
//       {
//         sum_nonz_U ++;
//         printf("U[%d][%d]=%e\n", j, i, U[j][i]);
// // 	fprintf(fp_U, "U[%d][%d]=%f\n", i, j, U[i][j]);
//       }
//       //printf("U[%d][%d]=%lf\n", i, j, U[i][j]);
//       
//     }
//   }
//   printf(" \n#nnz of L is: %d   #nnz of U is: %d\n", sum_nonz_L, sum_nonz_U);
//   y[0] = 1.0;  
//   for ( i = 1; i < n; i++ )
//   {
//     for ( k = 0; k < i; k++ )
//     {
//       sum_y += L[i][k]*y[k]; 
//     }
//     y[i] = 1.0 - sum_y;
//     sum_y = 0.0;
//   }
//   x[n-1] = y[n-1] / U[n-1][n-1];
//   for ( i = n-2; i >= 0; i-- )
//   {
//     for ( k = i + 1; k < n; k++ )
//     {
//       sum_x += U[i][k]*x[k];
//     }
//     x[i] = ( y[i] - sum_x ) / U[i][i];
//     sum_x = 0.0;
//   }
//   for ( i = 0; i < n; i++ )
//   {
//     //printf("x[%d] = %f\n", i, x[i]);
//     fprintf(result, "x[%d]=%lf\n", i, x[i]);
//   }
 // x[n] = sum_pivot;
  return x;
}

