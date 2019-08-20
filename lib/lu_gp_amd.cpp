# include <stdio.h>
# include <stdlib.h>
# include "lu.h"
# include "amd.h"

double* lu_gp_amd(double *a, int *asub, int *xa, int n)
{
  double **L, **U, *xx;
  double sum_y = 0.0, sum_x = 0.0;
  double *y, *x;
  double *check_sum, sum = 0.0;
  double check_max; 
  int i, j, r, k;
  int change_k, change_record_order;
  int *p, *row_num;
  int *P, *Pinv, *per_asub;
  int nnz = xa[n];
  double *Control, *Info;
  double *l, *u, l_value;
  double temp = 0.0;
  int record_order_temp;
  FILE *fp_L = fopen("result/Sparse_GP_L.txt", "w");
  FILE *fp_U = fopen("result/Sparse_GP_U.txt", "w");
  FILE *result = fopen("result/Sparse_GP_solution.txt", "w");
  FILE *pivot = fopen("result/Sparse_GP_pivot.txt", "w");
  FILE *per_array = fopen("per_array.txt", "w");
  FILE *pinv_file = fopen("pinv.txt", "w");
  L = ( double ** )malloc(sizeof(double *) * n);
  U = ( double ** )malloc(sizeof(double *) * n);
  P = ( int *)malloc(sizeof(int) * n );
  Pinv = ( int *)malloc(sizeof(int) * n );
  per_asub= ( int *)malloc(sizeof(int) * nnz );
  Control = ( double * )malloc(sizeof(double) * 5);
  Info = ( double * )malloc(sizeof(double) * 20);
  int sum_nonz_L = 0, sum_nonz_U = 0, sum_pviot_num = 0;
  double sum_element = n * n;
  double L_ratio, U_ratio;
  double start, finish, start_pivot, finish_pivot, sum_pivot = 0.0, start_cal, finish_cal, sum_cal = 0.0;
  int m = n + 1;
  amd_defaults (Control) ;
  amd_control  (Control) ;
  amd_order (n, xa, asub, P, Control, Info);
  printf(" asub[0] = %d\n ", asub[0]);
  printf(" asub[%d] = %d\n ", xa[n-1], asub[xa[n-1]]);
  for (k = 0 ; k < n ; k++)
  {
	j = P [k] ;
	Pinv [j] = k ;
  }
  printf(" xa[n] = %d\n ", xa[n]);
  for ( i = 0; i < n; i++ )
  {
      fprintf(pinv_file, "%d\n", Pinv[i]);
}
  for ( i = 0; i < nnz; i++ )
  {
      per_asub[i] = Pinv[asub[i]];
}
  for ( i = 0; i < nnz; i++ )
  {
      fprintf(per_array, "%d\n", per_asub[i]);
      printf("%d\n ", Pinv[asub[i]]+1);
   //   fprintf(per_array, " this is %d row\n ", i);
  }
  printf( " i = %d\n ", i);
  printf(" Info [AMD_LNZ] = %lf\n ", 2*Info [AMD_LNZ]);
//   for ( i = 0; i < n; i++ )
//   {
//     printf(" P[%d] = %d\n ", i, P[i]);
//   }
                                  
//   for ( i = 0; i < n; i++ )
//   {
//     L[i] = ( double * )malloc(sizeof(double) * n);
//     U[i] = ( double * )malloc(sizeof(double) * n);
//   }
//   for ( i = 0; i < n; i++ )
//     {
//       for ( j = 0; j < n; j++ )
//       {
// 	L[i][j] = 0.0;
// 	U[i][j] = 0.0;
//       }
//     }
//   for ( i = 0; i < n; i++ )
//   {
//     L[i][i] = 1.0;
//   }
//   xx = ( double *)malloc(sizeof(double) * n );
//   x = ( double *)malloc(sizeof(double) * m );
//   row_num = ( int *)malloc(sizeof(int) * n );
//   for ( i = 0; i < n; i++ )
//   {
//     row_num[i] = i;
//   }
//   // column-oriented G/P algorithm with partial pivoting
//   start = microtime();
//   for ( k = 0; k < n; k++ )
//   {
//     record_order_temp = k;
//     // xx = A[;k]
//     for ( j = xa[k]; j < xa[k+1]; j++ )
//     {
//       xx[row_num[asub[j]]] = a[j];
//     }
//     start_cal = microtime();
//     for ( j = 0; j < k; j++ )
//     {
//       for ( i = j+1; i < n; i++ )
//       {
// 	xx[i] -= xx[j]*L[i][j];
//       }
//     }
//     finish_cal = microtime() - start_cal;
//     sum_cal += finish_cal;
//     // find the max value in a[k:n,k]
//     temp = xx[k];
//     for ( i = k+1; i < n; i++ )
//     {
//       if ( Abs(xx[i]) > Abs(temp) )
//       {
// 	temp = xx[i];
// 	record_order_temp = i;
//       }
//     }
//     // change the k'th and max value
//     start_pivot = microtime();
//     if ( record_order_temp != k )
//     {
//       sum_pviot_num ++;
//       fprintf(pivot, " %d %d\n ", record_order_temp, k);
//       for ( i = 0; i < n; i++ )
//       {
// 	if ( row_num[i] == k )
// 	{
// 	  change_k = i;
// 	}
// 	if ( row_num[i] == record_order_temp )
// 	{
// 	  change_record_order = i;
// 	}
//       }
//       // exchange row_num[change_k] and row_num[change_record_order]
//       row_num[change_k] = record_order_temp;
//       row_num[change_record_order] = k;
//       // change the k'th and max value in xx
//       temp = xx[k];
//       xx[k] = xx[record_order_temp];
//       xx[record_order_temp] = temp;
//       // change the k'th and max value in L
//       for ( i = 0; i < k; i++ )
//       {
// 	temp = L[k][i];
// 	L[k][i] = L[record_order_temp][i];
// 	L[record_order_temp][i] = temp;
//       }
//     }
//     finish_pivot = microtime() - start_pivot;
//     sum_pivot += finish_pivot;
//     //printf(" pivot time is: %lf\n ", finish_pivot);
//     U[k][k] = xx[k];
//     for ( i = 0; i < k; i++ )
//     {
//       U[i][k] = xx[i];
//     }
//     for ( i = k+1; i < n; i++ )
//     {
//       L[i][k] = xx[i] / xx[k];
//     }
//     for ( i = 0; i < n; i++ )
//     {
//       xx[i] = 0;
//     }
//   }
//   finish = microtime() - start;
//   printf("The time of LU_column decomposition is %lf\n", finish);
//   //printf(" The number of pvioting is %d\n", sum_pviot_num);
//   printf(" sum pivot time is: %lf\n ", sum_pivot);
//   printf(" sum calculate time is: %lf\n ", sum_cal);	
//   y = ( double * )malloc(sizeof(double) * n);
//   for ( i = 0; i < n; i++ )
//   {
//     for( j = 0; j < n; j++ )
//     {
//       //printf("L[%d][%d]=%lf\n", i, j, L[i][j]);
//       if ( L[i][j] != 0 )
//       {
// 	sum_nonz_L ++;
// 	fprintf(fp_L, "L[%d][%d]=%lf\n", i, j, L[i][j]);
//       }
//       
//     }
//   }
//   for ( i = 0; i < n; i++ )
//   {
//     for( j = 0; j < n; j++ )
//     {
//       if ( U[i][j] != 0 )
//       {
// 	sum_nonz_U ++;
// 	fprintf(fp_U, "U[%d][%d]=%f\n", i, j, U[i][j]);
//       }
//       //printf("U[%d][%d]=%lf\n", i, j, U[i][j]);
//       
//     }
//   }
//   printf(" #nnz of L is: %d   #nnz of U is: %d\n", sum_nonz_L, sum_nonz_U);
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
//   x[n] = sum_pivot;
  return x;
}

