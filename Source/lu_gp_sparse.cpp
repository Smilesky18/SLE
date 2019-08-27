# include <stdio.h>
# include <stdlib.h>
# include "../Include/lu.h"
# include <sys/time.h>
# include <float.h>
# define MICRO_IN_SEC 1000000.00

/* Time Stamp */
double microtime()
{
        struct timeval tv;
        struct timezone tz;
        gettimeofday(&tv,&tz);

        return tv.tv_sec+tv.tv_usec/MICRO_IN_SEC;
}
/* a ?= b */
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
/* return |x| */
double Abs(double x)
{
  return x < 0 ? -x : x;
}

void* lu_gp_sparse(double *a, int *asub, int *xa, int n, int nzl, int nzu, int *perm_c, int *perm_r, int *asub_L, int *xa_L, int *asub_U, int *xa_U, double *l_data, double *u_data)
{
  int sum_l = 0, sum_u = 0, row;
  double *L, *U, *xx;
  int *row_index, row_column;
  double U_diag;
  int i, j, k, current_column;
  FILE *fp_L_Error = fopen("Error_Result/L_Error.txt", "w");
  FILE *fp_U_Error = fopen("Error_Result/U_Error.txt", "w");
  L = ( double * )malloc(sizeof(double *) * nzl);
  U = ( double * )malloc(sizeof(double *) * nzu);
  double start_cal, finish_cal, sum_cal = 0.0, start_xx_0, end_xx_0, sum_xx_0 = 0.0;
  xx = ( double *)malloc(sizeof(double) * n );
  int *sum_calculation, sum = 0;
  sum_calculation = ( int *)malloc(sizeof(int *) * n);
  memset(sum_calculation, 0, sizeof(int)*n);
  
  /* Array xx initialization*/
  for ( i = 0; i < n; i++ )
  {
    xx[i] = 0;
    L[xa_L[i]] = 1.0;
  }

  /* column-oriented G/P algorithm without partial pivoting */
  for ( k = 0; k < n; k++ )
  {
    current_column = perm_c[k];
    
    /* xx[] = A[:,current_column] */
    for ( j = xa[current_column]; j < xa[current_column+1]; j++ )
    {
      xx[perm_r[asub[j]]] = a[j];
    }
    row_column = xa_U[k+1] - xa_U[k] - 1;
    row_index = (int *)malloc( sizeof(int *) * row_column);
    for ( j = 0; j < row_column; j++ )
    {
        row_index[j] = asub_U[j+xa_U[k]];
    }
    
    /* L[:,0~k]*xx = A[:,current_column], solve for xx*/
    start_cal = microtime();
    for ( j = 0; j < row_column; j++ )
    {
        row = row_index[j];
        for ( i = xa_L[row]+1; i < xa_L[row+1]; i++ )
        {
            xx[asub_L[i]] -=  xx[row]*L[i];
            sum_calculation[k]++;
        }
    }
    finish_cal = microtime() - start_cal;
    sum_cal += finish_cal;

    /* solve for U[:,k]*/
//     U[xa_U[k]] = xx[k];
    for ( i = xa_U[k]; i < xa_U[k+1]; i++ )
    {
      U[i] = xx[asub_U[i]];
      xx[asub_U[i]] = 0;
    }

    /* solve for L[:,k] */
    U_diag = U[i-1];
//     if ( equal(U_diag, 0.0) ) printf("%d\n", k);
    for ( i = xa_L[k]+1; i < xa_L[k+1]; i++ )
    {
      L[i] = xx[asub_L[i]] / U_diag;
      xx[asub_L[i]] = 0;
    }
  }

  
  /* Check L value*/
  for ( i = 0; i < nzl; i++ )
  {
      if ( equal(L[i], l_data[i]) ) continue;
      else 
      {
          sum_l++; 
          fprintf(fp_L_Error, "Correct: %lf Error: %lf\n", l_data[i], L[i]);
          
    }
      
  }
  
  
  /* Check U value*/
  for ( i = 0; i < nzu; i++ )
  {
      if ( equal(u_data[i], U[i]) ) continue;
      else
      {
          sum_u++;
          fprintf(fp_U_Error, "Correct: %lf Error: %lf\n", u_data[i], U[i]);
          
      }
  }


  /* Print check results information*/
  if ( sum_l == 0 && sum_u == 0 ) printf("Correct Results!\n");
  else 
  {
      printf("Num of Errors in L: %d\n", sum_l);
      printf("Num of Errors in U: %d\n", sum_u);
      printf("Notice: For Error Results, please check the detailed information in Error_Result directory!\n");
  }
  
  
  /* Print Time Information */
  printf("\n********************MY LU Decomposition: Time Statisistics********************\n");
  printf("Calculation Time: %lf\n", sum_cal);
}

