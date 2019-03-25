#include "lib/lu.h"
#include <stdio.h>
#include <stdlib.h>


bool equal( double a, double b )
{
  if ( Abs(a-b) < 0.0001)
  {
    return true;
  }
  else
  {
    return false;
  }
}
int main( int argc, char *argv[] )
{
    FILE *fp;
    double *x;
    double *sparse_x, *dense_x, *super_x;
    double **A;
    double start_super, start_sparse, start_dense, finish_super, finish_sparse, finish_dense;
    double time_super, time_sparse, time_dense;
    double *h_val;
    int *h_cols;       
    int *h_rowDelimiters;
    // Number of non-zero elements in the matrix
    int nItems;
    int numRows;
    int numCols;
    int NUMCol, counter_num = 0;
    double   *a;
    int      *asub, *xa;
    int m, n, nnz;
    int i, j;
    int check_solution_right_num, check_solution_wrong_num;
    fp = fopen(argv[1], "r");
    dreadMM(fp, &m, &n, &nnz, &a, &asub, &xa);
    readMatrix(argv[1], &h_val, &h_cols, &h_rowDelimiters, &nItems, &numRows, &numCols);
    x = ( double * )malloc(sizeof(double) * numCols);
    sparse_x = ( double * )malloc(sizeof(double) * numCols);
    dense_x = ( double * )malloc(sizeof(double) * numCols);
    super_x = ( double * )malloc(sizeof(double) * numCols);
    A = ( double ** )malloc(sizeof(double *) * m);  
    int ncol = n + 1;                                                                                                                                                         
    for ( i = 0; i < m; i++ )                                                                                                                                                 
    {                                                                                                                                                                         
      A[i] = ( double * )malloc(sizeof(double) * ncol);                                                                                                                       
    }                                                                                                                                                                         
    for ( i = 0; i < m; i++ )                                                                                                                                                 
    {                                                                                                                                                                         
      for ( j = 0; j < n; j++ )                                                                                                                                               
      {                                                                                                                                                                        
	A[i][j] = 0.0;                                                                                                                                                            
      }                                                                                                                                                                       
    }              
    for ( j = 0; j < n; j++ )
    {
      NUMCol = xa[j+1] - xa[j];
      for ( i = 0; i < NUMCol; i++ )
      {
	A[asub[counter_num]][j] = a[counter_num];
	counter_num++;
      }
    }
    for ( i = 0; i < n; i++ )
    {
      A[i][n] = 1.0;
    }
    printf(" m = %d  ", m);
    printf(" n = %d  ", n);
    printf(" nnz = %d\n ", nnz);
    printf("-----------------------LU_Super----------------------\n");
    /*for ( i = 0; i < 10; i++ )
    {
      start_super = microtime();
      super_lu( argv[1], x );
      finish_super = microtime() - start_super;
      time_super += finish_super;
    }*/
    start_super = microtime();
    super_x = super_lu( argv[1], x );
    finish_super = microtime() - start_super;
    time_super += finish_super;
    /*for ( i = 0; i < n; i++ )
    {
      printf("super_x[%d] = %lf\n", i , super_x[i]);
    }*/
    printf("The average cost time of super_lu function is: %f seconds\n", time_super);
    printf("-----------------------LU_Sparse----------------------\n");
    /*for ( i = 0; i < 10; i++ )
    {
      start_sparse = microtime();
      lu_sparse( h_val, h_cols, h_rowDelimiters, x, numRows );
      finish_sparse = microtime() - start_sparse; 
      time_sparse += finish_sparse;
    }*/
    start_sparse = microtime();
    sparse_x = lu_sparse( h_val, h_cols, h_rowDelimiters, x, numRows );
    finish_sparse = microtime() - start_sparse;
    time_sparse += finish_sparse;
    /*for ( i = 0; i < n; i++ )
    {
      printf("sparse_x[%d] = %lf\n", i , sparse_x[i]);
    }*/
    printf("The average cost time of lu_SPARSE function is: %f seconds\n", time_sparse); 
    printf("-----------------------LU_Dense----------------------\n");
    /*for ( i = 0; i < 10; i++ )
    {
      start_dense = microtime();
      lu_dense( A, x, n );
      finish_dense = microtime() - start_dense; 
      time_dense += finish_dense;
    }*/
    start_dense = microtime();
    dense_x = lu_dense( A, x, n );
    finish_dense = microtime() - start_dense;
    time_dense += finish_dense;
    /*for ( i = 0; i < n; i++ )
    {
      printf("dense_x[%d] = %lf\n", i , dense_x[i]);
    }*/
    printf("The average cost time of lu_DENSE function is: %f seconds\n", time_dense);
    printf(" ------------------Verify the Result-----------------\n");
    for ( i = 0; i < n; i++ )
    {
      /*if ( sparse_x[i] == super_x[i] )
      {
	printf("sparse_x[%d] is right\n", i);
      }
      if ( dense_x[i] == super_x[i] )
      {
	printf("dense_x[%d] is right\n", i);
      }*/
      if ( equal( sparse_x[i], super_x[i] ) && equal( dense_x[i], super_x[i] ))
      {
	check_solution_right_num ++;
	continue;
      }
      else
      {
	//printf(" The %dth result is wrong\n", i);
	//printf("result[%d] = %lf  sparse_x[%d] = %lf  dense_x[%d] = %lf\n", i, super_x[i], i, sparse_x[i],
	//  i, dense_x[i]
	//);
	check_solution_wrong_num ++;
	continue;
      }
    }
    if ( check_solution_right_num == n )
    {
      printf("################The solution is right!!################\n");
    }
    else
    {
      printf("The wrong number is %d\n", check_solution_wrong_num);
    }
    printf("Normal end of execution");
    free(x); 
    return 0;
}
