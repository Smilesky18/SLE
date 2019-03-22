#include "lib/lu.h"
#include <stdio.h>
#include <stdlib.h>

int main( int argc, char *argv[] )
{
    FILE *fp;
    double *x;
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
    fp = fopen(argv[1], "r");
    dreadMM(fp, &m, &n, &nnz, &a, &asub, &xa);
    readMatrix(argv[1], &h_val, &h_cols, &h_rowDelimiters, &nItems, &numRows, &numCols);
    x = ( double * )malloc(sizeof(double) * numCols);
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
    for ( i = 0; i < 10; i++ )
    {
      start_super = microtime();
      super_lu( argv[1] );
      finish_super = microtime() - start_super;
      time_super += finish_super;
    }
    printf("The average cost time of super_lu function is: %f seconds\n", time_super/10);
    for ( i = 0; i < 10; i++ )
    {
      start_sparse = microtime();
      lu_sparse( h_val, h_cols, h_rowDelimiters, x, numRows);
      finish_sparse = microtime() - start_sparse; 
      time_sparse += finish_sparse;
    }
    printf("The average cost time of lu_SPARSE function is: %f seconds\n", time_sparse/10);  
    for ( i = 0; i < 10; i++ )
    {
      start_dense = microtime();
      lu_dense(A, x, n);
      finish_dense = microtime() - start_dense; 
      time_dense += finish_dense;
    }
    printf("The average cost time of lu_DENSE function is: %f seconds\n", time_dense/10); 
    printf("Normal end of execution");
    free(x); 
    return 0;
}
