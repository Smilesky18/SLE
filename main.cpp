#include "lib/lu.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <slu_ddefs.h>
#define MICRO_IN_SEC 1000000.00

double microtime(){
        int tv_sec,tv_usec;
        double time;
        struct timeval tv;
        struct timezone tz;
        gettimeofday(&tv,&tz);

        return tv.tv_sec+tv.tv_usec/MICRO_IN_SEC;
}

int main( int argc, char *argv[] )
{
    FILE *fp;
    double *x;
    double **A;
    double start, finish;
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
    double sum_time;
    fp = fopen(argv[1], "r");
    dreadMM(fp, &m, &n, &nnz, &a, &asub, &xa);
    readMatrix(argv[1], &h_val, &h_cols, &h_rowDelimiters,&nItems, &numRows, &numCols);
    printf(" nItems = %d\n ", nItems);
    printf(" numRows = %d\n ", numRows);
    printf(" numCols = %d\n ", numCols);
    x = ( double * )malloc(sizeof(double) * numCols);
    printf("\n%d %d %d\n", numRows, numCols, nItems);
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
    for ( i = 0; i < 10; i++ )
    {
      start = microtime();
      lu_sparse( h_val, h_cols, h_rowDelimiters, x, numRows);
      finish = microtime() - start; 
      sum_time += finish;
    }
    printf( "The average cost time of lu_DENSE function is: %f seconds\n", sum_time/10 );  
    for ( i = 0; i < 10; i++ )
    {
      start = microtime();
      lu_dense(A, x, n);
      finish = microtime() - start; 
      sum_time += finish;
    }
    printf( "The average cost time of lu_SPARSE function is: %f seconds\n", sum_time/10 ); 
    printf("Normal end of execution");
    free(x); 
    return 0;
}
