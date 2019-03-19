# include "slu_ddefs.h"
# include <stdio.h>
#include <time.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#define MICRO_IN_SEC 1000000.00

double Abs(double x)
{
  return x < 0 ? -x : x;
}
double microtime(){
        int tv_sec,tv_usec;
        double time;
        struct timeval tv;
        struct timezone tz;
        gettimeofday(&tv,&tz);

        return tv.tv_sec+tv.tv_usec/MICRO_IN_SEC;
}

struct Coordinate {
    int x;
    int y;
    double val;
};
int coordcmp(const void *v1, const void *v2)
{
    struct Coordinate *c1 = (struct Coordinate *) v1;
    struct Coordinate *c2 = (struct Coordinate *) v2;

    if (c1->x != c2->x)
    {
        return (c1->x - c2->x);
    }
    else
    {
        return (c1->y - c2->y);
    }
}
void readMatrix(char *filename, double **val_ptr, int **cols_ptr, 
                int **rowDelimiters_ptr, int *n, int *numRows, int *numCols) 
{
    std::string line;
    char id[128];
    char object[128]; 
    char format[128]; 
    char field[128]; 
    char symmetry[128]; 
    //FILE *mfs = fopen( filename, "r");

    std::ifstream mfs( filename );
    if( !mfs.good() )
    {
        //std::cerr << "Error: unable to open matrix file " << filename << std::endl;
        printf(" Error: unable to open matrix file ");
        exit( 1 );
    }

    int symmetric = 0; 
    int pattern = 0; 
    int field_complex = 0;
    int nRows, nCols, nElements;  
    struct Coordinate *coords;
    // read matrix header
    if( getline( mfs, line ).eof() )
    {
        //std::cerr << "Error: file " << filename << " does not store a matrix" << std::endl;
        printf(" Error: file does not store a matrix ");
        exit( 1 );
    }
    sscanf(line.c_str(), "%s %s %s %s %s", id, object, format, field, symmetry); 
    if (strcmp(object, "matrix") != 0) 
    {
        fprintf(stderr, "Error: file %s does not store a matrix\n", filename); 
        exit(1); 
    }
  
    if (strcmp(format, "coordinate") != 0)
    {
        fprintf(stderr, "Error: matrix representation is dense\n"); 
        exit(1); 
    } 

    if (strcmp(field, "pattern") == 0) 
    {
        pattern = 1; 
    }

	if (strcmp(field, "complex") == 0)
	{
        field_complex = 1;
	}

    if (strcmp(symmetry, "symmetric") == 0) 
    {
        symmetric = 1; 
    }
    while (!getline( mfs, line ).eof() )
    {
        if (line[0] != '%') 
        {
            break; 
        }
    } 
    sscanf(line.c_str(), "%d %d %d", &nRows, &nCols, &nElements); 
    int nElements_padding = (nElements%16 == 0) ? nElements : (nElements + 16)/ 16 * 16;
    int valSize = nElements * sizeof(struct Coordinate);

    if (symmetric) 
    {
        valSize*=2; 
    }                 
    coords = (struct Coordinate*)malloc(valSize);
    int index = 0; 
    float xx99 = 0;
    while (!getline( mfs, line ).eof() )
    {
        if (pattern) 
        {
            sscanf(line.c_str(), "%d %d", &coords[index].x, &coords[index].y); 
            coords[index].val = index%13;
        }
        else if(field_complex)
        {
            // read the value from file
            sscanf(line.c_str(), "%d %d %lf %f", &coords[index].x, &coords[index].y, 
                   &coords[index].val, xx99); 
        }
        else 
        {
            // read the value from file
            sscanf(line.c_str(), "%d %d %lf", &coords[index].x, &coords[index].y, 
                   &coords[index].val); 
        } 
        index++; 
        if (symmetric && coords[index-1].x != coords[index-1].y) 
        {
            coords[index].x = coords[index-1].y; 
            coords[index].y = coords[index-1].x; 
            coords[index].val = coords[index-1].val; 
            index++;
        }

    }  

    nElements = index; 
  
    nElements_padding = (nElements%16 == 0) ? nElements : (nElements + 16)/ 16 * 16;
 
  /*printf("===========================================================================\n");
  printf("=========*********  Informations of the sparse matrix   *********==========\n");
  printf(" Number of Rows is %d\n:", nRows);
  printf(" Number of Columns is %d\n:", nCols);
  printf(" Number of Elements is %d\n:", nElements);
  printf(" After Alignment : %d\n", nElements_padding);
  printf("===========================================================================\n");
  printf("............ Converting the Raw matrix to CSR .................\n");*/ 
  /*for (int qq = index; qq <nElements_padding; qq ++)
  {
            coords[qq].x = coords[index - 1].x;
            coords[qq].y = coords[index - 1].y;
            coords[qq].val = 0;      
  }*/
  for ( int i = 0; i < nElements; i++ )
  {
    coords[i].x = coords[i].x - 1;
    coords[i].y = coords[i].y - 1;
  }
  qsort(coords, nElements, sizeof(struct Coordinate), coordcmp); 
  /*for ( int i = 0; i < nElements; i++ )
  {
    printf(" coords[%d].x = %d\n", i, coords[i].x);
    printf(" coords[%d].y = %d\n", i, coords[i].y);
    printf(" coords[%d].val = %f\n", i, coords[i].val);
  }*/


    // create CSR data structures
    *n = nElements; 
    *numRows = nRows; 
    *numCols = nCols; 
//    *val_ptr = new floatType[nElements_padding]; 
//    *cols_ptr = new int[nElements_padding];
//    *rowDelimiters_ptr = new int[nRows+2]; 

    *val_ptr = (double *)malloc(sizeof(double) * nElements); 
    *cols_ptr = (int *)malloc(sizeof(int) * nElements);
    *rowDelimiters_ptr = (int *)malloc(sizeof(int) * (nRows+2));

    double *val = *val_ptr; 
    int *cols = *cols_ptr; 
    int *rowDelimiters = *rowDelimiters_ptr; 
    rowDelimiters[0] = 0; 
    int r=0; 
    int i=0;
    for (i=0; i<nElements; i++) 
    {
        while (coords[i].x != r) 
        {
            rowDelimiters[++r] = i; 
        }
        val[i] = coords[i].val; 
        cols[i] = coords[i].y;    
    }
    for(int k = r + 1; k<=(nRows+1); k++)
    {
       rowDelimiters[k] = i; 
    }
    r = 0; 
    free(coords); 
}
void LU(double *a, int *asub, int *xa, double *x, int n)
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
  FILE *fp_L = fopen("L_CSC.txt", "w");
  FILE *fp_U = fopen("U_CSC.txt", "w");
  FILE *result = fopen("result_CSC.txt", "w");
  double lu_start, lu_finish, trsv_start, trsv_finish;
  double lu_duration, trsv_duration;
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
    //printf(" row_num[%d] = %d\n", i, row_num[i]);
  }
  printf(" before the LU decomposition ");
  for ( r = 0; r < n; r++ )
  {
    //pivot_start = clock();
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
      //printf(" l[%d] = %lf\n", i, l[i]);
    }
    check_max = check_sum[r];
    record_order_temp = r;
    for ( k = r+1; k < n; k++ )
    {
      //if ( Abs(check_sum[k]) > Abs(check_max) )
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
    //pivot_finish = clock();
    //pivot_duration = pivot_duration + (double)(pivot_finish - pivot_start) / CLOCKS_PER_SEC;
    //printf("check_max = %f\n", U[r][r]);
    for ( j = 0; j < row_num[r]; j++ )
    {
      u[asub[p[r]+j]] = a[p[r]+j];
      //printf(" u[%d] = %lf\n", asub[p[r]+j], a[p[r]+j]);
    }
    for ( i = r+1; i < n; i++ ) 
    {
      L[i][r] = l[i] / U[r][r];
      //printf(" L[%d][%d] = %f\n", i, r, L[i][r]);
      for ( k = 0; k < r; k++ )
      { 
	sum += L[r][k]*U[k][i];
	//printf(" sum_U = %f\n", sum);
      }
      U[r][i] = u[i] - sum;
      sum = 0.0;
      //U[r][i] = -sum_U[i];
      /*for ( j = 0; j < row_num[r]; j++ )
      {
	if ( asub[p[r]+j] == i )
	{
	  U[r][i] = a[p[r]+j] - sum;
	  row_used = false;
	  break;
	}
      }
      if ( row_used )
      {
	U[r][i] = 0 - sum;
      }
      row_used = true;
      //printf(" U[%d][%d] = %f\n", r, i, U[r][i]);
      sum = 0.0;*/
    }
    for ( i = 0; i < n; i++ )
    {
      l[i] = 0.0;
      u[i] = 0.0;
    }
  }
  y = ( double * )malloc(sizeof(double) * n);
 // printf( "The cost time of LU decomposition is: %f seconds\n", lu_duration );  
  //printf("The element of L: ");
 // printf("the value of sum_nonz_L is: %d\n", sum_nonz_L);
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
  printf("sum_element = %f\n", sum_element);
  printf("sum_nonz_L = %d\n", sum_nonz_L);
  printf("The sparsity of L is: %lf%%\n", L_ratio);
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
  printf("sum_nonz_U = %d\n", sum_nonz_U);
  printf("The sparsity of U is: %lf%%\n", U_ratio);
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
  //printf( "The cost time of TRSV is: %f seconds\n", trsv_duration );
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

int main( int argc, char *argv[] )
{
    double *x;
    double start, finish;
    double *h_val;
    int *h_cols;       
    int *h_rowDelimiters;
    // Number of non-zero elements in the matrix
    int nItems;
    int numRows;
    int numCols;
    //printf(" This is only for testing! \n");
    //dreadMM(fp, &m, &n, &nnz, &a, &asub, &xa);
    readMatrix(argv[1], &h_val, &h_cols, &h_rowDelimiters,&nItems, &numRows, &numCols);
    /*printf(" vals: ");
    for ( i = 0; i < nItems; i++ )
    {
      printf(" %lf ", h_val[i]);
    }
    printf("\n");
    printf(" colIdx: ");
    for ( i = 0; i < nItems; i++ )
    {
      printf(" %d ", h_cols[i]);
    }
    printf("\n");
    printf(" offset: ");
    for ( i = 0; i < numRows+1; i++ )
    {
      printf(" %d ", h_rowDelimiters[i]);
    }
    printf("\n");*/
    printf(" nItems = %d\n ", nItems);
    printf(" numRows = %d\n ", numRows);
    printf(" numCols = %d\n ", numCols);
    //printf(" This is for testing!!! ");
    /*for ( i = 0; i < nnz; i++ )
    {
      printf("a[%d]=%f\n", i, a[i]);
      printf("asub[%d]=%d\n", i, asub[i]);
    }
    for ( i = 0; i < (n+1); i++ )
    {
      printf("xa[%d]=%d\n", i, xa[i]);
    }*/
    x = ( double * )malloc(sizeof(double) * numCols);
    printf("\n%d %d %d\n", numRows, numCols, nItems);
    start = microtime();
    LU( h_val, h_cols, h_rowDelimiters, x, numRows);
    finish = microtime() - start; 
    printf( "The cost time of LU function is: %f seconds\n", finish );  
    printf("Normal end of execution");
    free(x); 
    return 0;
}

