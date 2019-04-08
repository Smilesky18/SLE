

/*! @file lu.h
 * \brief Header file for real operations
 * 
 * <pre> 
 * </pre>
 */

/*
 * File name:		lu.h
 * Purpose:             Sparse matrix types and function prototypes
 * History:
 */

#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/time.h>
#include <slu_ddefs.h>


/*! \brief Driver routines */
extern double* lu_dense(double **, double *, int );
extern double* lu_sparse(double *, int *, int *, int );
extern double* lu_sparse_column(double *, int *, int *, int );
extern double* super_lu(char *, double *);
extern void readMatrix(char *, double **, int **, int **, int *, int *, int *);
extern double Abs(double );
extern double microtime();
//extern void super_lu(FILE *);



