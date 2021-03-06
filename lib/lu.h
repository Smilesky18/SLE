

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
extern void lu_dense(double **, double *, int );
extern void lu_sparse_1(double *, int *, int *, double *, int );
extern void readMatrix(char *, double **, int **, int **, int *, int *, int *);
extern double Abs(double );
extern int super_lu(char *);
extern double microtime();
//extern void super_lu(FILE *);



