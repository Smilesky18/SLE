

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
#include <cs.h>
// #ifndef csi
// #define csi ptrdiff_t
// #endif

typedef struct lu_data_stru    /* matrix in compressed-column or triplet form */
{
    double *l_data;
    double *u_data;
} lu_data ;


/*! \brief Driver routines */
extern double* lu_dense(double **, double *, int );
extern double* lu_sparse(double *, int *, int *, int, int, int );
extern double* lu_gp_sparse(double *, int *, int *, int, int, int, int *, int *, int *, int *, int *, int * , double *, double * );
extern void* sparse_load(FILE *, int, int, int *, int *, int *, int *, int *, int * );
extern double* lu_gp(double *, int *, int *, int );
extern double* lu_gp_v2(double *, int *, int *, int );
extern double* lu_gp_v3(double *, int *, int *, int );
extern double* lu_gp_amd(double *, int *, int *, int );
// extern double* super_lu(char *, double *);
extern int super_lu(char *, double *);
extern void readMatrix(char *, double **, int **, int **, int *, int *, int *);
extern double Abs(double );
extern double microtime();
//extern void super_lu(FILE *);



