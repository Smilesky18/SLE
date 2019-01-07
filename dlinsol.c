/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 */
#include "slu_ddefs.h"
#include <stdio.h>

int main(int argc, char *argv[])
{
    SuperMatrix A;
    NCformat *Astore;
    DNformat *Bstore;
    int *colptr_B;
    int *rowind_B;
    double *nzval_B, *nzval_A, *nzval_U;
    double   *a;
    int      *asub, *xa;
    int      *perm_c; /* column permutation vector */
    int      *perm_r; /* row permutations from partial pivoting */
    SuperMatrix L;      /* factor L */
    SCformat *Lstore;
    SuperMatrix U;      /* factor U */
    NCformat *Ustore;
    SuperMatrix B;
    int      nrhs, ldx, info, m, n, nnz;
    double   *xact, *rhs;
    mem_usage_t   mem_usage;
    superlu_options_t options;
    SuperLUStat_t stat;
    FILE      *fp = stdin;
    FILE *fp_U = fopen("U_superlu.txt", "w");
    FILE *result = fopen("result_superlu.txt", "w");
    int i, j;
    int NUMCol, counter_num = 0;
    double A_Me[800][800];
    double x[800];
    
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC("Enter main()");
#endif

    /* Set the default input options:
	options.Fact = DOFACT;
        options.Equil = YES;
    	options.ColPerm = COLAMD;
	options.DiagPivotThresh = 1.0;
    	options.Trans = NOTRANS;
    	options.IterRefine = NOREFINE;
    	options.SymmetricMode = NO;
    	options.PivotGrowth = NO;
    	options.ConditionNumber = NO;
    	options.PrintStat = YES;
     */
    set_default_options(&options);

#if 0
    /* Read the matrix in Harwell-Boeing format. */
    dreadhb(fp, &m, &n, &nnz, &a, &asub, &xa);
#else
    /* Read the matrix in Matrix Market format. */
    dreadMM(fp, &m, &n, &nnz, &a, &asub, &xa);
#endif

    dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
    Astore = A.Store;
    nzval_A = Astore->nzval;
    printf("Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, Astore->nnz);
    
    nrhs   = 1;
    if ( !(rhs = doubleMalloc(m * nrhs)) ) ABORT("Malloc fails for rhs[].");
    //for ( i = 0; i < m; i++ )
    //{
    //  rhs[i] = 1.0;
    //}
    dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);
    /*for ( i = 0; i < m; i++ )
      printf(" the element in rhs is: %f\n", rhs[i]);*/
    Bstore = (DNformat *)B.Store;
    nzval_B = (double *)Bstore->nzval;
    /*for (i=0; i < Bstore->lda; i++)
      printf("B的第%d个元素是: %lf\n",i, *(nzval_B+i));*/
    xact = doubleMalloc(n * nrhs);
    ldx = n;
    dGenXtrue(n, nrhs, xact, ldx);
    dFillRHS(options.Trans, nrhs, xact, ldx, &A, &B);

    if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
    if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");

    /* Initialize the statistics variables. */
    StatInit(&stat);
    Bstore = (DNformat *)B.Store;
    nzval_B = (double *)Bstore->nzval;
    /*for (i=0; i < Bstore->lda; i++)
      printf("before dgssv B的第%d个元素是: %lf\n",i, *(nzval_B+i));*/
    dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
    Bstore = (DNformat *)B.Store;
    nzval_B = (double *)Bstore->nzval;
    //printf("after dgssv the lda of B is: %d\n", Bstore->lda);
    //for (i=0; i < Bstore->lda; i++)
    //  printf("after dgssv B的第%d个元素是: %lf\n",i, *(nzval_B+i));
    Ustore = (NCformat *)U.Store;
    nzval_U = (double *)Ustore->nzval;
    for ( i = 0; i < m*n; i++ )
    {
      //fprintf("U的第%d个元素是: %lf", i, nzval_U[i]);
      fprintf(fp_U, "U[%d]=%lf\n", i, *(nzval_U+i));
      /*for( j = 0; j < n; j++ )
      {
	fprintf(fp_U, "U[%d][%d]=%lf\n", i, j, nzval_U[i]);
      }*/
    }
    double *sol = (double*) ((DNformat*) B.Store)->nzval;
    for (i=0; i < Bstore->lda; i++)
    {
      fprintf(result, "x[%d]=%lf\n", i, sol[i]);
    }
      //printf("解元素为: %lf\n", sol[i]);
  
    if ( options.PrintStat ) StatPrint(&stat);
    StatFree(&stat);

    SUPERLU_FREE (rhs);
    SUPERLU_FREE (xact);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC("Exit main()");
#endif
}

