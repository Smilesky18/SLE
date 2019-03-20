#include "lu.h"
#include <stdio.h>

void super_lu( int m, int n, int nnz, int *xa, int *asub, double *a )
{
    SuperMatrix A;
    NCformat *Astore;
    int      *perm_c; /* column permutation vector */
    int      *perm_r; /* row permutations from partial pivoting */
    SuperMatrix L;      /* factor L */
    SCformat *Lstore;
    SuperMatrix U;      /* factor U */
    NCformat *Ustore;
    SuperMatrix B;
    int      nrhs, ldx, info;
    double   *xact, *rhs;
    mem_usage_t   mem_usage;
    superlu_options_t options;
    SuperLUStat_t stat;
    int i;
    FILE *result = fopen("result/SuperLU_solution.txt", "w");
    
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC("Enter super_lu()");
#endif
    set_default_options(&options);

    dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
    //Astore = A.Store;
    //printf("Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, Astore->nnz);
    
    nrhs   = 1;
    if ( !(rhs = doubleMalloc(m * nrhs)) ) ABORT("Malloc fails for rhs[].");
    for ( i = 0; i < m; i++ )
    {
      rhs[i] = 1.0;
    }
    dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);
    //xact = doubleMalloc(n * nrhs);
    //ldx = n;
    //dGenXtrue(n, nrhs, xact, ldx);
    //dFillRHS(options.Trans, nrhs, xact, ldx, &A, &B);

    if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
    if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");

    /* Initialize the statistics variables. */
    StatInit(&stat);
    dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
    
    if ( info == 0 ) {

	/* This is how you could access the solution matrix. */
        double *sol = (double*) ((DNformat*) B.Store)->nzval;
	for ( i = 0; i < ((DNformat*) B.Store)->lda; i++ )
	{
	  fprintf(result, "x[%d]=%lf\n", i, sol[i]);
	  //printf("The solution of the LS is %f\n", sol[i]);
	}
	 /* Compute the infinity norm of the error. */
	dinf_norm_error(nrhs, &B, xact);

	Lstore = (SCformat *) L.Store;
	Ustore = (NCformat *) U.Store;
	double *Lsol = (double*) (Lstore->nzval);
	double *Usol = (double*) (Ustore->nzval);
	/*for ( i = 0; i < Lstore->nnz; i++ )
	{
	  printf("The non-zero value of L is %f\n", Lsol[i]);
	}
	for ( i = 0; i < Ustore->nnz; i++ )
	{
	  printf("The non-zero value of U is %f\n", Usol[i]);
	}*/
    	/**printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
    	printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
    	printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - n);
    	printf("FILL ratio = %.1f\n", (float)(Lstore->nnz + Ustore->nnz - n)/nnz);*/
	
	dQuerySpace(&L, &U, &mem_usage);
	//printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
	//       mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
	
    } else {
	//printf("dgssv() error returns INFO= %d\n", info);
	if ( info <= n ) { /* factorization completes */
	    dQuerySpace(&L, &U, &mem_usage);
	    //printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
	//	   mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
	}
    }

    if ( options.PrintStat ) StatPrint(&stat);
    StatFree(&stat);

    SUPERLU_FREE (rhs);
    //SUPERLU_FREE (xact);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC("Exit super_lu()");
#endif
}

