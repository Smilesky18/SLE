#include "Lib/lu.h"
#include <stdio.h>
#include <stdlib.h>
# include "cs.h"
# include <cs_demo.h>

// bool equal( double a, double b )
// {
//   if ( Abs(a-b) < 0.0001)
//   {
//     return true;
//   }
//   else
//   {
//     return false;
//   }
// }

void rhs (double *x, double *b, int m)
{
    int i ;
    for (i = 0 ; i < m ; i++) b [i] = 1 + ((double) i) / m ;
    for (i = 0 ; i < m ; i++) x [i] = b [i] ;
}

int is_sym (cs *A)
{
    int is_upper, is_lower, j, p, n = A->n, m = A->m, *Ap = A->p, *Ai = A->i ;
    if (m != n) return (0) ;
    is_upper = 1 ;
    is_lower = 1 ;
    for (j = 0 ; j < n ; j++)
    {
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            if (Ai [p] > j) is_upper = 0 ;
            if (Ai [p] < j) is_lower = 0 ;
        }
    }
    return (is_upper ? 1 : (is_lower ? -1 : 0)) ;
}
int dropdiag (int i, int j, double aij, void *other) { return (i != j) ;}

cs *make_sym (cs *A)
{
    cs *AT, *C ;
    AT = cs_transpose (A, 1) ;          /* AT = A' */
    cs_fkeep (AT, &dropdiag, NULL) ;    /* drop diagonal entries from AT */
    C = cs_add (A, AT, 1, 1) ;          /* C = A+AT */
    cs_spfree (AT) ;
    return (C) ;
}

problem *free_problem (problem *Prob)
{
    if (!Prob) return (NULL) ;
    cs_spfree (Prob->A) ;
    if (Prob->sym) cs_spfree (Prob->C) ;
    cs_free (Prob->b) ;
    cs_free (Prob->x) ;
    cs_free (Prob->resid) ;
    return (cs_free (Prob)) ;
}

problem *get_problem (FILE *f, double tol)
{
    cs *T, *A, *C ;
    int sym, m, n, mn, nz1, nz2 ;
    problem *Prob ;
    Prob = cs_calloc (1, sizeof (problem)) ;
    if (!Prob) return (NULL) ;
    T = cs_load (f) ;                   /* load triplet matrix T from a file */
    Prob->A = A = cs_compress (T) ;     /* A = compressed-column form of T */
    cs_spfree (T) ;                     /* clear T */
    if (!cs_dupl (A)) return (free_problem (Prob)) ; /* sum up duplicates */
    Prob->sym = sym = is_sym (A) ;      /* determine if A is symmetric */
    m = A->m ; n = A->n ;
    mn = CS_MAX (m,n) ;
    nz1 = A->p [n] ;
    cs_dropzeros (A) ;                  /* drop zero entries */
    nz2 = A->p [n] ;
    if (tol > 0) cs_droptol (A, tol) ;  /* drop tiny entries (just to test) */
    Prob->C = C = sym ? make_sym (A) : A ;  /* C = A + triu(A,1)', or C=A */
    if (!C) return (free_problem (Prob)) ;
    printf ("\n--- Matrix: %g-by-%g, nnz: %g (sym: %g: nnz %g), norm: %8.2e\n",
            (double) m, (double) n, (double) (A->p [n]), (double) sym,
            (double) (sym ? C->p [n] : 0), cs_norm (C)) ;
    if (nz1 != nz2) printf ("zero entries dropped: %g\n", (double) (nz1 - nz2));
    if (nz2 != A->p [n]) printf ("tiny entries dropped: %g\n",
            (double) (nz2 - A->p [n])) ;
    Prob->b = cs_malloc (mn, sizeof (double)) ;
    Prob->x = cs_malloc (mn, sizeof (double)) ;
    Prob->resid = cs_malloc (mn, sizeof (double)) ;
    return ((!Prob->b || !Prob->x || !Prob->resid) ? free_problem (Prob) : Prob) ;
}



int main( int argc, char *argv[] )
{
    FILE *fp, *fpp;
    double   *a, *xx;
    lu_data *lu;
    int      *asub, *xa;
    int mm, nn, nnz, nnn, temp;
    int i, j, superlu_nnz_U, sum_l = 0, sum_u = 0;
    int nzl, nzu, *perm_c, *perm_r, *iperm_r, *asub_L, *asub_U, *xa_L, *xa_U;
    double csparse_start, csparse_end, sparse_start, sparse_end; 
    fp = fopen(argv[1], "r");
    fpp = fopen(argv[2], "r");
//     printf("argv[1] = %s\n", argv[1]);
    dreadMM(fp, &mm, &nn, &nnz, &a, &asub, &xa);
    nnn = nn + 1;
    perm_c = (int *)malloc(sizeof(int) * nn );
    perm_r = (int *)malloc(sizeof(int) * nn );
    iperm_r = (int *)malloc(sizeof(int) * nn );
    xa_L = (int *)malloc(sizeof(int) * nnn );
    xa_U = (int *)malloc(sizeof(int) * nnn );
    csparse_start = microtime();
    problem *Prob = get_problem (fpp, 1e-14) ;
    cs *A, *C ;
    double *b, *x, *resid,  t, tol ;
    int k, m, n, order, nb, ns, *r, *s, *rr, sprank ;
    csd *D ;
    if (!Prob) return (0) ;
    A = Prob->A ; C = Prob->C ; b = Prob->b ; x = Prob->x ; resid = Prob->resid;
    m = A->m ; n = A->n ;
    tol = Prob->sym ? 0.001 : 1 ;               /* partial pivoting tolerance */
    D = cs_dmperm (C, 1) ;                      /* randomized dmperm analysis */
    if (!D) return (0) ;
    nb = D->nb ; r = D->r ; s = D->s ; rr = D->rr ;
    sprank = rr [3] ;
    for (ns = 0, k = 0 ; k < nb ; k++)
    {
        ns += ((r [k+1] == r [k]+1) && (s [k+1] == s [k]+1)) ;
    }
    printf ("blocks: %g singletons: %g structural rank: %g\n",
        (double) nb, (double) ns, (double) sprank) ;
    cs_dfree (D) ;
    if (m != n || sprank < n) return (1) ;      /* return if rect. or singular*/
    order = 1;
    rhs (x, b, m) ;                         /* compute right-hand side */
//     cs_lusol (order, C, x, tol) ;      /* solve Ax=b with LU */
//     printf("C.nzmax = %d C.m = %d\n", C->nzmax, C->m);
    
    css *S ;
    csn *N ;
    if (!CS_CSC (A) || !b) return (0) ;     /* check inputs */
//     csparse_start = microtime();
    S = cs_sqr (order, C, 0) ;              /* ordering and symbolic analysis */
    N = cs_lu (C, S, tol) ;                 /* numeric LU factorization */
    csparse_end = microtime() - csparse_start;
    int lnz = N->L->nzmax;
    int unz = N->U->nzmax;
    printf("# of L is %d  # of U is %d\n", lnz, unz);
    printf("# of L is %d  # of U is %d\n", N->L->nz, N->U->nz);
    printf(" the time of csparse is %lf\n", csparse_end);
    
//     printf("row permutation is: ");
//     for ( i = 0; i < n; i++ ) printf("%d ", N->pinv[i]);
//     printf("\ncolumn permutation is: ");
//     for ( i = 0; i < n; i++ ) printf("%d ", S->q[i]);
// //     printf("\nnnz of U is: %d\n", U->nz) ;
// //     printf("nnz of L is: %d\n", L->nz) ;
//     printf("\ninformation of L is: \n");
//     printf("\nrow index of L is: ");
//     for ( i = 0; i < lnz; i++ ) printf("%d ", N->L->i[i]);
//     printf("\noffset of column is: ");
//     for ( i = 0; i <= n; i++) printf("%d ", N->L->p[i]);
//     
//     printf("\ninformation of U is: \n");
//     printf("\nrow index of U is: ");
//     for ( i = 0; i < unz; i++ ) printf("%d ", N->U->i[i]);
//     printf("\noffset of column is: ");
//     for ( i = 0; i <= n; i++) printf("%d ", N->U->p[i]);
//     printf("\n");
    
    perm_c = S->q;
    iperm_r = N->pinv;
    asub_L = N->L->i;
    xa_L = N->L->p;
    asub_U = N->U->i;
    xa_U = N->U->p;
    nzl = lnz;
    nzu = unz;
    for ( i = 0; i < n; i++ )
    {
        temp = iperm_r[i];
        perm_r[temp] = i;
    }
//     printf("row permutation is: ");
//     for ( i = 0; i < n; i++ ) printf("%d ", iperm_r[i]);
//     printf("\ncolumn permutation is: ");
//     for ( i = 0; i < n; i++ ) printf("%d ", perm_c[i]);
//     printf("\nthe information of L is: \n");
//     printf("row indices of L is: ");
//     for ( i = 0; i < nzl; i++ ) printf("%d ", asub_L[i]);
//     printf("\ncolumn offset of L is: ");
//     for ( i = 0; i <= n; i++ ) printf("%d ", xa_L[i]);
//     printf("\nthe information of U is: \n");
//     printf("row indices of U is: ");
//     for ( i = 0; i < nzu; i++ ) printf("%d ", asub_U[i]);
//     printf("\ncolumn offset of U is: ");
//     for ( i = 0; i <= n; i++ ) printf("%d ", xa_U[i]);
    
    sparse_start = microtime();
    lu_gp_sparse(a, asub, xa, n, nzl, nzu, perm_c, iperm_r, asub_L, xa_L, asub_U, xa_U, N->L->x, N->U->x);
    sparse_end = microtime() - sparse_start;
    printf(" the time of sparse is %lf\n", sparse_end);
   
//     printf("N.L.nzmax = %d\n", N->L->nzmax);
    printf("\nNormal end of execution");
    return 0;
}






























































//     sparse_load(fpp, &nzl, &nzu, &perm_c, &perm_r, &asub_L, &asub_U, &xa_L, &xa_U);
//     printf("-----------------------super_LU----------------------\n");
//     start_lu = microtime();
//     x = ( double *)malloc(sizeof(double) * n );
//     x = super_lu( argv[1], x );
//     finish_lu = microtime() - start_lu;
//     for ( i = 0; i < n; i++ )
//     {
//       printf(" x[%d] = %lf\n", i, x[i]);
//     }
//     test = x + n;
//     for ( i = 0; i < 3*n; i++ )
//     {
//       printf(" test[%d] = %lf\n ", i, test[i]);
//     }
//     printf("The average cost time of super_LU function is: %f seconds\n", finish_lu);
//     printf("-----------------------LU_GP----------------------\n");
//     start_GP = microtime();
//     GP_x = lu_gp( a, asub, xa, n );
//     finish_GP = microtime() - start_GP;
//     pivot_ratio = ( GP_x[n] / finish_GP ) * 100;
//     printf(" the pivot time ration is: %lf\n ", pivot_ratio);
//     printf("The average cost time of LU_GP function is: %f seconds\n", finish_GP);
//     printf("-----------------------LU_GP_v2----------------------\n");
//     start_GP_v2 = microtime();
//     GP_x_v2 = lu_gp_v2( a, asub, xa, n );
//     finish_GP_v2 = microtime() - start_GP_v2;
//     printf("The average cost time of LU_GP_v2 function is: %f seconds\n", finish_GP_v2);
//     printf("-----------------------LU_GP_v3----------------------\n");
//     start_GP_v3 = microtime();
//     GP_x_v3 = lu_gp_v3( a, asub, xa, n );
//     finish_GP_v3 = microtime() - start_GP_v3;
//     for ( i = 0; i < n; i++ )
//     {
//       printf(" x[%d] = %lf\n", i, GP_x_v3[i]);
//     }
//     printf("The average cost time of LU_GP_v23function is: %f seconds\n", finish_GP_v3);
//     printf("-----------------------LU_SPARSE----------------------\n");
//     start_GP_amd = microtime();
//     GP_x_amd = lu_sparse( a, asub, xa, n );
//     finish_GP_amd = microtime() - start_GP_amd;
//     //pivot_ratio = ( GP_x[n] / finish_GP ) * 100;
//     //printf(" the pivot time ration is: %lf\n ", pivot_ratio);
//     printf("The average cost time of LU_SPARSE function is: %f seconds\n", finish_GP_amd);
        
    
//     lu_gp(a, asub, xa, n);
