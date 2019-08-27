#include "Include/lu.h"
#include <stdio.h>
#include <stdlib.h>
# include "cs.h"
# include <cs_demo.h>

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
    Prob->b = cs_malloc (mn, sizeof (double)) ;
    Prob->x = cs_malloc (mn, sizeof (double)) ;
    Prob->resid = cs_malloc (mn, sizeof (double)) ;
    return ((!Prob->b || !Prob->x || !Prob->resid) ? free_problem (Prob) : Prob) ;
}

int main( int argc, char *argv[] )
{
    /* Read the file1 and store it by CSC: <a, asub, xa> */
    FILE *file1, *file2;
    double   *a;
    int      *asub, *xa;
    int mm, nn, nnz, nnn;
    double csparse_start, csparse_end, sparse_start, sparse_end; 
    nnn = nn + 1;
    file1 = fopen(argv[1], "r");
//     file2 = fopen(argv[2], "r");
//     dreadMM(file1, &mm, &nn, &nnz, &a, &asub, &xa);
    
//     printf("\n********************Matrix Information********************\n");
//     printf("Num of Rows: %d\n", mm);
//     printf("Num of Columns: %d\n", nn);
//     printf("Num of non-zeros: %d\n", nnz);
    
    /* CSparse LU decomposition*/
    problem *Prob = get_problem (file1, 1e-14) ;
    cs *A, *C ;
    double *b, *x, *resid,  t, tol ;
    int k, m, n, order, nb, ns, *r, *s, *rr, sprank ;
    csd *D ;
    
    /* CSparse: data pre-processing */
    if (!Prob) return (0) ;
    A = Prob->A ; C = Prob->C ; b = Prob->b ; x = Prob->x ; resid = Prob->resid;
    m = A->m ; n = A->n ;
    tol = Prob->sym ? 0.001 : 1 ;               
    D = cs_dmperm (C, 1) ;                      
    if (!D) return (0) ;
    nb = D->nb ; r = D->r ; s = D->s ; rr = D->rr ;
    sprank = rr [3] ;
    for (ns = 0, k = 0 ; k < nb ; k++)
    {
        ns += ((r [k+1] == r [k]+1) && (s [k+1] == s [k]+1)) ;
    }
    cs_dfree (D) ;
    if (m != n || sprank < n) return (1) ;     
    order = 1;
    rhs (x, b, m) ;                         
    
    /* CSparse: ordering and symbolic analysis */
    css *S ;
    csn *N ;
    if (!CS_CSC (A) || !b) return (0) ;    
    S = cs_sqr (order, C, 0) ;       
    
    /* CSparse: numeric LU factorization */
    csparse_start = microtime();
    N = cs_lu (C, S, tol) ;                
    csparse_end = microtime() - csparse_start;
   
    printf("\n********************CSparse LU Decomposition********************\n");
    printf("Num of non-zeros in L: %d\n", N->L->nzmax);
    printf("Num of non-zeros in U: %d\n", N->U->nzmax);
    printf("Time of CSparse LU Decomposition: %lf\n", csparse_end);
    
    /* MY LU Decomposition */
    printf("\n********************MY LU Decomposition: Check Results********************\n");
    sparse_start = microtime();
    lu_gp_sparse(C->x, C->i, C->p, n, N->L->nzmax, N->U->nzmax, S->q, N->pinv, N->L->i, N->L->p, N->U->i, N->U->p, N->L->x, N->U->x);
    sparse_end = microtime() - sparse_start;
    printf("Time of MY LU Decomposition: %lf\n", sparse_end);
    
    return 0;
}







