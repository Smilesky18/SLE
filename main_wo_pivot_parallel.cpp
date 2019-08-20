#include "lib/lu.h"
#include <stdio.h>
#include <stdlib.h>

int max ( int a, int b )
{
    if ( a > b ) return a;
    else return b;
}

int main( int argc, char *argv[] )
{
    FILE *fp;
    double   *a;
    int      *asub, *xa;
    int *index_row, *level, *depend_list;
    int m, n, nnz, sum = 0, num_row, num_row_selected;
    int i, j, k;
    fp = fopen(argv[1], "r");
    dreadMM(fp, &m, &n, &nnz, &a, &asub, &xa);
    level = ( int * ) malloc (sizeof( int *) * nnz);
    depend_list = ( int * ) malloc (sizeof( int *) * nnz);
    double start, finish;
    
    start = microtime();
    
    for ( i = 0; i < nnz; i++ )
    {
        level[i] = 0;
        depend_list[i] = 0;
    }
    
//     for ( i = 0; i < nnz; i++ )
//     {
//         printf("level[%d] = %d\n", i, level[i]);
// //         printf("depend_list[%d] = %d\n", i, depend_list[i]);
//     }

//     for ( i = 0; i < nnz; i++ )
//     {
//         printf("asub[%d] = %d\n", i, asub[i]);
//     }
//     for ( i = 0; i <=n; i++ )
//     {
//         printf("xa[%d] = %d\n", i, xa[i]);
//     }
    
    for ( i = 0; i < n; i++ )
    {
        num_row = xa[i+1] - xa[i];
        if ( num_row == 1 ) continue;
        if ( asub[xa[i]] == i )
        {
            for ( j = xa[i]+1; j < xa[i+1]; j++ )
            {
                level[j] = level[xa[i]]+1;
                depend_list[j] = sum;
                //union_list[sum] = xa[i];
                sum++;
            }
        }
        else
        {
            index_row = ( int *) malloc (sizeof(int *) * num_row);
            for ( j = 0; j < num_row; j++ )
            {
                index_row[j] = asub[xa[i]+j];
//                 if ( index_row[j] == 30 )
//                 {
//                     printf(" j = %d\n ", j);
//                 }
            }
            for ( j = 0; j < num_row; j++ )
            {
                if ( index_row[j] < i )
                {
                    for ( k = xa[index_row[j]]; k < xa[index_row[j]+1]; k++ )
                    {
                        if ( asub[k] <= index_row[j] ) continue;
                        for ( m = j+1; m < num_row; m++ )
                        {
                            if (asub[k] == index_row[m])
                            {
//                                 printf(" level[%d] = %d\n ", xa[i]+j, level[xa[i]+j]);
//                                 printf(" level[%d] = %d\n ", xa[k], level[xa[k]]);
                                level[xa[i]+m] = max(level[xa[i]+j], level[k]);
//                                 printf(" level[%d] = %d\n ", xa[i]+m, level[xa[i]+m]);
                                depend_list[xa[i]+m] = sum;
                                sum = sum+2;
                            }
                        }
                    }
                }
                else
                    break;
            }
        }

    }
    finish = microtime() - start;
    printf("The cost time is %lf\n", finish);
    return 0;
}
