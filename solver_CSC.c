# include "slu_ddefs.h"
# include <stdio.h>
#include <time.h>

double Abs(double x)
{
  return x < 0 ? -x : x;
}
typedef struct ColPointer {
  int row_index;
  int nonz_num;
  double nonz_value;
  struct ColPointer *next;
} LinkList;

void LU(double *a, int *asub, int *xa, double *x, int n)
{
  double **L, **U;
  double sum_U = 0.0, sum_L = 0.0, sum_y = 0.0, sum_x = 0.0;
  double *y;
  double *check_sum, check_max, sum = 0.0;
  int i, j, r, k, col_change_order;
  int *change_order;
  double temp = 0.0;
  int record_order_temp;
  int col_num;
  FILE *fp_L = fopen("L_CSC.txt", "w");
  FILE *fp_U = fopen("U_CSC.txt", "w");
  FILE *result = fopen("result_CSC.txt", "w");
  clock_t lu_start, lu_finish, trsv_start, trsv_finish;
  double lu_duration, trsv_duration;
  L = ( double ** )malloc(sizeof(double *) * n);
  U = ( double ** )malloc(sizeof(double *) * n);
  int sum_nonz_L = 0, sum_nonz_U = 0;
  double sum_element = n * n;
  double L_ratio, U_ratio;
  LinkList col_list[n];
  for ( i = 0; i < n; i++ )
  {
    L[i] = ( double * )malloc(sizeof(double) * n);
    U[i] = ( double * )malloc(sizeof(double) * n);
  }
  for ( i = 0; i < n; i++ )
  {
    L[i][i] = 1.0;
  }
  // from the right-looking to create the linklist
  for ( i = 0; i < n; i++ )
  {
      j = xa[i+1] - xa[i];
      LinkList *head, *node, *end;//head node, normal node, end node
      head = (LinkList*)malloc(sizeof(LinkList));
      end = head; //if the link is null, then end = head
      for ( r = 0; r < j; r++ )
      {
	node = (LinkList*)malloc(sizeof(LinkList));
	node->row_index = asub[xa[i]+r];
	node->nonz_value = a[xa[i]+r];
	node->nonz_num = j - r;
	end->next = node;
	end = node;
      }
      end->next = NULL;
      col_list[i] = *head;
  }
  /*for ( i = 0; i < n; i++ ) 
  {
    printf("the first value of col_list is: %d\n", col_list[i].next->next->nonz_num);
    printf("the second value of col_list is: %f\n", col_list[i].next->next->nonz_value);
    printf("the third value of col_list is: %d\n", col_list[i].next->next->row_index);
  }*/
  y = ( double * )malloc(sizeof(double) * n);
  //sum = ( double * )malloc(sizeof(double) * n);
  check_sum = ( double * )malloc(sizeof(double) * n);
  change_order = ( int * )malloc(sizeof(int) * n);
  for ( i = 0; i < n; i++ )
  {
    change_order[i] = i;
  }
  LinkList *col_0 = col_list[0].next;
  col_num = col_0->nonz_num;
  temp = col_0->nonz_value;
  for ( i = 1; i < col_num; i++ )
  {
    col_0 = col_0->next;
    if ( Abs(col_0->nonz_value) > Abs(temp) )
    {
      temp = col_0->nonz_value;
      change_order[0] = i;
      change_order[i] = 0;
    }
  } 
  U[0][0] = temp;
  // compute the first column of L
  col_0 = col_list[0].next;
  for ( i = 0; i < col_num; i++ )
  {
    /*if ( change_order[i] < i )
    {
      col_0 = col_0->next;
      continue;
    }*/                                                                                                                                                                                                                                                                                                                                                                                    
    j = change_order[col_0->row_index];
    L[j][0] = col_0->nonz_value / temp;
    col_0 = col_0->next;
  }
  // compute the first row of U
  for ( i = 1; i < n; i++ )
  {
    col_0 = col_list[i].next;
    col_num = col_0->nonz_num;
    for ( j = 0; j < col_num; j++ )
    {
      if ( change_order[col_0->row_index] == 0 )
      {
	U[0][i] = col_0->nonz_value;
	break;
      }
      col_0 = col_0->next;
    }
  }
  // compute the second... column/row of the L and the U
  for ( r = 1; r < n; r++ ) 
  {
    LinkList *col = col_list[r].next;                                                                    
    col_num = col->nonz_num;
    for ( i = r; i < n; i++ ) 
    {
      for ( k = 0; k < r; k++ )
      {
	sum += L[i][k]*U[k][r];
      }
      /*if ( col != NULL )
      {
	col_change_order = change_order[col->row_index];
	check_sum[col_change_order] = col->nonz_value - sum;
      }
      else
      {
	check_sum[i] = sum;
      }*/
      check_sum[i] = 0.0 - sum; 
      sum = 0.0;        
      printf("the first loop: check_sum[%d] = %f\n", i, check_sum[i]);
      //col = col->next;
      //check_sum[i] = A[i][r] - sum;
      //A[i][r] = A[i][r] - sum; 
    } 
    //printf(" col->nonz_num = %d\n ", col->nonz_num);
    //printf(" col->row_index = %d\n ", col->next->row_index);
    //printf(" col->nonz_value = %f\n", col->next->nonz_value);
    for ( j = 0; j < col_num; j++ )
    {
      //printf(" testing!! \n");
      //printf(" col->row_index = %d\n ", col->row_index);
      //printf(" col->nonz_value = %f\n", col->nonz_value);
      check_sum[change_order[col->row_index]] = col->nonz_value + check_sum[change_order[col->row_index]];
      //printf(" check_sum[%d] = %f\n", change_order[col->row_index], check_sum[change_order[col->row_index]]);
      col = col->next;
    }
    for ( i = 0; i < n; i++ )
    {
      printf(" check_sum[%d] = %f\n", i, check_sum[i]);
    }
    check_max = check_sum[r];
    for ( k = r+1; k < n; k++ )
    {
      if ( Abs(check_sum[k]) > Abs(check_max) )
      {
	printf(" abs compare \n");
	check_max = check_sum[k];
	change_order[r] = k;
	change_order[k] = r;
      }
    }
    if ( change_order[r] > r )
    {
      printf(" row interchanges!! ");
      record_order_temp = change_order[r];
      temp = check_sum[r];
      check_sum[r] = check_sum[record_order_temp];
      check_sum[record_order_temp] = temp;
      for ( i = 0; i < r; i++ )
      {
	temp = L[r][i];
	L[r][i] = L[record_order_temp][i]; 
	L[record_order_temp][i] = temp;
      }
    }
    // compute the value of r row and r column
    U[r][r] = check_max;
    printf("check_max = %f\n", check_max);
    //printf(" k = %d\n", k);
    
    for ( i = r+1; i < n; i++ ) 
    {
      LinkList *col = col_list[i].next;
      col_num = col->nonz_num;
      L[i][r] = check_sum[i] / U[r][r];
      printf(" L[%d][%d] = %f\n", i, r, L[i][r]);
      for ( k = 0; k < r; k++ )
      { 
	sum_U += L[r][k]*U[k][i];
	printf(" sum_U = %f\n", sum_U);
      }
      for ( j = 0; j < col_num; j++ )
      {
	if ( change_order[col->row_index] == r )
	{
	  U[r][i] = col->nonz_value - sum_U;
	  break;
	}
	col = col->next;
      }
      if ( U[r][i] == 0 )
      {
	U[r][i] = -sum_U;
      }
      printf(" U[%d][%d] = %f\n", r, i, U[r][i]);
      sum_U = 0.0;
    }  
    
  }
  lu_finish = clock();
  lu_duration = (double)(lu_finish - lu_start) / CLOCKS_PER_SEC;  
  printf( "The cost time of LU decomposition is: %f seconds\n", lu_duration );  
  //printf("The element of L: ");
  printf("the value of sum_nonz_L is: %d\n", sum_nonz_L);
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
  printf(" the value of sum_nonz_U is: %d\n", sum_nonz_U);
  for ( i = 0; i < n; i++ )
  {
    for( j = 0; j < n; j++ )
    {
      //printf("U[%d][%d]=%lf\n", i, j, U[i][j]);
      fprintf(fp_U, "U[%d][%d]=%lf\n", i, j, U[i][j]);
      if ( U[i][j] != 0 )
	sum_nonz_U++;
    }
  }
  U_ratio = (1 - sum_nonz_U/sum_element)*100;
  printf("sum_nonz_U = %d\n", sum_nonz_U);
  printf("The sparsity of U is: %lf%%\n", U_ratio);
  //solve Ly=b && Ux=y
  trsv_start = clock();
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
  trsv_finish = clock();
  trsv_duration = (double)(trsv_finish - trsv_start) / CLOCKS_PER_SEC;  
  printf( "The cost time of TRSV is: %f seconds\n", trsv_duration );
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
  free(change_order);
}*/

void LU(double *a, int *asub, int *xa, double *x, int n)
{
  double **L, **U;
  double sum_U = 0.0, sum_L = 0.0, sum_y = 0.0, sum_x = 0.0;
  double *y;
  double *check_sum, check_max, sum = 0.0;
  int i, j, r, k, col_change_order;
  int *change_order;
  double temp = 0.0;
  int record_order_temp;
  int col_num;
  FILE *fp_L = fopen("L_CSC.txt", "w");
  FILE *fp_U = fopen("U_CSC.txt", "w");
  FILE *result = fopen("result_CSC.txt", "w");
  clock_t lu_start, lu_finish, trsv_start, trsv_finish;
  double lu_duration, trsv_duration;
  L = ( double ** )malloc(sizeof(double *) * n);
  U = ( double ** )malloc(sizeof(double *) * n);
  int sum_nonz_L = 0, sum_nonz_U = 0;
  double sum_element = n * n;
  double L_ratio, U_ratio;
  LinkList col_list[n];
  for ( i = 0; i < n; i++ )
  {
    L[i] = ( double * )malloc(sizeof(double) * n);
    U[i] = ( double * )malloc(sizeof(double) * n);
  }
  for ( i = 0; i < n; i++ )
  {
    L[i][i] = 1.0;
  }
  // from the right-looking to create the linklist
  for ( i = 0; i < n; i++ )
  {
      j = xa[i+1] - xa[i];
      LinkList *head, *node, *end;//head node, normal node, end node
      head = (LinkList*)malloc(sizeof(LinkList));
      end = head; //if the link is null, then end = head
      for ( r = 0; r < j; r++ )
      {
	node = (LinkList*)malloc(sizeof(LinkList));
	node->row_index = asub[xa[i]+r];
	node->nonz_value = a[xa[i]+r];
	node->nonz_num = j - r;
	end->next = node;
	end = node;
      }
      end->next = NULL;
      col_list[i] = *head;
      /*printf("the first value of col_list is: %d\n", col_list[i].next->nonz_num);
      printf("the second value of col_list is: %f\n", col_list[i].next->nonz_value);
      printf("the third value of col_list is: %d\n", col_list[i].next->row_index);*/
  }
  /*for ( i = 0; i < n; i++ ) 
  {
    printf("the first value of col_list is: %d\n", col_list[i].next->next->nonz_num);
    printf("the second value of col_list is: %f\n", col_list[i].next->next->nonz_value);
    printf("the third value of col_list is: %d\n", col_list[i].next->next->row_index);
  }*/
  y = ( double * )malloc(sizeof(double) * n);
  //sum = ( double * )malloc(sizeof(double) * n);
  check_sum = ( double * )malloc(sizeof(double) * n);
  change_order = ( int * )malloc(sizeof(int) * n);
  for ( i = 0; i < n; i++ )
  {
    change_order[i] = i;
  }
  LinkList *col_0 = col_list[0].next;
  col_num = col_0->nonz_num;
  temp = col_0->nonz_value;
  for ( i = 1; i < col_num; i++ )
  {
    col_0 = col_0->next;
    if ( Abs(col_0->nonz_value) > Abs(temp) )
    {
      temp = col_0->nonz_value;
      change_order[0] = i;
      change_order[i] = 0;
    }
  } 
  U[0][0] = temp;
  // compute the first column of L
  col_0 = col_list[0].next;
  for ( i = 0; i < col_num; i++ )
  {
    /*if ( change_order[i] < i )
    {
      col_0 = col_0->next;
      continue;
    }*/                                                                                                                                                                                                                                                                                                                                                                                    
    j = change_order[col_0->row_index];
    L[j][0] = col_0->nonz_value / temp;
    col_0 = col_0->next;
  }
  // compute the first row of U
  for ( i = 1; i < n; i++ )
  {
    col_0 = col_list[i].next;
    col_num = col_0->nonz_num;
    for ( j = 0; j < col_num; j++ )
    {
      if ( change_order[col_0->row_index] == 0 )
      {
	U[0][i] = col_0->nonz_value;
	break;
      }
      col_0 = col_0->next;
    }
  }
  // compute the second... column/row of the L and the U
  for ( r = 1; r < n; r++ ) 
  {
    LinkList *col = col_list[r].next;                                                                    
    col_num = col->nonz_num;
    for ( i = r; i < n; i++ ) 
    {
      for ( k = 0; k < r; k++ )
      {
	sum += L[i][k]*U[k][r];
      }
      /*if ( col != NULL )
      {
	col_change_order = change_order[col->row_index];
	check_sum[col_change_order] = col->nonz_value - sum;
      }
      else
      {
	check_sum[i] = sum;
      }*/
      check_sum[i] = 0.0 - sum; 
      sum = 0.0;        
      printf("the first loop: check_sum[%d] = %f\n", i, check_sum[i]);
      //col = col->next;
      //check_sum[i] = A[i][r] - sum;
      //A[i][r] = A[i][r] - sum; 
    } 
    //printf(" col->nonz_num = %d\n ", col->nonz_num);
    //printf(" col->row_index = %d\n ", col->next->row_index);
    //printf(" col->nonz_value = %f\n", col->next->nonz_value);
    for ( j = 0; j < col_num; j++ )
    {
      //printf(" testing!! \n");
      //printf(" col->row_index = %d\n ", col->row_index);
      //printf(" col->nonz_value = %f\n", col->nonz_value);
      check_sum[change_order[col->row_index]] = col->nonz_value + check_sum[change_order[col->row_index]];
      //printf(" check_sum[%d] = %f\n", change_order[col->row_index], check_sum[change_order[col->row_index]]);
      col = col->next;
    }
    for ( i = 0; i < n; i++ )
    {
      printf(" check_sum[%d] = %f\n", i, check_sum[i]);
    }
    check_max = check_sum[r];
    for ( k = r+1; k < n; k++ )
    {
      if ( Abs(check_sum[k]) > Abs(check_max) )
      {
	printf(" abs compare \n");
	check_max = check_sum[k];
	change_order[r] = k;
	change_order[k] = r;
      }
    }
    if ( change_order[r] > r )
    {
      printf(" row interchanges!! ");
      record_order_temp = change_order[r];
      temp = check_sum[r];
      check_sum[r] = check_sum[record_order_temp];
      check_sum[record_order_temp] = temp;
      for ( i = 0; i < r; i++ )
      {
	temp = L[r][i];
	L[r][i] = L[record_order_temp][i]; 
	L[record_order_temp][i] = temp;
      }
    }
    // compute the value of r row and r column
    U[r][r] = check_max;
    printf("check_max = %f\n", check_max);
    //printf(" k = %d\n", k);
    
    for ( i = r+1; i < n; i++ ) 
    {
      LinkList *col = col_list[i].next;
      col_num = col->nonz_num;
      L[i][r] = check_sum[i] / U[r][r];
      printf(" L[%d][%d] = %f\n", i, r, L[i][r]);
      for ( k = 0; k < r; k++ )
      { 
	sum_U += L[r][k]*U[k][i];
	printf(" sum_U = %f\n", sum_U);
      }
      for ( j = 0; j < col_num; j++ )
      {
	if ( change_order[col->row_index] == r )
	{
	  U[r][i] = col->nonz_value - sum_U;
	  break;
	}
	col = col->next;
      }
      if ( U[r][i] == 0 )
      {
	U[r][i] = -sum_U;
      }
      printf(" U[%d][%d] = %f\n", r, i, U[r][i]);
      sum_U = 0.0;
    }  
    
  }
  lu_finish = clock();
  lu_duration = (double)(lu_finish - lu_start) / CLOCKS_PER_SEC;  
  printf( "The cost time of LU decomposition is: %f seconds\n", lu_duration );  
  //printf("The element of L: ");
  printf("the value of sum_nonz_L is: %d\n", sum_nonz_L);
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
  printf(" the value of sum_nonz_U is: %d\n", sum_nonz_U);
  for ( i = 0; i < n; i++ )
  {
    for( j = 0; j < n; j++ )
    {
      //printf("U[%d][%d]=%lf\n", i, j, U[i][j]);
      fprintf(fp_U, "U[%d][%d]=%lf\n", i, j, U[i][j]);
      if ( U[i][j] != 0 )
	sum_nonz_U++;
    }
  }
  U_ratio = (1 - sum_nonz_U/sum_element)*100;
  printf("sum_nonz_U = %d\n", sum_nonz_U);
  printf("The sparsity of U is: %lf%%\n", U_ratio);
  //solve Ly=b && Ux=y
  trsv_start = clock();
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
  trsv_finish = clock();
  trsv_duration = (double)(trsv_finish - trsv_start) / CLOCKS_PER_SEC;  
  printf( "The cost time of TRSV is: %f seconds\n", trsv_duration );
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
  free(change_order);
}

int main( int argc, char *argv[] )
{
    FILE *fp = stdin;
    double   *a;
    int      *asub, *xa;
    int m, n, nnz;
    int i, j;
    int NUMCol, counter_num = 0;
    double *x;
    clock_t start, finish, sum_start, sum_finish;
    double duration, sum_duration;
    sum_start = clock();
    printf(" This is only for testing! ");
    dreadMM(fp, &m, &n, &nnz, &a, &asub, &xa);
    printf(" This is for testing!!! ");
    /*for ( i = 0; i < nnz; i++ )
    {
      printf("a[%d]=%f\n", i, a[i]);
      printf("asub[%d]=%d\n", i, asub[i]);
    }
    for ( i = 0; i < (n+1); i++ )
    {
      printf("xa[%d]=%d\n", i, xa[i]);
    }*/
    x = ( double * )malloc(sizeof(double) * n);
    printf("\n%d %d %d\n", m, n, nnz);
    start = clock();
    LU( a, asub, xa, x, n);
    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;  
    printf( "The cost time of LU function is: %f seconds\n", duration );  
    printf("Normal end of execution");
    free(x);
    sum_finish = clock();
    sum_duration = (double)(sum_finish - sum_start) / CLOCKS_PER_SEC;  
    printf( "The cost time of the program is: %f seconds\n", sum_duration);  
    return 0;
}

