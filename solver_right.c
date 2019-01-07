# include "slu_ddefs.h"
# include <stdio.h>

double Abs(double x)
{
  return x < 0 ? -x : x;
}
void LU(double A[][800], double x[800], int n)
{
  double L[n][n], U[n][n];
  double sum_U = 0.0, sum_L = 0.0, sum_y = 0.0, sum_x = 0.0;
  double y[n];
  double check_sum[n], sum = 0.0, check_max;
  int i, j, r, k, check;
  int record_order[n];
  double temp = 0.0;
  double max = 0.0;
  int record_order_temp;
  FILE *fp_L = fopen("L_pviot.txt", "w");
  FILE *fp_U = fopen("U_pviot.txt", "w");
  FILE *result = fopen("result.txt", "w");
  for ( i = 0; i < n; i++ )
  {
    record_order[i] = i;
  }
  temp = A[0][0];
  for ( i = 1; i < n; i++ )
  {
    if ( Abs(temp) < Abs(A[i][0]) )
    {
      record_order[0] = i;
      temp = A[i][0];
    }
  }
  if ( record_order[0] != 0 )
  {
    record_order_temp = record_order[0];
    for ( i = 0; i < n; i++ )
    {
      temp = A[0][i];
      A[0][i] = A[record_order_temp][i];
      A[record_order_temp][i] = temp;
      //printf("the %dth time to change the pviot!", record_order_temp);
    }
    temp = 0.0;
  }
  // 将调换的次序记录在record_order中，然后调换矩阵中的数值
  //record_order[0] = record_order_temp;
  
  //U的第一行的解向量
  for ( i = 0; i < n; i++ )
  {
    U[0][i] = A[0][i];
  }
  //L的第一列的解向量
  for ( i = 1; i < n; i++ )
  {
    L[i][0] = A[i][0]/U[0][0];
  }
  //对L的对角线元素赋值为1
  for ( i = 0; i < n; i++ )
  {
    L[i][i] = 1;
  }
  // 从第1步开始，进行选主元的LU分解法
  for ( r = 1; r < n; r++ )
  {
    for ( i = r; i < n; i++ )
    {
      for ( k = 0; k < r; k++ )
      {
	sum += L[i][k]*U[k][r];
      }
      check_sum[i] = A[i][r] - sum;
      A[i][r] = A[i][r] - sum; 
      //printf("check_sum[%d] = %lf\n", r, check_sum[r]);
      sum = 0.0;
    }
    check_max = check_sum[r];
    for ( k = r+1; k < n; k++ )
    {
      if ( Abs(check_sum[k]) > Abs(check_max) )
      {
	check_max = check_sum[k];
	//printf("check_max = %lf\n", check_max);
	record_order[r] = k;
	//printf("record_order[%d] = %d\n", r, k);
      }
    }
    // 交换第r行和第record_order[r]行的数据
    if ( record_order[r] != r )
    {
      record_order_temp = record_order[r];
      for ( i = 0; i < n; i++ )
      {
	temp = A[r][i];
	A[r][i] = A[record_order_temp][i];
	A[record_order_temp][i] = temp;
      }
      for ( i = 0; i < r; i++ )
      {
	temp = L[r][i];
	L[r][i] = L[record_order_temp][i];
	L[record_order_temp][i] = temp;
      }
    }
    // 计算U的第r行元素和L的第r列元素
    U[r][r] = A[r][r];
    //                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    max = 0.0;
    for ( i = r+1; i < n; i++ )
      
    {
      L[i][r] = A[i][r] / U[r][r];
      for ( k = 0; k < r; k++ )
      { 
	sum_U += L[r][k]*U[k][i];
      }
      U[r][i] = A[r][i] - sum_U;
      sum_U = 0.0;
    }    
  }
  printf("The element of L: ");
  for ( i = 0; i < n; i++ )
  {
    for( j = 0; j < n; j++ )
    {
      //if ( (char)L[i][j] == nan)
      //{
	fprintf(fp_L, "L[%d][%d]=%lf\n", i, j, L[i][j]);
      //}
      //printf("L[%d][%d]=%f ", i, j, L[i][j]);
     // fprintf(fp_L, "%lf ", L[i][j]);
    }
   // fprintf(fp_L, "\n");
  }
  printf("The element of U: ");
  for ( i = 0; i < n; i++ )
  {
    for( j = 0; j < n; j++ )
    {
      //printf("U[%d][%d]=%f ", i, j, U[i][j]);
      //if ( (char)U[i][j] == nan)
      //{
	fprintf(fp_U, "U[%d][%d]=%lf\n", i, j, U[i][j]);
      //}
      //fprintf(fp_U, "%lf ", U[i][j]);
    }
   // fprintf(fp_U, "\n");
  }
  //求解Ly=b, Ux=y的计算公式
  y[0] = A[0][n];  
  for ( i = 1; i < n; i++ )
  {
    for ( k = 0; k < i; k++ )
    {
      sum_y += L[i][k]*y[k]; 
    }
    y[i] = A[i][n] - sum_y;
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
  for ( i = 0; i < n; i++ )
  {
    //printf("解元素为: %lf\n", x[i]);
    fprintf(result, "x[%d]=%lf\n", i, x[i]);
  }
}

int main( int argc, char *argv[] )
{
    FILE *fp = stdin;
    FILE *fp_A_ORI = fopen("A_ORI.txt", "w");
    FILE *fp_A_LU = fopen("A_LU.txt", "w");
    double   *a;
    int      *asub, *xa;
    int m, n, nnz;
    //const int numarray;
    int i, j;
    int NUMCol, counter_num = 0;
    double A[800][800];
    double x[800];
    dreadMM(fp, &m, &n, &nnz, &a, &asub, &xa);
    /*for ( i = 0; i < nnz; i++ )
    {
      printf("%d\n", asub[i]);
    }*/
    printf("\n%d %d %d", m, n, nnz);
    for ( j = 0; j < n; j++ )
    {
      NUMCol = xa[j+1] - xa[j];
      for ( i = 0; i < NUMCol; i++ )
      {
	A[asub[counter_num]][j] = a[counter_num];
	counter_num++;
      }
    }
    for ( i = 0; i < n; i++ )
    {
       for ( j = 0; j < n; j++ )
       {
	  fprintf(fp_A_ORI, "%lf ",A[i][j]);
       }
       fprintf(fp_A_ORI, "\n");
    }
    for ( i = 0; i < n; i++ )
    {
      A[i][n] = 1.0;
    }
    /*for ( i = 0; i < n; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
	printf("A[%d][%d]=%f ", i, j, A[i][j]);
      }
    }*/
    LU( A, x, n);
    /*for ( i = 0; i < n; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
	printf("A[%d][%d]=%f ", i, j, A[i][j]);
      }
    }*/
    for ( i = 0; i < n; i++ )
    {
       for ( j = 0; j < n; j++ )
       {
	  fprintf(fp_A_LU, "%lf	",A[i][j]);
       }
       fprintf(fp_A_LU, "\n");
    }
     
    /*for ( i = 0; i < n; i++ )
    {
      printf("在main函数中的解元素为: %lf\n", x[i]);
    }*/	
    return 0;
}

