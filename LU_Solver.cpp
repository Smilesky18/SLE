#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
//#include <random>
#include <boost/concept_check.hpp>
#include <boost/timer.hpp>
#include <time.h>
#include <fstream>
#include <iostream>

using namespace std;

void RandInit ( vector < vector <double> > & A, int & n)
{
    //cout << " 111 " << endl;
    srand((unsigned)time(NULL));
    for ( int i = 0; i < n; i++ )
      for ( int j = 0; j < n+1; j++ )
      {
	A[i][j] = rand();
      }
    //cout << " 222 " << endl;
}

void Init ( vector < vector <double> > & A, int & n)
{
    cout << " 请手动输入系数 " << endl;
    for ( int i = 0; i < n; i++ )
      for ( int j = 0; j < n+1; j++ )
      {
	cin >> A[i][j];
      }
}

bool Gauss(vector < vector <double> > & A, vector <double> & x, int & n)
{
  vector < vector<double> > Atemp(n);
  double l;
  for ( int k = 0; k < n; k++ )
    {
      Atemp[k].resize(n+1);
    }
  /*复制系数矩阵，防止修改原矩阵
  for ( int i =0; i < n; i++ )
  {
    vector <double> temp(n);
    for ( int j = 0; j< n+1; j++ )
    {
      temp[j] = A[i][j];
    }
    Atemp[i] = temp;
    temp.clear();
  }*/
  //复制系数矩阵，防止修改原矩阵
  for ( int i =0; i < n; i++ )
  {
    for ( int j = 0; j< n+1; j++ )
    {
      Atemp[i][j] = A[i][j];
    }
  }
  if ( Atemp[0][0]!= 0 )
  {
    //进行第一次消元
    for ( int i = 1; i < n; i++ )
    {
      l = Atemp[i][0]/Atemp[0][0];
      for ( int j = 0; j < n+1; j++ )
      {
	Atemp[i][j] = Atemp[i][j] - l*Atemp[0][j];
      }
    }
    //剩余的消元
    for ( int k = 1; k < n-1; k++ )
    {
      for ( int i =k+1; i < n; i++ )
      {
	if ( Atemp[k][k] == 0 )
	{
	  x.clear();
	  return false;
	}
	double l = Atemp[i][k] / Atemp[k][k];
	for ( int j = k; j < n+1; j++ )
	{
	  Atemp[i][j] = Atemp[i][j] - l*Atemp[k][j]; 
	}
      }
    }
    //回代
    x[n-1] = Atemp[n-1][n] / Atemp[n-1][n-1];
    for ( int k = n - 2; k >= 0; k-- )
    {
      double s = 0.0;
      for ( int j = k+1; j < n; j++ )
      {
	s += Atemp[k][j] * x[j];
      }
      x[k] = (Atemp[k][n] - s) / Atemp[k][k];
    }
  }
  vector < vector<double> >().swap(Atemp);
  return true;
}

bool LU(vector < vector <double> > & A, vector <double> & x, int n)
{
  vector < vector <double> > L(n), U(n);
  double sum_U = 0.0, sum_L = 0.0, sum_y = 0.0, sum_x = 0.0;
  vector < double > y(n);
  for ( int k = 0; k < n; k++ )
  {
      L[k].resize(n);
      U[k].resize(n);
  }
  //U的第一行的解向量
  for ( int i = 0; i < n; i++ )
  {
    U[0][i] = A[0][i];
  }
  //L的第一列的解向量
  for ( int i = 1; i < n; i++ )
  {
    L[i][0] = A[i][0]/U[0][0];
  }
  //对L的对角线元素赋值为1
  for ( int i = 0; i < n; i++ )
  {
    L[i][i] = 1;
  }
  //计算U的第r行，L的第r列元素（r=2,...n-1）
  for (int r = 1; r < n; r++ )
  {
    for ( int i = r; i < n; i++ )
    {
      for ( int k = 0; k < r; k++ )
      {
	sum_U += L[r][k]*U[k][i];
      }
      U[r][i] = A[r][i] - sum_U;
      sum_U = 0.0;
    }
    for ( int i = r+1; i < n; i++ )
    {
      for ( int k = 0; k < r; k++ )
      {
	sum_L += L[i][k]*U[k][r];
      }
      L[i][r] = ( A[i][r] - sum_L ) / U[r][r];
      sum_L = 0.0;
    }
  }
 /* cout << "打印L矩阵元素:" << " ";
  for ( int i = 0; i < n; i++ )
    for( int j = 0; j < n; j++ )
    {
      cout << L[i][j] << " ";
    }
  cout << endl << "打印U矩阵元素:" << " ";
  for ( int i = 0; i < n; i++ )
    for ( int j = 0; j < n; j++ )
    {
      cout << U[i][j] << " ";
    }
  cout << endl;*/
  //求解Ly=b, Ux=y的计算公式
  y[0] = A[0][n];  
  for ( int i = 1; i < n; i++ )
  {
    for ( int k = 0; k < i; k++ )
    {
      sum_y += L[i][k]*y[k]; 
    }
    y[i] = A[i][n] - sum_y;
    sum_y = 0.0;
  }
  x[n-1] = y[n-1] / U[n-1][n-1];
  for ( int i = n-2; i >= 0; i-- )
  {
    for ( int k = i + 1; k < n; k++ )
    {
      sum_x += U[i][k]*x[k];
    }
    x[i] = ( y[i] - sum_x ) / U[i][i];
    sum_x = 0.0;
  }
}

int main()
{
    int n;
    int select;
    cout << "请输入方程组的阶：" ;
    cin >> n;
    vector < vector<double> > A(n);
    cout << A.size() << endl;
    for ( int k = 0; k < A.size(); k++ )
    {
      A[k].resize(n+1);
    }
    vector <double> x;
    x.resize(n);
    if ( n <=1 )
    {
        cout << "方程组的阶数必须大于1" << endl;
        exit(1);
    }
    cout << " 请选择是自动产生系数矩阵(1)还是手动输入(0)?" << endl;
    cin >> select;
    if ( select == 1 )
    {
      RandInit(A, n);
    }
    else
    {
      Init(A, n);
    }
    boost::timer time;
    Gauss(A, x, n);
    cout << "高斯消去法占用的时间为: " << time.elapsed() << endl;
    cout << "高斯消去法的解向量为: " << endl;
    for ( int i = 0; i < n; i++ )
    {
      cout << x[i] << "  "; 
    }
    cout << "清空解向量的元素!" << endl;
    for ( int i = 0; i < n; i++ )
    {
      x[i] = 0;
    }
    for ( int i = 0; i < n; i++ )
    {
      cout << x[i] << "  "; 
    }
   /* ofstream solver_result;
    solver_result.open( "result.txt", ios::app );
    for ( int i = 0; i < n; i++ )
    {
      solver_result << x[i] << endl;
    } 
    solver_result.close();
    x.clear();*/
    boost::timer time_LU;
    LU(A, x, n);
    cout << " LU分解占用的时间为: " << time_LU.elapsed() << endl;
    cout << "LU分解法的解向量为: " << endl;
    for ( int i = 0; i < n; i++ )
    {
      cout << x[i] << "  "; 
    }
    return 0;
}
