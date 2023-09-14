#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

//polinom Lagranja
double
Polinom_Lagranja (double *x, double *y, int N, double x0)
{
  double result = 0.0;
  double P = 1.0;
  for (int i = 0; i < N; i++)
    {
      P = 1.0;
      for (int j = 0; j < N; j++)
	if (i != j)
	  P = P * (x0 - x[j]) / (x[i] - x[j]);
      result = result + P * y[i];
    }
  return result;
}

// Variante 19
double
myFunction (double x)
{
  return (x * x + 2) + sin (sqrt (x * x + 2));
}

int
main ()
{
  int N;
  double x1 = -1;
  double x2 = 1;

  cout << "Enter the depth of the partition N(N > 1) :";
  cin >> N;

  int M = 3 * N;		// smaller partitioning

  double step = (x2 - x1) / N;	    // step for points
  double *xn = new double[N + 1];	//array of x values
  double *yn = new double[N + 1];	//array of y values
  double *xm = new double[M + 1];	//array of x values
  double *ym = new double[M + 1];	//array of y values

  // filling the array
  for (int i = 0; i <= N; i++)
    {
      xn[i] = x1 + i * step;
      yn[i] = myFunction (xn[i]);
    }
  for (int i = 0; i <= M; i++)
    {
      xm[i] = x1 + i * step;
      ym[i] = myFunction (xm[i]);
    }

  double xmm = x1;
  double xnn = x1;
  cout << endl << "Lagrangian interpolation" << endl;
  cout << endl
    << setw (5) << "   x\t\t"
    << setw (5) << "   F(x)\t\t"
    << setw (5) << "   L(x)\t\t" << setw (5) << "   Delta\t\t" << endl;

  //calculation and output of the result
  for (int i = 0; i <= N; i++)
    {
      cout << fixed << setprecision (5)
	<< setw (10) << xnn << " \t "
	<< setw (10) << myFunction (xnn) << " \t "
	<< setw (10) << Polinom_Lagranja (xn, yn, N, xnn) << " \t "
	<< setw (10) << abs (myFunction (xnn) -
			     Polinom_Lagranja (xn, yn, N, xnn)) << endl;
      xnn = xnn + step;
    }

  step = (x2 - x1) / M;		// changing the step
  cout << endl;
  for (int i = 0; i <= M; i++)
    {
      cout << fixed << setprecision (5)
	<< setw (10) << xmm << " \t "
	<< setw (10) << myFunction (xmm) << " \t "
	<< setw (10) << Polinom_Lagranja (xm, ym, M, xmm) << " \t "
	<< setw (10) << abs (myFunction (xmm) -
			     Polinom_Lagranja (xm, ym, M, xmm)) << endl;
      xmm = xmm + step;
    }
  return 0;
}
