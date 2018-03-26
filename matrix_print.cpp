//the print the matrix
#include <iostream>
#include <cstdlib>
#include <complex>
#include<fstream>
#include<cstring>
#include "matrix.h"
#define pi 3.1415926535898 
using namespace std;
int main()
{
	int N;
	double mag;
	cout << "please input the N and mag:";
	cin >> N >> mag;
	getchar();
	ofstream out;

	out.open("matrix_new.txt", ios::out);
	double kx, ky;
	int i, j;
	int jp1, jm1;
	double n1, m1;
	double t1 = 1;
	double t2 = 0.75;
	for (int n = 0; n < 101; n++)
	{
		for (int m = 0; m < 101; m++)
		{
			kx = (-pi / double(N)) + 2.0*pi* double(n) / (double(N)*100.0);
			ky = -pi + 2.0*pi* double(m) / 100.0;
			complex<double> **a = new complex<double>*[2 * N];
			for (i = 0; i < 2 * N; i++)
			{
				a[i] = new complex<double>[2 * N];
			}
			createm(N, a, mag, kx, ky, t1, t2);
			for (i = 0; i < 2 * N; i++)
			{
				for (j = 0; j < 2 * N; j++)
				{
					out << a[i][j];
				}
				out << endl;

			}
			out << endl;//write the eigenvalue
			delete[]a;

		}
	}
	out.close();
}
