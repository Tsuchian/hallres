#ifndef lapack_complex_double 
#define lapack_complex_double MKL_Complex16
#endif
#include <iostream>
#include <cstdlib>
#include <complex>
#define MKL_Complex16 std::complex<double>
#include "mkl.h"
#include<fstream>
#include<cstring>
#include "creatematrix.h"
#define pi 3.1415926535898 
using namespace std;
int main()
{
	int N = 1;
	double mag = 0.0;
	FILE *out, *matrix, *pm;
	out = fopen("eigenvalue.txt", "w");
	double kx, ky;
	int i, j;
	int jp1, jm1;
	double n1, m1;
	double t1 = 20.0;
	double t2 = 10.0;
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
			complex<double> *b = new complex<double>[2 * N * 2 * N];
			for (int k = 0; k<2 * N; k++)
			{
				for (int l = 0; l<4 * N; l++)
				{
					b[k * 2 * N + l] = a[k][l];
				}
			}
			char jobz = 'V';
			char uplo = 'U';
			double w[2 * N];
			lapack_int lda = 2 * N;
			LAPACKE_zheev(LAPACK_ROW_MAJOR, jobz, uplo, 2 * N, b, lda, w); //diagonalization
			fprintf(out, "%f %f ", kx, ky);
			for (i = 0; i < 2 * N; i++)
			{
				//fprintf(out, "%f %f ",kx,ky);
				fprintf(out, " %f ", w[i]);
			}
			fprintf(out, "\n");//write the eigenvalue
			delete[]a;

		}
	}
}
