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
#include"matrix_p_der.h"
#include"OptCond.h"
#include"dos.h"
#define pi 3.14159265358979323846
using namespace std;
int main()
{
	int N = 1;
	double mag = 1.0;
	FILE *out, *matrix, *pm;
	out = fopen("sigma.txt", "w");
	double kx, ky;
	int i, j;
	double kT, mu, omega;
	double temp[2] = { 0 };
	double sigmaxy, sigmaxx; //optical conductivity
	double t1 = 33.5;
	double t2 = 25.8;
	int t = 0;
	kT = 5.00;
	mu = dos(N, mag, 0.075, kT, t1, t2,0.0001);
	//cout << mu;
    //getchar();
    for (omega = 0.1; omega < 300; omega = omega + 0.1)
	{
		sigmaxx = sigmaxy = 0.0;
		for (int n = 0; n < 201; n++)
		{
			for (int m = 0; m < 201; m++)
			{
				//kx = (-pi / double(N)) + 2.0*pi* double(n) / (double(N)*100.0);
				//initial the Brillouin zone
				temp[0]=temp[1]=0;
				kx = -pi + 2.0*pi* double(n) / 100.0;
				ky = -pi + 2.0*pi* double(m) / 100.0;
				complex<double> **a = new complex<double>*[4 * N];//initial the matrix
				complex<double> **dx = new complex<double>*[4 * N];
				complex<double> **dy = new complex<double>*[4 * N];
				for (i = 0; i < 4 * N; i++)
				{
					a[i] = new complex<double>[4 * N];
					dx[i] = new complex<double>[4 * N];
					dy[i] = new complex<double>[4 * N];
				}
				createm(N, a, dx, dy, mag, kx, ky, t1, t2);//create the Hamiltonian and its partial derivative
				complex<double> *b = new complex<double>[4 * N * 4 * N];
				for (int k = 0; k < 4 * N; k++)
				{
					for (int l = 0; l < 4 * N; l++)
					{
						b[k * 4 * N + l] = a[k][l];
					}
				}
				char jobz = 'V';
				char uplo = 'U';
				double *w = new double[4 * N];
				lapack_int lda = 4 * N;
				LAPACKE_zheev(LAPACK_ROW_MAJOR, jobz, uplo, 4 * N, b, lda, w); //diagonalization
				Optcond(4 * N, kT, mu, kx, ky, omega, w, dx, dy, b, temp);
				sigmaxx = sigmaxx + temp[0];
				sigmaxy = sigmaxy + temp[1];
				for (i=0;i<4 * N;i++)
				{
					delete[]a[i];
					delete[]dx[i];
					delete[]dy[i];
				}
				delete []w;
				delete []b;
			}
		}
        sigmaxx = sigmaxx/sqrt(201*201.0);
        sigmaxy = sigmaxy/sqrt(201*201.0);
		double Rxy = sigmaxy / pow(sigmaxx, 2);
		fprintf(out, "%f %f %f %f\n", omega, sigmaxx, sigmaxy, Rxy);
	}
	fclose(out);
}
