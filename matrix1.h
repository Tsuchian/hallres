//creatematrix.h
//hpp file to create Hamiltonian under strong magnetic field for Wei's boson metal model, written by Laurence
//If you have any questiones, please inform me to fix it.
//
#ifndef __CREATEMATRIX_H__
#define __CREATEMATRIX_H__
#include<iostream>
#include<complex>
#include<cmath>
#include<cstring>
#define pi 3.1415926535898 
using namespace std;
void createm(int N, complex <double> **a, double mag, double kx, double ky, double t1, double t2)
{
	//std::complex<double> **a = new std::complex<double>* [2 * N];
	int i, j, im1, jm1, jp1, ip1;
	double theta, theta1, theta2, theta3;

	//for orbital 1 to orbital 1
	//[0][0]-->[N][N]
	//
	for (i = 0; i < N; i++)
	{
		j = i;
		jm1 = j - 1;
		jp1 = j + 1;
		double i1 = i;
		if (i == N - 1)
		{
			jp1 = 0;
		}

		if (i == 0)
		{
			jm1 = N - 1;
		}
		a[i][jp1] = { t2*cos(kx),-sin(kx)*t2 };
		a[i][jm1] = { t2*cos(kx),t2*sin(kx) };

	}

	// for orbital 1 to orbital 2
	//
	//[0][N]-->[N][2N]
	//
	for (i = 0; i < N; i++)
	{
		j = i + N;
		jm1 = j - 1;
		jp1 = j + 1;
		double i1 = i;
		if (i == N - 1)
		{
			jp1 = N;
		}
		theta1 = mag*2.0*pi*(0.125 + i1 / 2.0) / double(N);
		a[i][j] = { t1*(cos(theta1) + cos(ky + theta1)),t1*(sin(ky + theta1) - sin(theta1)) };
		theta3 = mag*2.0*pi*(0.375 + i1 / 2.0) / double(N);
		a[i][jp1] = { t1*(cos(kx + theta3) + cos(ky - kx + theta3)),t1*(sin(ky - kx + theta3) - sin(kx + theta3)) };
	}
	//for orbital 2 to orbital 1
	//[N][0]-->[2N][N]
	//
	for (i = N; i < 2 * N; i++)
	{
		j = i - N;
		jm1 = j - 1;
		jp1 = j + 1;
		if (i == N)
		{
			jm1 = N - 1;

		}
		double i1 = j;
		theta1 = mag*2.0*pi*(0.125 + i1 / 2.0) / double(N);
		a[i][j] = { t1*(cos(theta1) + cos(ky + theta1)),t1*(sin(theta1) - sin(ky + theta1)) };
		theta2 = mag*2.0*pi*(i1 / 2.0 - 0.125) / double(N);
		a[i][jm1] = { t1*(cos(kx + theta2) + cos(ky - kx + theta2)), t1*(sin(kx + theta2) - sin(ky - kx + theta2)) };

	}
	//for orbital 2 to orbital 2
	//
	//[N][N]-->[2N][2N]
	//
	for (i = N; i < 2 * N; i++)
	{
		j = i;
		double i1 = i - N;
		a[i][j] = { t2*2.0*cos(ky - mag*2.0*pi*i1 / double(N)),0.0 };

	}
}
#endif
