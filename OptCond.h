//To calculate the optical conductivity
// Writen by Long and Laurence
//
//N		:	order of matrix
//kT		:	Energy scale of temperature
//mu		:	chemical potential
//E		:	Matrix store Eigenvalue
//dx		:	Partial H/ Partial kx(Partial H/Partial ky for dy)
#ifndef __OPTCOND_H__
#define __OPTCOND_H__
#include<iostream>
#include<complex>
#include<cmath>
using namespace std;

void Optcond(int N, double kT, double mu, double kx, double ky, double omega, double *E, \
	complex<double> **dx, complex<double> **dy, complex<double> *Eigenvector, double temp[2])
{
	complex<double> delta;
	delta = { 0,1.0 };
	double *Nb = new double[N];
	complex <double> vx, vy;//velocity at x or y direction <psi| partial H/partial kx |psi>
	for (int i = 0; i < N; i++)
	{
		Nb[i] = 1 / (exp((E[i] - mu) / kT) - 1);
		for (int j = i + 1; j < N; j++)
		{
			vx = vy = { 0.0,0.0 };
			for (int n = 0; n < N; n++)
			{
				for (int m = 0; m < N; m++)
				{
					vx = vx + conj(Eigenvector[N*i + n]) * dx[n][m] * Eigenvector[N*j + m];
					vy = vy + conj(Eigenvector[N*j + n]) * dy[n][m] * Eigenvector[N*i + m];
				}
			}
			vx = vx / double(N*N);
			vy = vy / double(N*N);
			Nb[j] = 1 / (exp((E[j] - mu) / kT) - 1);
			if (E[i] != E[j])
			{
				temp[1] = temp[1] - ((Nb[i] - Nb[j]) / (E[j] - E[i]))*(imag(vx*vy / ((omega + E[i] - E[j]) + delta)) + imag(vx*vy / ((omega + E[j] - E[i]) + delta)));
				temp[0] = temp[0] - ((Nb[i] - Nb[j]) / (E[j] - E[i]))*(imag(vx*conj(vx) / ((omega + E[i] - E[j]) + delta)) + imag(vx*conj(vx) / ((omega + E[j] - E[i]) + delta)));
			}
		}
	}
	delete[]Nb;
}
#endif // !__OPTCOND_H__
