#pragma once
#include<iostream>
#include<complex>
#include<cmath>
#include<cstring>
#define pi 3.1415926535898 
using namespace std;
void createm(int N, complex <double> **a, double mag, double kx, double ky, double t1, double t2)
{
	//std::complex<double> **a = new std::complex<double>* [4 * N];
	int i, j, im1, jm1, jp1, ip1, ip2, jp2, ip3, jp3;
	double theta, theta1, theta2, theta3;
	//memset(a, 0.0, sizeof(a));
	for (i = 0; i < 4 * N; i = i + 2)
	{
		j = i;
		ip1 = i + 1;
		ip2 = i + 2;
		ip3 = i + 3;
		jp1 = j + 1;
		jp2 = j + 2;
		jp3 = j + 3;
		theta = mag*2.0*pi*(double(i / 2.0)) / double(N);
		theta1 = mag*2.0*pi*(0.125 + double(i / 2.0) / 2.0) / double(N);
		theta2 = mag*2.0*pi*(0.375 + double(i / 2.0) / 2.0) / double(N);
		theta3 = mag*2.0*pi*(0.375 + double(N - 1) / 2.0 + double(i / 2.0) / 2.0) / double(N);
		
		//for i
		a[jp1][i] = { t1*cos(theta1) + t1*cos(ky + theta1),-t1*sin(theta1) + t1*sin(ky + theta1) };
		if (i < 4 * N - 2)
		{
			a[j][ip2] = { t2,0 };
			a[jp2][i] = { t2,0 };
			a[jp3][i] = { t1*cos(theta2) + t1*cos(ky + theta2),-t1*sin(theta2) + t1*sin(ky + theta2) };
		}
		
		if (i == 0)
		{
			a[4 * N - 2][i] = { t2*cos(kx), t2*sin(kx) };
		}
		if (i == 4 * N - 2)
		{
			a[0][i] = { t2*cos(kx)+real(a[0][i]),-t2*sin(kx)+imag(a[0][i]) };
			a[1][i] = { t1*cos(kx + theta3) + t1*cos(theta3 + ky - kx)+real(a[1][i]),-t1*sin(kx + theta3) + t1*sin(theta3 + ky - kx)+imag(a[1][i]) };
		}
		
		//for ip1
		a[j][ip1] = { t1*cos(theta1) + t1*cos(ky + theta1)+real(a[j][ip1]),t1*sin(theta) - t1*sin(ky + theta1)+imag(a[j][ip1]) };
		a[ip1][ip1] = { 2 * cos(ky + theta),0 };
		if (i < 4 * N - 2)
		{
			a[j][ip3] = { t1*cos(theta2) + t1*cos(ky + theta2),t1*sin(theta2) - t1*sin(ky + theta2) };
		}

		if (i == 0)
		{
			a[4 * N - 2][ip1] ={ t1*cos(kx - theta1) + t1*cos(ky - kx - theta1)+real(a[4 * N - 2][ip1]),\
			t1*sin(kx - theta1) - t1*sin(ky - kx - theta1)+imag(a[4 * N - 2][ip1]) };
		}
	}
}
