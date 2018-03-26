#ifndef lapack_complex_double 
#define lapack_complex_double MKL_Complex16
#endif
#include <iostream>
#include <cstdlib>
#include<complex>
#define MKL_Complex16 std::complex<double>
#include "mkl.h"
#define N 2
#define mag 0.0
#include<fstream> 
#define pi 3.14159265358979323846
using namespace std;
int main()
{
		FILE *out;
		out=fopen("out.txt","w");
/*	ofstream out;
	out.open("out.txt",ios::app);
	if(!out.is_open())  
    {  
        cout<<"Could Not Open File!!!";  
        exit(EXIT_FAILURE);  
    }*/
	double kx,ky;
	double i1,j1;
	int i,j;
	int n , m;
	int jp1,jm1;
	double n1,m1;
	for (n=0;n<101;n++)
	{
		for (m=0;m<101;m++)
		{
			n1=double(n);
			m1=double(m);
		    kx = (-pi)+2.0*pi* n1/100.0;
		    //kx = -pi+2.0*pi* n1/100.0;
			ky = (-pi)+2.0*pi* m1/100.0;
		    //cout << kx << " " << ky <<endl;
            complex<double> a[N][N]={{0,0}};


/*------------------------------------------------- create the matrix------------------------------------------*/        
		
		for(i=0;i<N;i++)
        {
        	j=i;
        	jm1=j-1;
			jp1=j+1;
        	i1=double(i);
			a[i][j] ={2.0*cos(ky-mag*2.0*pi*i1/double(N)),0.0};
			//cout<<2.0*cos(ky-2.0*pi*i1/double(N))<<a[i][j];
			if (i==N-1)
			{
				jp1=0;
			}
				
			if (i==0)
			{
				jm1=N-1;
			}
        	a[i][jp1] = {cos(kx),-sin(kx)};
        	a[i][jm1] = {cos(kx),sin(kx)};	
		}
/*------------------------------------------------ create the matrix over----------------------------------------*/
//next is to diagonalize the matrix a[][]
//----------------------------------------------------------------------------------------------------------------
//

/*for (i=0;i<N;i++)
{
		for (j=0;j<N;j++)
		{
				cout << a[i][j] << " ";
		}
		cout << endl;
		}
getchar();*/
		
		char jobz='N' ;
        char uplo='U' ;
        //float wr[N];
        //float wi[N];

		double w[N];
        //lapack_complex_double vl[N];
        //lapack_complex_double vr[N];
        lapack_int lda = N;
        //lapack_int ldvl = N;
        //lapack_int ldvr = N;

        LAPACKE_zheev(LAPACK_ROW_MAJOR, jobz, uplo, N, *a, lda, w); //diagonalization
/*for (i=0;i<N;i++)
{
		for (j=0;j<N;j++)
		{
				cout << real(a[i][j]) << " ";
		}
		cout << endl;
		}
getchar();*/
		fprintf(out,"%f %f ",kx,ky);
		for (i=0;i<N;i++)
		{
				//fprintf(out, "%f %f ",kx,ky);
				fprintf(out, " %f ",w[i]);
				//cout << w[i] << " ";
		}
		fprintf(out, "\n");
		//cout << endl;

	}

	}
		fclose(out);
}
