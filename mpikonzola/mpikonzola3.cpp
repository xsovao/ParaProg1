// mpikonzola.cpp : Defines the entry point for the console application.
//

#include <mpi.h>
#include "stdafx.h"
#include <chrono>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#define M_PI         3.141592653589793238462643383279502884

double R = 3671, D = 800, GM = 398600;

double ss(double *a,double *b) {
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

double ss(double *a, double *b, int n) {
	double rr = 0;
	for (int i = 0; i < n; i++)rr += a[i] * b[i];
	return rr;
}

double* d(double *a, double *b) {
	//std::cout << a[0] << "," << a[1] << "," << a[2] << std::endl;
	//std::cout << b[0] << "," << b[1] << "," << b[2] << std::endl;
	double f[3];
	for(int i=0; i<3; i++){
		f[i] = a[i] - b[i];}
	return f;
}

double s(double *a) { return ss(a, a); }
double dtr(double a) { return a /180.0 * M_PI; }


void map(double *X[],double R, double B, double L) {
	*X=new double[3];
	(*X)[0] = R * cos(dtr(B)) * cos(dtr(L));
	(*X)[1] = R * cos(dtr(B)) * sin(dtr(L));
	(*X)[2] = R * sin(dtr(B));
}

double norm(double *a,int n) {
	double r = 0;
	for (int i = 0; i < n; i++)r += a[i] * a[i];
	return sqrt(r);
}
double *scale(double a,double *x) {
	double d[3];
	d[0] = a*x[0]; d[1] = a*x[1]; d[2] = a*x[2];
	return d;
}
/*
double *map(double R, double B, double L) {
	double X[3];
	X[0] = R * cos(dtr(B)) * cos(dtr(L));
	X[1] = R * cos(dtr(B)) * sin(dtr(L));
	X[2] = R * sin(dtr(B));
	return X;
}*/

void makeM(double **A, int n, double **X) {
	double r[3], size,nn[3], sc = (R-D)/R;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			r[0] = X[i][0] - sc*X[j][0];
			r[1] = X[i][1] - sc*X[j][1];
			r[2] = X[i][2] - sc*X[j][2];
			nn[0] = X[j][0] / R;
			nn[1] = X[j][1] / R;
			nn[2] = X[j][2] / R;
			size = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
			//std::cout << r[0] << " " << r[1] << " " << r[2] << std::endl;
			//if(i==j)std::cout << size << std::endl;
			A[i][j] = 1.0 / (4 * M_PI * size * size * size) * ss(r, nn);
		}
	}
}

void printM(double **A,int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < 1; j++) {
			std::cout << A[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

int main(int argc, char** argv) {

	MPI_Init(&argc, &argv);

	int nprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	int start, end, tp;
	std::chrono::time_point<std::chrono::system_clock> tstart, tend;
	std::chrono::duration<double> elapsed_seconds;
	int n = 902;
	int maxit = 100;
	double **M;

	double *Q1,*Q2;
	double **X;

	double *x,q,tol=1e-9;

	double rh,rhold, bet, aph, omg,nr,ss1,ss2;
	double *r, *rold, *p, *s,*v,*t;

	Q1 = new double[n];
	Q2 = new double[n];
	X = new double*[n];
	M = new double*[n];
	x = new double[n];
	r = new double[n];
	rold = new double[n];
	p = new double[n];
	s = new double[n];
	v = new double[n];
	t = new double[n];

	q = GM/(R*R);

	if (myrank == 0) {
		std::ifstream file;
		file.open("D:\\sova\\PAR\\MPI_Z1\\BL-902.dat", std::ios::in);
		if (file.is_open()) {
			for (int i = 0; i < n; i++) {

				X[i] = new double[3];
				M[i] = new double[n];
				
				file >> X[i][0] >> X[i][1] >> X[i][2] >> Q1[i] >> Q2[i];

				map(&(X[i]), R, X[i][0], X[i][1]);
				//std::cout <<std::setprecision(6) << X[i][0] << "\t" << X[i][1] << "\t" << X[i][2] << "\t" << Q1[i] << "\t" << Q2[i] << std::endl;
			}
		}
		else {
			std::cout << "not open";
		}

		makeM(M, n, X);
		//std::cout << 1 / (4 * M_PI * D *D) << "\n";;
		std::cout << "dun matica" << "\n";
		//printM(M, n);
		
		std::cout << "____" << std::endl;
		for (int i = 0; i < n; i++) {
			x[i] = 0;
			r[i] = q;
			rold[i] = q;
		}
		for (int k = 0; k < maxit; k++) {

			std::cout << "it " << k << std::endl;
			
			ss1 = 0;
			for (int i = 0; i < n; i++) {
				ss1 += rold[i] * r[i];
			}
			rh = ss1;
			if (rh == 0) {
				std::cout << "fail" << std::endl;
				break;
			}
			if (k == 0)
				for (int i = 0; i < n; i++)p[i] = r[i];
			else {
				bet = (rh / rhold)*(aph / omg);
				for (int i = 0; i < n; i++) {
					p[i] = r[i] + bet*(p[i] - omg*v[i]);
				}
			}
			for (int i = 0; i < n; i++) {
				v[i] = 0;
				for (int j = 0; j < n; j++)v[i] += M[i][j] * p[j];
			}


			ss1 = 0;
			for (int i = 0; i < n; i++) {
				ss1 += rold[i] * v[i];
			}
			aph = rh / ss1;

			for (int j = 0; j < n; j++)s[j] = r[j] - aph * v[j];
			if (norm(s, n) < tol) {
				for (int j = 0; j < n; j++)x[j] = x[j] + aph*p[j];
				break;
			}
			for (int i = 0; i < n; i++) {
				t[i] = 0;
				for (int j = 0; j < n; j++)t[i] += M[i][j] * s[j];
			}

			ss1 = 0;
			ss2 = 0;
			for (int i = 0; i < n; i++) { 
				ss1 += t[i] * s[i];
				ss2 += t[i] * t[i];
			}
			omg = ss1/ss2;
			for (int j = 0; j < n; j++)x[j] = x[j] + aph*p[j] + omg * s[j];
			for (int j = 0; j < n; j++)r[j] = s[j] - omg * t[j];
			
			rhold = rh;
			nr = norm(r, n);
			std::cout << "norm[r] = "<< nr << "\n";
			if (nr < tol)break;
		}



	
	}
	//BROADCAST 
	/*
	int a[4];
	for (int i = 0; i < 4; i++) {
		if (myrank == 0)a[i] = i;
		else a[i] = 0;}

	std::cout << "proc " << myrank << " original a = " << a[0] << "," <<a[1] << "," << a[2] << "," << a[3] << std::endl;
	MPI_Bcast(a, 4, MPI_INT, 0, MPI_COMM_WORLD);;
	std::cout << "proc " << myrank << " received a = " << a[0] << "," << a[1] << "," << a[2] << "," << a[3] << std::endl;
	
	MPI_Barrier(MPI_COMM_WORLD);
	if (myrank == 0)std::cout << std::endl;
	*/
	//GATHER int[]
	/*
	int isenda[2];
	isenda[0] = myrank + 1;
	isenda[1] = myrank * myrank;
	MPI_Gather(isenda, 2, MPI_INT, irecv, 2, MPI_INT, 0, MPI_COMM_WORLD);
	if (myrank == 0) {
		for (int i = 0; i < nprocs*2; i++)std::cout << irecv[i] << " ";
		std::cout << std::endl;
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	if (myrank == 0)std::cout << std::endl;
	*/

if (myrank == 0) {

	//write
	/*
	std::cout << "M = {" << std::endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			//std::cout<< M[i][j] << " ";
		}
		std::cout << std::endl;
	}*/
	for (int i = 0; i < n; i++) {
	//	std::cout << a[i] << " ";
	}
	//std::cout << std::endl << std::endl;

	
	//serial
	tstart = std::chrono::system_clock::now();//clockstart
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
		}
	}
	tend = std::chrono::system_clock::now();//clockstop
	elapsed_seconds = tend - tstart;

	//for (int i = 0; i < n; i++)std::cout << s[i] << " ";
	std::cout << std::endl << std::endl;

	std::cout << elapsed_seconds.count() << " s\n"; 
	std::cout << std::endl << std::endl;

}

	//bcast
/*
MPI_Bcast(a, n, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast(M, n*n, MPI_FLOAT, 0, MPI_COMM_WORLD);
*/
	//parallel
/*
tstart = std::chrono::system_clock::now();//clockstart	
		int nlocal = n / nprocs +1;
		int nlast;
		start = myrank * nlocal;
		end = start + nlocal -1;
		if (myrank == nprocs - 1) {
			end = n - 1;
			nlocal = end - start + 1;
		}

		std::cout << start << "-" << end << " nloc:" << nlocal << std::endl;
		Mloc = new float *[nlocal];
		r = new float[nlocal];
		for (int i = 0; i < nlocal; i++) {
			Mloc[i] = new float[n];
			r[i] = 0;
		}

			for (int i = 0; i < nlocal; i++) {
			for (int j = 0; j < n; j++) {
				Mloc[i][j] = M[i+start][j];
				//std::cout << Mloc[i][j] << " ";
			}
			//std::cout << std::endl;

		}

		for (int i = 0; i < nlocal; i++) {
			for (int j = 0; j < n; j++) {				
		r[i] += (float)a[j] * Mloc[i][j];
		r[i] += pow((float)a[j] * Mloc[i][j],5); //ubersucin
	}}

float irecv[n];
MPI_Gather(r, nlocal, MPI_FLOAT, irecv, nlocal, MPI_FLOAT, 0, MPI_COMM_WORLD);

tend = std::chrono::system_clock::now();//clockstop
elapsed_seconds = tend - tstart;
if (myrank == 0) {
	for (int i = 0; i < n; i++)std::cout << irecv[i] << " ";
	std::cout << std::endl;
	std::cout << elapsed_seconds.count() << " s\n";
}
*/
//end

	MPI_Finalize();
	getchar();
	return 0;
}