// mpikonzola.cpp : Defines the entry point for the console application.
//

#include <mpi.h>
#include "stdafx.h"
#include <chrono>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
# define M_PI         3.141592653589793238462643383279502884L

double ss(double *a,double *b) {
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

double* d(double *a, double *b) {
	double f[3];
	for(int i=0; i<3; i++){
		f[i] = a[i] - b[i];}
	return f;
}

double s(double *a) { return ss(a, a); }
double dtr(double a) { return a / 180.0 * M_PI; }

double* n(double R, double B, double L) {
	double x[3];
	x[0] = std::cos(dtr(B)) * std::cos(dtr(L));
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
	float **M;
	float **Mloc;
	float *r;
	char *buf;
	double R = 3671,D=800;
	double *Q1,*Q2;
	double **X;

	Q1 = new double[n];
	Q2 = new double[n];
	X = new double*[n];
	std::ifstream file;
	file.open("D:\\sova\\PAR\\MPI_Z1\\BL-902.dat", std::ios::in);
	if (file.is_open()) {
		for (int i = 0; i < n; i++) {
			X[i] = new double[3];
			file >> X[i][0] >> X[i][1] >> X[i][2] >> Q1[i] >> Q2[i];
			std::cout << X[i][0] << " " << X[i][1] << " " << X[i][2] << " " << Q1[i] << " " << Q2[i] << std::endl;
		}
	}
	else {
		std::cout << "not open";
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
	
		
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