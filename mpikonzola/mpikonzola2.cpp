// mpikonzola.cpp : Defines the entry point for the console application.
//

#include <mpi.h>
#include "stdafx.h"
#include <chrono>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#define n 150

int main(int argc, char** argv) {

	MPI_Init(&argc, &argv);

	int nprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

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


int start, end, tp;
std::chrono::time_point<std::chrono::system_clock> tstart, tend;
std::chrono::duration<double> elapsed_seconds;
int a[n];
float s[n];
float M[n][n];
float **Mloc;
float r[n];

if (myrank == 0) {

	//generate
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			M[i][j] = (float)rand() / (RAND_MAX + 1) * 20.0 - 10.0;
		}
	}

	for (int i = 0; i < n; i++) {
		r[i] = 0;
		s[i] = 0;
		a[i] = rand() % 100;
	}

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
	std::cout << std::endl << std::endl;


	//serial
	tstart = std::chrono::system_clock::now();//clockstart
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			s[i] += (float)a[j] * M[i][j];
			s[i] += pow((float)a[j] * M[i][j], 5); //ubersucin
		}
	}
	tend = std::chrono::system_clock::now();//clockstop
	elapsed_seconds = tend - tstart;

	for (int i = 0; i < n; i++)std::cout << s[i] << " ";
	std::cout << std::endl << std::endl;

	std::cout << elapsed_seconds.count() << " s\n"; 
	std::cout << std::endl << std::endl;

}

	//bcast
MPI_Bcast(r, n, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast(a, n, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast(M, n*n, MPI_FLOAT, 0, MPI_COMM_WORLD);

	//parallel

tstart = std::chrono::system_clock::now();//clockstart	
		int nlocal = n / nprocs +1;
		int nlast;
		start = myrank * nlocal;
		end = start + nlocal -1;
		if (myrank == nprocs - 1) {
			end = n - 1;
			nlast = end - start + 1;
		}

for (int i = start; i <= end; i++) {
	for (int j = 0; j < n; j++) {
		r[i] += (float)a[j] * M[i][j];
		r[i] += pow((float)a[j] * M[i][j],5); //ubersucin
	}}

float irecv[n];
MPI_Gather(r+start, nlocal, MPI_FLOAT, irecv, nlocal, MPI_FLOAT, 0, MPI_COMM_WORLD);

tend = std::chrono::system_clock::now();//clockstop
elapsed_seconds = tend - tstart;
if (myrank == 0) {
	for (int i = 0; i < n; i++)std::cout << irecv[i] << " ";
	std::cout << std::endl;
	std::cout << elapsed_seconds.count() << " s\n";
}

//end

	MPI_Finalize();
	getchar();
	return 0;
}