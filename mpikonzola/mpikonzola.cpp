// mpikonzola.cpp : Defines the entry point for the console application.
//

#include <mpi.h>
#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#define n 100
int main(int argc, char** argv) {

	MPI_Init(&argc, &argv);

	int nprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
/*
#ifdef _DEBUG
	if (myrank == 0)
	{
		system("pause");
	}
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	
	char hostname[MPI_MAX_PROCESSOR_NAME];
	int name_len;
	MPI_Get_processor_name(hostname, &name_len);
	*/


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

	//GATHER int
	/*
	int isend, irecv[32];
	isend = myrank + 1;
	isend = myrank * myrank;
	MPI_Gather(&isend, 1, MPI_INT, irecv, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (myrank == 0) {
		for (int i = 0; i < nprocs * 2; i++)std::cout << irecv[i] << " ";
		std::cout << std::endl;
	}

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

	//ALLGATHER int[]
	/*
	isenda[0] = myrank + 1;
	isenda[1] = myrank * myrank;
	MPI_Allgather(isenda, 2, MPI_INT, irecv, 2, MPI_INT, MPI_COMM_WORLD);
		for (int i = 0; i < nprocs * 2; i++)std::cout << irecv[i] << " ";
		std::cout << std::endl;
	

	MPI_Barrier(MPI_COMM_WORLD);
	if (myrank == 0)std::cout << std::endl;
	*/

	//GATHERV (varied gather) int[] 3x proc
	/*
	int isenda[3],irecv[32];
	int sendc; //send count
	int recvc[3] = {1,2,3};//received counts from senders
	int disp[3] = {7,0,3}; //displace received
	for (int i = 0; i < myrank; i++) {
		isenda[i] = myrank + 1;
	}
	sendc = myrank + 1;
	MPI_Gatherv(isenda, sendc, MPI_INT, irecv, recvc, disp, MPI_INT, 0, MPI_COMM_WORLD);
	if (myrank == 0) {
		for (int i = 0; i < 8; i++)std::cout << irecv[i] << " ";
		std::cout << std::endl;
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	if (myrank == 0)std::cout << std::endl;
	*/

	//GATHERV (varied gather) int[] #2 4x proc
	/*
	int isenda[4], irecv[32];
	for (int i = 0; i < 32; i++)irecv[i] = -1;//mark empty
	
	int sendc; //send count
	int recvc[4] = { 3,2,6,4 };//received counts from senders
	int disp[4] = { 14,5,7,1 }; //displace received
	
	sendc = recvc[myrank];
	for (int i = 0; i < sendc; i++) {
		isenda[i] = myrank;
	}
	
	MPI_Gatherv(isenda, sendc, MPI_INT, irecv, recvc, disp, MPI_INT, 0, MPI_COMM_WORLD);
	if (myrank == 0) {
		for (int i = 0; i < 17; i++)std::cout << irecv[i] << " ";
		std::cout << std::endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if (myrank == 0)std::cout << std::endl;
	*/

	//REDUCE SUM 3x proc
/*
int start, end, sum, tp;
int a[9] = { 1,2,3,4,5,6,7,8,9 };

start = myrank * 3;
end = start + 2;
sum = 0;
for (int i = start; i <= end; i++) {
	sum +=a[i];
}

std::cout << "proc " << myrank << " subsum = " << sum << std::endl;
MPI_Reduce(&sum, &tp, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
sum = tp;
if (myrank == 0)std::cout << "totalsum = " << sum << std::endl;
*/

//REDUCE SUM 3x proc #2
/*
int start, end, sum, tp;
int a[16];
int range[3] = { 5,5,6 };
start = myrank * 5;
end = start + range[myrank];
sum = 0;
for (int i = start; i <= end; i++) {
	a[i] = i;
	sum += a[i];
}

std::cout << "proc " << myrank << " "<< range[myrank] << " elems; " << " subsum = " << sum << std::endl;
MPI_Reduce(&sum, &tp, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
sum = tp;
if (myrank == 0)std::cout << "totalsum = " << sum << std::endl;
*/

//REDUCE MAXLOC 3x proc
/*
int isend[2], irecv[2];
int start, end, max, loc,tp;
int a[9] = {2,5,8,4,25,48,5,4,1};
start = myrank * 3;
end = start + 2;
max = 0;
for (int i = start; i <= end; i++) {
	if (a[i] > max) {
		max = a[i];
		loc = i;
	}
}
isend[0] = max;
isend[1] = loc;
std::cout << "proc " << myrank << " " << " max = " << max << " loc = " << loc << std::endl;
MPI_Reduce(isend, irecv, 1, MPI_2INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
if (myrank == 0) {
	max = irecv[0];
	loc = irecv[1];
	std::cout << "totalmax = " << max << " loc = "<< loc<< std::endl;
}
*/
//REDUCE MAXLOC 3x proc #2

int isend[4], irecv[4];
int start, end, sum, max, min, locx, locn, tp;
double avg, sig;
int a[n];
float M[n][n];
sum = 0;
sig = 0;
max = -100000;
min = 100000;

if (myrank == 0) {
//generate
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			M[i][j] = (float)rand() / (RAND_MAX + 1) * 10.0 - 10.0;
		}
	}

for (int i = 0; i < n; i++) {
	a[i] = rand() % 100;
}

	//serial




	for (int i = 0; i < n; i++) {
		sum += a[i];
		if (a[i] > max) { max = a[i]; locx = i; }
		if (a[i] < min) { min = a[i]; locn = i; }
	}
	avg = (double)sum / n;
	std::cout << " *sum = " << sum << std::endl;
	std::cout << std::setprecision(20) << " *avg = " << avg << std::endl;

	for (int i = 0; i < n; i++) {
		sig += (double)((a[i] - avg)*(a[i] - avg));
	}
	sig = std::sqrt(sig / n);
	std::cout << std::setprecision(20) << " *stddev = " << sig << std::endl;
}

	//parallel
	sum = 0;
	sig = 0;
	max = -100000;
	min = 100000;

	
	int nlocal = n / nprocs;
start = myrank * nlocal;
end = start + nlocal -1;
if (myrank==nprocs-1)end = n-1;

for (int i = start; i <= end; i++) {

		//a[i] = i;
		sum += a[i];
		if (a[i] > max) { max = a[i]; locx = i; }
		if (a[i] < min) { min = a[i]; locn = i; }
	
}

avg = (double)sum / nlocal;
isend[0] = max;
isend[1] = locx;
isend[2] = min;
isend[3] = locn;
isend[4] = sum;

std::cout << "proc " << myrank << " (" << start << " - " << end << ") \n" 
<< " max = " << max << " [" << locx << "] \n" 
<< " min = " << min << " [" << locn << "] \n"
<< " sum = " << sum << "\n"
<< " avg = " << std::setprecision(20) << avg << std::endl;

MPI_Reduce(isend, irecv, 1, MPI_2INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
if (myrank == 0) {
	max = irecv[0]; locx = irecv[1];
	std::cout << "totalmax = " << max << " [" << locx << "]" << std::endl;
}
MPI_Reduce(isend+2, irecv, 1, MPI_2INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
if (myrank == 0) {
	min = irecv[0]; locn = irecv[1];
	std::cout << "totalmin = " << min << " [" << locn << "]" << std::endl;
}
MPI_Reduce(isend + 4, irecv, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
if (myrank == 0) {
	sum = irecv[0];
	std::cout << "totalsum = " << sum << std::endl;
	avg = (double)sum / n;
	std::cout << std::setprecision(20) << "totalavg = " << avg << std::endl;
}

MPI_Bcast(&avg, 1, MPI_INT, 0, MPI_COMM_WORLD);
sig = 0;
for (int i = start; i <= end; i++) {
	sig += (a[i] - avg)*(a[i] - avg);
}

isend[0] = sig;
MPI_Reduce(&sig,irecv,1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
if (myrank == 0) {
	sig = std::sqrt((double)irecv[0] / n);
}
std::cout << std::setprecision(20) << " par stddev = "<< sig << std::endl;


//end

	//MPI_Finalize();
	getchar();
	return 0;
}